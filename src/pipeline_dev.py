#!/usr/bin/env python
# Consensus error correction for nanopore sequencing
# bilgenurb@gmail.com

# Loading libraries
import numpy as np
import pandas as pd

# For interfacing with the file system
import glob
import subprocess
import os
import time
import sys

# Loading the custom libraries
import bilge_pype as bpy

#########################################################################################################

def generate_pseudo_refdb(primers, reads, block = 10000, f_sizes = '[[500,1000]]',
                          workspace = './pseudo_refdb/', config = '-k5 -w1 --score-N 0 -s 20 -P',
                          log_file = 'log.txt'):
    '''
    Build reference database using paired end search of some primer sequences in the uncorrected reads
    primers   = dataframe containing primer information with columns:
                [fwd_id, fwd_seq, rev_id, rev_seq]
    reads     = dataframe containing read information with columns:
                [id, sequence] or [id, sequence, quality]
    
    block     = block size of reads to work on each time. This is for saving on RAM.

    workspace = workspace for the aligner
    config    = configs to pass to minimap2 aligner
    pseudodb_file = csv file to save pseudo reference database
    log_file  = where to log progress
    
    returns location of pseudodb_file
    '''
    def parse_fsize(text):
        text = text.split('[[')[1].split(']]')[0].split('],[')
        data = []
        for i in text:
            data.append(i.split(','))
        return np.array(data).astype(int)

    def split_sequence(df):
        '''
        Function to get sequence info for each frag
        df = pandas dataframe with columns:
             [id, orientation, start, end, sequence, quality, fs_lower, fs_upper]
        '''
        data = []
        for r_id, ot, t1, t2, seq, qual in df.values:
            seq = seq[t1:t2]
            qual = qual[t1:t2]
            if ot == '-':
                seq = bpy.dna_revcomp(seq)
                qual = qual[::-1]
            data.append([r_id+'_rs'+str(t1)+'_'+str(t2), seq, qual])
        return pd.DataFrame(data, columns = ['id', 'sequence', 'quality'])

    # look for the primers
    data = []
    for s1,s2 in parse_fsize(f_sizes):
        # filter for some read sizes to reduce search time
        s = (reads['length'] > s1)&(reads['length'] < s2)
        # chunk up the query so we can see progress
        N_chunks = np.ceil(len(reads[s])/block)
        for i, df_i in enumerate(np.array_split(reads[s], N_chunks)):
            bpy.append_to_log('Working on '+str(i)+'/'+str(N_chunks)+' block size = '+str(block)+
                              ' read size = '+str(s1)+' to '+str(s2)+' bp', log_file)
            data.append(match_primer_to_read(primers, df_i, workspace = workspace, log_file = log_file))

    if len(data)==0:
        bpy.append_to_log('generate_pseudo_refdb: no primers found in reads. Exiting.', log_file)
        sys.exit(1)
    data = pd.concat(data)
    
    # get only reads with both fwd and rev sequences
    data = data[(data['fwd_primer']!='') & (data['rev_primer']!='')]
    
    # set the correct start and end
    x = []
    for i in data[['id','fwd_orientation','fwd_q_start','fwd_q_end','rev_q_start','rev_q_end']].values:
        if i[1] == '+':
            x.append([i[0], i[1], i[2], i[5]])
        else:
            x.append([i[0], i[1], i[4], i[3]])
    data = pd.DataFrame(x, columns = ['id','orientation','start','end'])

    # merge and get pseudo reference sequences
    data = data.merge(reads[['id','sequence','quality']], on = 'id', how = 'left')
    df_d = split_sequence(data[['id','orientation','start','end','sequence','quality']]).drop_duplicates()
    bpy.append_to_log('Reference sequences found '+str(len(df_d)), log_file)
    return df_d
    
def find_RCA_frags(ref, reads, block = 10000, folder = './frags/', config = '-k10 -w1', log_file = 'log.txt'):
    '''
    Function to search sequencing primers and identify RCA reads from fastq info
    ref            = similarily related sequences to search repeated RCA frags with columns = [<id>,<sequence>]
    reads          = fastq info
    block          = number of fastq reads to process each time
    folder         = where to save data
    '''
    aligner_folder = folder+'aligner/'

    N_chunks = np.ceil(len(reads)/block)
    build_index = True
    for i, reads in enumerate(np.array_split(reads, N_chunks)):
        bpy.append_to_log('Working on '+str(i)+'/'+str(N_chunks)+' block size = '+str(block), log_file)
        # Find frag via minimap2
        frags = bpy.run_minimap2(reads[['id','sequence','quality']], ref,
                                 workspace = aligner_folder, config = config, cigar = True,
                                 build_index = build_index, use_index = True)
        build_index = False
        # Save the info
        bpy.append_to_log('Saving information', log_file)
        frags.to_csv(folder+'fraginfo_'+str(i)+'.csv.gz', index=False, compression='infer')
        bpy.append_to_log('frags found in this block = '+str(len(frags)), log_file)
        
    # clear the work folder
    subprocess.run(['rm','-r',aligner_folder]) # remove the aligner folder after work is done

def annotate_reads(df, log_file = 'log.txt'):
    '''
    Annotates some information about each read
    
    df = pandas dataframe with the following columns:
    [id, q_len, q_start, q_end, t_len, t_start, t_end, orientation, AS, match_score]
    
    return df_frags and df_1D2

    df_frags has the following columns: 
    [id, N frags, 1D2, q_len, q_unmapped, db_fwd_cover, db_rev_cover, avg_match, avg_AS, std_AS]
    query_id     = read id
    N frags      = number of repeated sequences from rca library prep
    1D2          = if read orientations are inverting, count how many times they invert
    q_len        = read length
    q_unmapped   = how many nucleotides did not map to anything in pseudo reference database
    db_fwd_cover = space between 5' of db and 5' of mapped  or t_start
    db_rev_cover = space between 3' of db and 3' of mapped  or t_len - t_end
    match_size   = nucleotides of match to pseudo reference or t_end - t_start
    avg_match    = average match score of pseudo reference to read
    avg_AS       = average alignment score AS
    std_AS       = standard deviation of alignment score AS
    
    df_1D has the following columns:
    [id, f_start, f_end]
    id           = read id + 1D2num
    f_start      = where the frag starts in the read
    f_end        = where the frag ends in the read
    
    '''
    df = df.sort_values(by = ['query_id','q_start']) # check if sorting on numpy array will be faster
    data = df[['query_id','q_len','q_start','q_end',
               't_len','t_start','t_end','orientation','AS','match_score']].values
    
    # storage variables
    df_1D = []
    df_frag = []
    # initial variables
    N_frags = 0
    N_1D2 = 0
    f_start = 0
    q_unmapped = 0
    db_fwd_cover = 0
    db_rev_cover = 0
    match_size = 0
    match_score = []
    AS = []
    
    # iterate on each read
    for i in range(0,len(data)):
        # count frags present
        N_frags+=1

        # 1D2 read inversion found
        if i+1 < len(data) and data[i,7]!=data[i+1,7] and data[i,0]==data[i+1,0]:
            df_1D.append([data[i,0]+'_1D2num'+str(N_1D2),f_start,data[i,3]])
            N_1D2+=1
            f_start = data[i+1,2]

        # q_unmapped computation by first summing mapped sequence space
        q_unmapped+= data[i,3]-data[i,2]
        # sum db coverage, match size, AS
        db_fwd_cover+= data[i,5]
        db_rev_cover+= data[i,4]-data[i,6]
        match_size+= data[i,6]-data[i,5]
        match_score.append(data[i,9])
        AS.append(data[i,8])
        
        # when end of read reached
        if i+1 >= len(data) or data[i,0]!=data[i+1,0]:
            # compute q_unmapped
            q_unmapped = data[i,1]-q_unmapped
            
            # record data for df_frag
            df_frag.append([data[i,0],N_frags,N_1D2,data[i,1],q_unmapped,
                            db_fwd_cover,db_rev_cover,
                            match_size,np.mean(match_score),np.mean(AS),np.std(AS)])
            
            # record data for df_1D2
            df_1D.append([data[i,0]+'_1D2num'+str(N_1D2),f_start,data[i,1]])
            
            # reset initial variables
            N_frags = 0
            N_1D2 = 0
            f_start = 0
            q_unmapped = 0
            db_fwd_cover = 0
            db_rev_cover = 0
            match_size = 0
            match_score = []
            AS = []

    # format data
    df_frag = pd.DataFrame(df_frag, columns = ['query_id','N frags','1D2','q_len','q_unmapped',
                                            'db_fwd_cover','db_rev_cover',
                                            'match_size','avg_match','avg_AS','std_AS'])
    df_1D = pd.DataFrame(df_1D, columns = ['query_id','f_start','f_end'])
    
    # save the info
    df_1D.to_csv('df_1D.csv.gz', index = False, compression = 'infer')
    df_frag.to_csv('df_frags.csv.gz', index = False, compression = 'infer')

    return df_frag

def perform_spoa_MSA(df, frags, batch_size=100, folder = './msa/',
                     thread_lock = True, config = '-l 0 -r 2', log_file = 'log.txt'):
    '''
    This function performs multisequence alignment on RCA reads
    df     = dataframe containing the original read fastq information
    frags  = RCA fragments with start and stop locations found by the previous aligner
    output = folder containing file of aligned RCA reads
    '''
    # make workspace folder if it does not exist
    bpy.check_dir(folder)
    # clear out the folder if anything is in it
    files = glob.glob(folder+'*')
    bpy.batch_file_remove(files)
    
    # merge the input data
    data = frags[['query_id','q_len','q_start','q_end','orientation']]
    data = data.rename(columns={'query_id':'id'})
    data = data.merge(df[['id','sequence','quality']], on = 'id', how = 'left')
    data = data.sort_values(by=['id','q_start']).values

    bpy.append_to_log('Writing frags to fastq files', log_file)
    MSA_infile = []
    seqlist = []
    k = 0
    for i in range(0,len(data)):
        # for first entries, set s1 = 0
        if i==0 or data[i,0]!=data[i-1,0]:
            s1 = 0
            s2 = int((data[i,3]+data[i+1,2])/2)
        # last entry reached, s2 to end of read
        elif i+1 >= len(data) or data[i,0]!=data[i+1,0]:
            s1 = int((data[i,2]+data[i-1,3])/2)
            s2 = data[i,1]
        else:
            s1 = int((data[i,2]+data[i-1,3])/2)
            s2 = int((data[i,3]+data[i+1,2])/2)
        
        # check that averaging did not create weird boundaries
        if s2 - s1 > 10:
            # slice out the sequences
            seq = data[i,5][s1:s2]
            q = data[i,6][s1:s2]

            # in case we have reverse sequences
            if data[i,4] == '-': 
                seq = bpy.dna_revcomp(seq)
                q = q[::-1]

            # record it
            seqlist.append([data[i,0]+'_frag'+str(k)+'_start'+str(s1)+'_end'+str(s2), seq, q])
        k+=1

        # end of list or end of read reached
        if i+1 >= len(data) or data[i,0]!=data[i+1,0]:
            if len(seqlist) > 1:
                fname = folder + data[i,0] + '_Nfrags'+str(len(seqlist))+'.fq'
                bpy.write_fastq(fname, np.array(seqlist))
                MSA_infile.append(fname)
            # reset the variables
            seqlist = []
            k = 0
        
        # submit to spoa every 100 files or whatever batch_size is
        if (len(MSA_infile) >= batch_size or i+1 >= len(data)) and len(MSA_infile) > 0:
            bpy.append_to_log('MSA progress = '+str(i)+'/'+str(len(data)), log_file)
            upy.run_msa(MSA_infile, aligner = 'spoa', config = config, thread_lock = thread_lock)
            MSA_infile = []
            
def get_spoa_consensus(MSA_outfile, log_file = 'log.txt'):
    '''
    Reads spoa outfile for the consensus from multi-sequence alignment
    MSA_outfile = multi-sequence alignment output from spoa

    output = consensus sequence information
    '''
    bpy.append_to_log('Reading consensus info from spoa output files', log_file)
    data = []
    for i in range(0,len(MSA_outfile)):
        fname = MSA_outfile[i]
        f = open(fname,'r')
        text = f.read().split('\n')
        f.close()
        rcaID = fname.split('/')[-1].split('_Nfrags')[0]
        nfrags = int(fname.split('/')[-1].split('_Nfrags')[1].split('.')[0])
        if len(text) < 2:
            bpy.append_to_log('file read error on '+fname, log_file)
        else:
            data.append([rcaID, fname.split('.out')[0], nfrags, text[1]])
        if i%20000 == 0:
            bpy.append_to_log('Consensus progress = '+str(i)+'/'+str(len(MSA_outfile)), log_file)
    return pd.DataFrame(data, columns=['id','msa_input_file','N frags','consensus'])

def match_primer_to_read(primers, reads, thresh = 10,
                         config = '-k5 -w1 --score-N 0 -s 20 -P', workspace = './match_pmr_to_rd/',
                         compact = True, verbose = True, log_file = 'log.txt'):
    '''
    Takes a list of primer pairs and reads and matches primer pairs to the reads
    primers   = pandas dataframe with columns = [fwd_id, fwd_sequence, rev_id, rev_sequence]
    reads     = pandas dataframe with columns = [id, sequence]
    config    = config passed to minimap2
    workspace = workspace used by the aligner
    output    = pandas dataframe with primer mapping info and sequence id as the first column
    '''
    # search for forward primers
    bpy.append_to_log('searching forward primers', log_file)
    database = primers.rename(columns={'fwd_id':'id','fwd_seq':'sequence'})
    df1 = bpy.run_minimap2(reads, database[['id','sequence']], workspace, config,
                           build_index = True, use_index = True)
    # search for reverse primers
    bpy.append_to_log('searching reverse primers', log_file)
    database = primers.rename(columns={'rev_id':'id','rev_seq':'sequence'})
    df2 = bpy.run_minimap2(reads, database[['id','sequence']], workspace, config,
                           build_index = True, use_index = True)
    
    # remove the working directory
    subprocess.run(['rm','-r',workspace])

    # merge the primers that were found
    if len(df1) == 0 and len(df2)==0:
        bpy.append_to_log('warning no fwd and rev primers found', log_file)
        return []
    elif len(df1) == 0:
        bpy.append_to_log('warning no fwd primers found', log_file)
        df2['pdir'] = 'fwd'
        df = df2
    elif len(df2) == 0:
        bpy.append_to_log('warning no rev primers found', log_file)
        df1['pdir'] = 'rev'
        df = df1
    else:
        bpy.append_to_log('fwd and rev primers found', log_file)
        df1['pdir'] = 'fwd'
        df2['pdir'] = 'rev'
        df = pd.concat([df1,df2]) # concatenate
  
    # finding best matches
    df = bpy.remove_overlaps(df, 'AS', thresh, log_file)
    
    bpy.append_to_log('finding best primer match for each read', log_file)
    # do a sort by query and position of the primers
    df = df.sort_values(by=['query_id','q_start']).reset_index(drop=True)
    data = df[['query_id','database_id','pdir','orientation',
               'q_start','q_end','AS']].values
    # find fwd primer and then rev primer
    pairs = []
    cur_pair = [-1,-1,-np.inf] # [index of fwd primer, index of rev primer, match score of combined]
    best_pair = cur_pair
    i = 0
    while i < len(data):
        # rev- only
        if data[i,2] == 'rev' and data[i,3] == '-':
            cur_pair[1] = i
            cur_pair[2] = data[i,6]
        # fwd- only
        elif data[i,2] == 'fwd' and data[i,3] == '-':
            cur_pair[0] = i
            cur_pair[2] = data[i,6]
        # fwd+ detected
        elif data[i,2] == 'fwd' and data[i,3] == '+':
            cur_pair[0] = i
            cur_pair[2] = data[i,6]
            # fwd+ and rev-
            if i+1 < len(data) and data[i,0] == data[i+1,0] and data[i+1,2] == 'rev' and data[i+1,3] == '-':
                cur_pair[1] = i+1
                cur_pair[2]+= data[i+1,6]
                i+=1
        # rev+ detected
        elif data[i,2] == 'rev' and data[i,3] == '+':
            cur_pair[1] = i
            cur_pair[2] = data[i,6]
            # rev+ and fwd-
            if i+1 < len(data) and data[i,0] == data[i+1,0] and data[i+1,2] == 'fwd' and data[i+1,3] == '-':
                cur_pair[0] = i+1
                cur_pair[2]+= data[i+1,6]
                i+=1
        
        # check for the best pair
        if best_pair[2] < cur_pair[2]:
            best_pair = cur_pair
        cur_pair = [-1,-1,-np.inf]
        
        # record final pick only when end of list or end of read reached
        if i+1 >= len(data) or data[i,0] != data[i+1,0]:
            pairs.append(best_pair)
            best_pair = [-1,-1,-np.inf]
        i+=1
    
    # split fwd and rev primers. then combine the info
    pairs = np.array(pairs)
    col1 = ['fwd_'+i for i in df.columns]
    col2 = ['rev_'+i for i in df.columns]
    
    tmp = []
    # fwd and reverse primers found
    s = pairs[(pairs[:,0] > -1)&(pairs[:,1] > -1)]
    df_f1 = []
    df_r1 = []
    if len(s) > 0:
        df_f1 = df.iloc[s[:,0]]
        df_r1 = df.iloc[s[:,1]]
        df_f1.columns = col1
        df_r1.columns = col2
        df_f1.loc[:,'id'] = df_f1['fwd_query_id'].values
        df_r1.loc[:,'id'] = df_r1['rev_query_id'].values
        tmp.append(df_f1.merge(df_r1, on = 'id', how = 'left'))
    
    # fwd primer only
    s = pairs[(pairs[:,0] > -1)&(pairs[:,1] == -1)]
    df_f2 = []
    if len(s) > 0:
        df_f2 = df.iloc[s[:,0]]
        df_f2.columns = col1
        df_f2.loc[:,'id'] = df_f2['fwd_query_id'].values
        # format the empty datatypes
        for i in range(0,len(col2)):
            if type(df_f2[col1[i]].values[0]) == np.int64:
                df_f2.loc[:,col2[i]] = 0
            elif type(df_f2[col1[i]].values[0]) == np.float64:
                df_f2.loc[:,col2[i]] = 0.0
            elif type(df_f2[col1[i]].values[0]) == str:
                df_f2.loc[:,col2[i]] = ''
        tmp.append(df_f2)
    
    # rev primer only
    s = pairs[(pairs[:,0] == -1)&(pairs[:,1] > -1)]
    df_r3 = []
    if len(s) > 0:
        df_r3 = df.iloc[s[:,1]]
        df_r3.columns = col2
        df_r3.loc[:,'id'] = df_r3['rev_query_id'].values
        # format the empty datatypes
        for i in range(0,len(col2)):
            if type(df_r3[col2[i]].values[0]) == np.int64:
                df_r3.loc[:,col1[i]] = 0
            elif type(df_r3[col2[i]].values[0]) == np.float64:
                df_r3.loc[:,col1[i]] = 0.0
            elif type(df_r3[col2[i]].values[0]) == str:
                df_r3.loc[:,col1[i]] = ''
        tmp.append(df_r3)

    # merge the data
    if len(tmp) > 0:
        df = pd.concat(tmp)
        tmp = []
    else:
        bpy.append_to_log('warning no primers found', log_file)
        return []
    
    cols = ['fwd_query_id','rev_query_id','fwd_pdir','rev_pdir']
    df.drop(columns=cols, inplace = True) # drop some useless info
    df = df.rename(columns={'fwd_database_id':'fwd_primer','rev_database_id':'rev_primer'})
    
    # only output the following info
    if compact:
        cols = ['id','fwd_primer','fwd_orientation','fwd_q_start','fwd_q_end',
                     'fwd_match','fwd_tot','fwd_AS','fwd_CIGAR',
                     'rev_primer','rev_orientation','rev_q_start','rev_q_end',
                     'rev_match','rev_tot','rev_AS','rev_CIGAR']
        df = df[cols]
    
    # print some info about primer matching
    if verbose:
        bpy.append_to_log(str(len(df_f1))+'/'+str(len(reads))+' reads with fwd and rev primer', log_file)
        bpy.append_to_log(str(len(df_f2))+'/'+str(len(reads))+' reads with only fwd primer found', log_file)
        bpy.append_to_log(str(len(df_r3))+'/'+str(len(reads))+' reads with only rev primer found', log_file)
        bpy.append_to_log(str(len(reads) - len(df_f1) - len(df_f2) - len(df_r3))+'/'+str(len(reads))+
                                ' reads with no primer matches', log_file)
    return df

def trim_consensus(df_pmatch, df_cons, log_file = 'log.txt'):
    '''
    This function trims the primers off the consensus reads

    df_pmatch = output dataframe from match_primer_to_reads()
    df_cons = output from get_consensus()

    returns dataframe with columns [id, N frags, fwd_pmr_seq, sequence, rev_pmr_seq]
    '''
    data = df_pmatch.merge(df_cons, on = 'id', how = 'left')
    c1 = data['fwd_primer'].notna() # has fwd primer
    c2 = data['rev_primer'].notna() # has rev primer
    
    out = []
    # trimming for fwd and rev primer
    df = data[c1 & c2][['id','sequence','fwd_orientation','fwd_q_start','fwd_q_end',
                                        'rev_orientation','rev_q_start','rev_q_end']].values
    for i in range(0,len(df)):
        s1 = int(df[i,3])
        s2 = int(df[i,4])
        fwd_pmr = df[i,1][s1:s2]
        s3 = int(df[i,6])
        s4 = int(df[i,7])
        rev_pmr = df[i,1][s3:s4]
        # do slicing for forward orientation
        if df[i,2] == '+':
            seq = df[i,1][s2:s3]
        # process reverse orientation
        else:
            seq = bpy.dna_revcomp(df[i,1][s4:s1])
            fwd_pmr = bpy.dna_revcomp(fwd_pmr)
            rev_pmr = bpy.dna_revcomp(rev_pmr)
        if len(seq) > 10:
            out.append([df[i,0], fwd_pmr, seq, rev_pmr])

    # fwd primer only
    df = data[c1 & ~c2][['id','sequence','fwd_orientation','fwd_q_start','fwd_q_end',
                                         'rev_orientation','rev_q_start','rev_q_end']].values
    for i in range(0,len(df)):
        s1 = int(df[i,3])
        s2 = int(df[i,4])
        fwd_pmr = df[i,1][s1:s2]
        # do slicing for forward orientation
        if df[i,2] == '+':
            seq = df[i,1][s2:]
        # process reverse orientation
        else:
            seq = bpy.dna_revcomp(df[i,1][:s1])
            fwd_pmr = bpy.dna_revcomp(fwd_pmr)
        # only append if seq > 10 nt
        if len(seq) > 10:
            out.append([df[i,0], fwd_pmr, seq, ''])
    
    # rev primer only
    df = data[~c1 & c2][['id','sequence','fwd_orientation','fwd_q_start','fwd_q_end',
                                         'rev_orientation','rev_q_start','rev_q_end']].values
    for i in range(0,len(df)):
        s3 = int(df[i,6])
        s4 = int(df[i,7])
        rev_pmr = df[i,1][s3:s4]
        # do slicing for rev orientation
        if df[i,2] == '+':
            seq = df[i,1][:s3]
        # process reverse orientation
        else:
            seq = bpy.dna_revcomp(df[i,1][s4:])
            rev_pmr = bpy.dna_revcomp(rev_pmr)
    if len(seq) > 10:
        out.append([df[i,0], '', seq, rev_pmr])
    
    df = pd.DataFrame(out, columns = ['id','fwd_primer_seq','consensus','rev_primer_seq'])
    return df.merge(df_cons[['id','N frags','msa_input_file']])
    #return df.merge(df_cons[['id','N frags']])

def clust_pseudoref(df_d, pw_config = '-k15 -w10 -p 0.9 -D', msa_config = '-n -15 -g -10 -l 0 -r 2',
                    metric = 'match_score', linkage = 'complete', thresh = 0.7,
                    workspace = './clust_pseudoref/', pw_output = './pw_out/', log_file = 'log.txt'):
    '''
    Do hierarchical clustering on pseudoreference sequences to reduce search time
    thresh = linkage threshold to cluster. If thresh > 1, this becomes n_clusters.
    '''
    # get pairwise distance info
    files = bpy.get_pairwise_distance(df_d, block = 500, output_folder = pw_output, workspace = workspace,
                                  config = pw_config, symmetric = True)
    pw = []
    for f in files:
        pw.append(pd.read_csv(f))
    subprocess.call(['rm','-r',pw_output])
    pw = pd.concat(pw).reset_index()

    # get cosimilarity matrix
    pw = bpy.get_best(pw, ['query_id','database_id'], metric)
    dct = np.unique(df_d['id'])
    dct = {dct[i]:i for i in range(0,len(dct))}
    pw = bpy.get_feature_vector(pw[['query_id','database_id',metric]], d_dct = dct, q_dct = dct)
    pw = bpy.get_symmetric_matrix(pw, sym_larger = False)
    # add back diagonal values if -D was used in config
    y = pw[pw.columns[pw.columns!='id']].values
    for i in range(0,len(y)):
        y[i,i] = 1
    y = bpy.dist_cosine(y,y)
    y = pd.DataFrame(y, columns = ['d'+str(i) for i in range(0,len(y))])
    y['id'] = pw['id'].values

    # do clustering
    if thresh > 1:
        df_cl = bpy.cluster_hierarchical(y, metric = 'precomputed', linkage = linkage,
                                         thresh = None, n_clusters = thresh)
    else:
        df_cl = bpy.cluster_hierarchical(y, metric = 'precomputed', linkage = linkage,
                                         thresh = thresh, n_clusters = None)
    # compress via multi-seq alignment
    N_cl = len(np.unique(df_cl['cluster_id']))-1
    N_uc = np.sum(df_cl['cluster_id']==-1)
    tot = N_cl + N_uc
    bpy.append_to_log('clustered = '+str(N_cl)+' unclustered='+str(N_uc)+
                      ' total_seq='+str(len(df_cl))+' thresh='+str(thresh),log_file)
    
    df_cl = df_cl.merge(df_d, on = 'id', how = 'left')
    c = df_cl['cluster_id'] > -1
    if len(df_cl[c]) > 0:
        df_m = bpy.cluster_spoa_merge(df_cl[c], config = msa_config, workspace = workspace,
                                      batch_size = 100, cleanup = True, log_file = log_file)
        if len(df_m) > 0:
            df_d = pd.concat([df_cl[~c][['id','sequence']], df_m])
        '''
        # sample instead of merging
        df_m = []
        for i in np.unique(df_cl[c]['cluster_id']):
            t = df_cl[df_cl['cluster_id']==i].values
            ridx = np.random.permutation(len(t))
            df_m.append(t[ridx[0]])
        out = pd.DataFrame(out, columns=df_cl.columns)
        df_d = pd.concat([df_cl[~c][['id','sequence']], out[['id','sequence']]])
        '''
    else:
        bpy.append_to_log('clustering did not reduce sequence space', log_file)
    return df_d

def check_folder(folder):
    # for error checking of folder name to make sure what program expects is consistent
    folder = folder.split('/')
    out = ''
    for j in folder:
        if j!='':
            out+=j+'/'
    return out

def print_help(argv):
    print('Usage: '+argv[0]+' [options] -i fastq_data_folder primers.csv -o output.csv')
    print('''
    Options:
    -h or --help           print help
    
    -i <folder>            folder for input fastq data

       <file>              csv file of forward and reverse primers used
                           primer csv must have the following columns:
                           [fwd_id,fwd_seq,rev_id,rev_seq]
    
    -o <file>              where consensus csv information is written
    
    -db <file>             provides a file containing custom reference sequences to search for.
                           valid inputs:
                           <output.csv> with columns [id, sequence]
                           
                           The script will try generating its own reference database via primer info
                           if this info is not given

    -ff <folder>           folder where frag info are saved
    
    -mf <folder>           folder where multi-sequence fasta and phred info is stored
    
    -log <file>            file to save run logs info

    -prfg                  Run pseudo reference database generator
                           if db <file> is provided, this option is overridden and 
                           the database file is used as the pseudo reference database for read mapping

    -fs <array>            array of sizes ranges to check for forward and reverse primers when generating the
                           pseudo reference database.
                           Default = [[500,1000]]
                           Examples:
                           -fs [[500,1000],[1200,2000]]

    -prf_raw               Do not do clustering on generated pseudo reference sequences
                           Output raw sequences only
                           
    -prf_th x              if x > 1, x = number of clusters for hierarchical complete linkage clustering
                           if x < 1, x = threshold linkage distance used

    -fgs                   run fragment finder, requires fastq data


    -msa                   run multi-sequence alignment, requires frag data

    -cons                  run consensus generation, requires msa data

    -pmf <file>            output for primer matching info

    -fpmr                  run matching of primers to consensus reads, requires consensus read <output.csv>

    -trmc                  trim primers from consensus reads, requires consensus read <output.csv> and <primer match file>
    
    -clst                  run clustering on consensus sequences
    
    -cf <folder>           folder for clustering data
    
    -cin <file>            input file for clustering, otherwise it defaults to trimmed_consensus.csv.gz
    
    -cout <file>           output file for cluster center sequences

    -cl_param a b c d
                        a = min_cluster_size
                        b = thresh_split = thresh accuracy of lower quartile for splitting during clustering
                        c = thresh_merge = thresh cosimilarity for merging during clustering
                        d = N_rsamp = how many random samples to take from subcluster
                            Lower number means faster run time at the cost of sensitivity for
                            rare clusters
    ''')
    sys.exit(1)
#########################################################################################################
# Actual sequence of code in pipeline begins here
def main(argv):
    if len(argv)<2:
        print_help(argv)

    # Important file and folder defaults
    fastq_folder = './fastq/'
    primer_file = 'primers.csv'
    db_file = ''
    pseudodb_file = 'pseudodb.csv.gz'
    pseudo_config = '[[500,1200]]'
    frag_folder = './frags/'
    frag_config = '-k10 -w1'
    msa_folder = './msa/'
    msa_config = '-n -15 -g -10 -l 0 -r 2'
    cons_file = 'consensus.csv.gz'
    pmatch_file = 'primer_match.csv.gz'
    cout_file = 'clusters.csv.gz'
    cluster_folder = './clusters/'
    log_file = 'log.txt'
    run_pseudoref_gen = False
    run_pseudoref_clust = True
    run_frag = False
    run_msa = False
    run_consensus = False
    run_primer_match = False
    run_trim_primers = False
    run_cluster = False
    min_cluster_size = 5
    prf_thresh = 0.50
    thresh_merge = 0.9
    thresh_split = 0.9
    N_rsamp = 2000
    
    # parse some arguments from user
    print(argv)
    bpy.append_to_log('Running the pipeline', log_file)
    for i in range(0,len(argv)):
        if argv[i] == '-h' or argv[i] == '--help': print_help(argv)
        elif argv[i] == '-i':
            fastq_folder = check_folder(argv[i+1])
            primer_file = argv[i+2]
        elif argv[i] == '-o': cons_file = argv[i+1]
        elif argv[i] == '-pmf': pmatch_file = argv[i+1]
        elif argv[i] == '-db': db_file = argv[i+1]
        elif argv[i] == '-log': log_file = argv[i+1]
        elif argv[i] == '-ff': frag_folder = check_folder(argv[i+1])
        elif argv[i] == '-mf': msa_folder = check_folder(argv[i+1])
        elif argv[i] == '-prfg': run_pseudoref_gen = True
        elif argv[i] == '-prf_raw': run_pseudoref_clust = False
        elif argv[i] == '-prf_th': prf_thresh = argv[i+1]
        elif argv[i] == '-fs': pseudo_config = argv[i+1]
        elif argv[i] == '-fgs': run_frag = True
        elif argv[i] == '-msa': run_msa = True
        elif argv[i] == '-cons': run_consensus = True
        elif argv[i] == '-fpmr': run_primer_match = True
        elif argv[i] == '-trmc': run_trim_primers = True
        elif argv[i] == '-clst': run_cluster = True
        elif argv[i] == '-cin': cin_file = argv[i+1]
        elif argv[i] == '-cout': cout_file = argv[i+1]
        elif argv[i] == '-cf': cluster_folder = argv[i+1]
        elif argv[i] == '-clst_param':
            min_cluster_size = int(argv[i+1])
            thresh_merge = float(argv[i+2])
            thresh_split = float(argv[i+3])
            N_rsamp = int(argv[i+4])

    # if nothing is select, then default is to run everything
    if (run_msa + run_frag + run_pseudoref_gen
        + run_consensus + run_primer_match + run_trim_primers + run_cluster)==False:
        run_pseudoref_gen = True
        run_frag = True
        run_msa = True
        run_consensus = True
        run_primer_match = True
        run_trim_primers = True
        run_cluster = True
    
    # Print run settings
    bpy.append_to_log('fastq_folder     = '+fastq_folder, log_file)
    bpy.append_to_log('primer_file      = '+primer_file, log_file)
    if db_file!='':
        run_pseudoref_gen = False
        bpy.append_to_log('db_file          = '+db_file, log_file)
    else:
        bpy.append_to_log('pseudodb_file    = '+pseudodb_file, log_file)
        bpy.append_to_log('pseudo_fs        = '+pseudo_config, log_file)
    bpy.append_to_log('log_file         = '+log_file, log_file)
    bpy.append_to_log('frag_data_folder = '+frag_folder, log_file)
    bpy.append_to_log('msa_data_folder  = '+msa_folder, log_file)
    bpy.append_to_log('cons_file        = '+cons_file, log_file)
    cin_file = 'trimmed_'+cons_file
    bpy.append_to_log('cin_file         = '+cin_file, log_file)
    bpy.append_to_log('cout_file        = '+cout_file, log_file)
    bpy.append_to_log('cluster_folder   = '+cluster_folder, log_file)
    bpy.append_to_log('Run Pseudoref generator = '+str(run_pseudoref_gen), log_file)
    bpy.append_to_log('Cluster pseudo ref      = '+str(run_pseudoref_clust), log_file)
    bpy.append_to_log('Run Fragment Search     = '+str(run_frag), log_file)
    bpy.append_to_log('Run Multi-Seq Alignment = '+str(run_msa), log_file)
    bpy.append_to_log('Run Consensus           = '+str(run_consensus), log_file)
    bpy.append_to_log('Run primer match        = '+str(run_primer_match), log_file)
    bpy.append_to_log('Run trim primers        = '+str(run_trim_primers), log_file)
    bpy.append_to_log('Run clustering          = '+str(run_cluster), log_file)
    bpy.append_to_log('prf: thresh             = '+str(prf_thresh), log_file)    
    bpy.append_to_log('clust: min_cluster_size = '+str(min_cluster_size), log_file)
    bpy.append_to_log('clust: thresh_merge     = '+str(thresh_merge), log_file)
    bpy.append_to_log('clust: thresh_split     = '+str(thresh_split), log_file)
    bpy.append_to_log('clust: N_rsamp          = '+str(N_rsamp), log_file)
    
    # Load the fastq data
    if run_frag or run_msa or run_pseudoref_gen:
        bpy.append_to_log('loading fastq data', log_file)
        files = glob.glob(fastq_folder+'*.fastq')
        reads = bpy.load_basecalled_data(files, log_file)
    
    # use reference library sequences if it exists
    if db_file!='':
        bpy.append_to_log('Loading '+db_file+' as reference database', log_file)
        ftype = db_file.split('.')[1]
        if ftype == 'csv':
            ref = pd.read_csv(db_file)
            ref = ref[['id','sequence']]
        elif ftype == 'fa' or ftype == 'fasta':
            ref = bpy.read_fasta(db_file)
        elif ftype == 'fq' or ftype == 'fastq':
            ref = bpy.read_fastq(db_file)
    else:
        db_file = pseudodb_file

    # run pseudo reference database generator if selected
    if run_pseudoref_gen:
        # load primers
        primers = bpy.load_primers(primer_file)

        # generate pseudo reference database using primer info
        bpy.append_to_log('Generating pseudo reference database', log_file)
        ref = generate_pseudo_refdb(primers, reads, block=100000, f_sizes = pseudo_config, log_file = log_file)
        ref.to_csv(db_file, index=False, compression='infer')
        db_file = pseudodb_file

        if run_pseudoref_clust:
            ref = pd.read_csv(db_file)
            n_iter = 1
            for k in range(0,n_iter):
                ref = clust_pseudoref(ref, metric = 'match_score', linkage = 'complete', thresh = prf_thresh,
                            pw_config = '-k15 -w10 -p 0.9 -D', msa_config = '-n -15 -g -10 -l 0 -r 2',
                            workspace = './clust_pseudoref/', pw_output = './pw_out/', log_file = log_file)
                #ref.to_csv(db_file+'.'+str(k), index=False, compression='infer')
            ref.to_csv(db_file, index=False, compression='infer')

    # Find RCA fragments with the aligner
    if run_frag:
        # load reference database
        if os.path.exists(db_file):
            ref = pd.read_csv(db_file)
        else:
            bpy.append_to_log('Reference database '+db_file+' not found', log_file)
            sys.exit(1)
            
        # make data directory
        bpy.check_dir(frag_folder)
        
        # search or RCA fragments in reads
        bpy.append_to_log('Searching for RCA frags in reads', log_file)
        find_RCA_frags(ref, reads, block=10000, folder = frag_folder, config = frag_config, log_file = log_file)
        bpy.append_to_log('Aligner done', log_file)

    if run_msa:
        # Loading the data for multi-sequence alignment
        bpy.append_to_log('Loading frags from '+frag_folder, log_file)
        files = glob.glob(frag_folder+'fraginfo_*.csv.gz')
        out = []
        for i in files:
            out.append(pd.read_csv(i))
        frags = pd.concat(out)
        
        # filtering out overlapping frags
        frags = bpy.remove_overlaps(frags, metric = 'AS', thresh = 50, log_file = 'log.txt')
        
        # annotating information about the reads
        bpy.append_to_log('annotating reads', log_file)
        df_frags = annotate_reads(frags, log_file)
        frags = frags.merge(df_frags[['query_id','N frags']], on = 'query_id', how = 'left')
        bpy.append_to_log('annotation complete', log_file)
        
        # Run multi-sequence alignment on fragments
        bpy.append_to_log('Running multi-sequence aligner', log_file)
        perform_spoa_MSA(reads, frags[frags['N frags'] > 1], batch_size = 100,
                         config = msa_config, thread_lock = False,
                         folder = msa_folder, log_file = log_file)
        bpy.append_to_log('multi-sequence alignment done', log_file)
    
    # Generate consensus sequence and save the info
    if run_consensus:
        MSA_outfile = glob.glob(msa_folder+'*.out')
        bpy.append_to_log('Number of MSA_outfiles='+str(len(MSA_outfile)), log_file)        
        df = get_spoa_consensus(MSA_outfile, log_file)

        # Write the data to file
        bpy.append_to_log('Writing consensus output to: '+cons_file, log_file)
        df.to_csv(cons_file, index=False, compression='infer')
        bpy.append_to_log('Consensus generation done', log_file)
        
    # match primers to read consensus
    if run_primer_match:
        # load primer info
        primers = bpy.load_primers(primer_file)

        # load consensus info
        bpy.append_to_log('Loading cons_file: '+cons_file, log_file)
        df = pd.read_csv(cons_file)
        
        # match the primers to consensus reads        
        bpy.append_to_log('Matching primers to reads', log_file)
        df = df.rename(columns={'consensus':'sequence'})
        df = match_primer_to_read(primers, df[['id','sequence']], log_file = log_file)
       
        # Write the data to file
        bpy.append_to_log('Writing primer_match to '+pmatch_file, log_file)
        df.to_csv(pmatch_file, index=False, compression='infer')
        bpy.append_to_log('primer matching done', log_file)

    if run_trim_primers:
        # load data
        bpy.append_to_log('Loading pmatch_file: '+pmatch_file, log_file)
        df_pmatch = pd.read_csv(pmatch_file)
        bpy.append_to_log('Loading cons_file: '+cons_file, log_file)
        df_cons = pd.read_csv(cons_file)
        df_cons = df_cons.rename(columns = {'consensus':'sequence'})

        # run the trimming operation
        bpy.append_to_log('Trimming the reads', log_file)
        df = trim_consensus(df_pmatch, df_cons, log_file = 'log.txt')
        df = df.merge(df_cons[['id','N frags']], on = 'id', how = 'left')

        # Write the data to file
        bpy.append_to_log('Writing trimmed consensus sequence', log_file)
        df.to_csv(cin_file, index=False, compression='infer')
        bpy.append_to_log('primer matching done', log_file)
    
    if run_cluster:
        # load data
        bpy.append_to_log('Loading cin_file: '+cin_file, log_file)
        df = pd.read_csv(cin_file).rename(columns = {'consensus':'sequence'})

        # only work on sequences which have both forward and reverse primers
        if 'fwd_primer_seq' in df.columns and 'rev_primer_seq' in df.columns:
            df = df[df['fwd_primer_seq'].notna() & df['rev_primer_seq'].notna()]

        # refine the clusters
        data = bpy.cluster_refine(df, config = '-k15 -w10 -p 0.9 -D', ofile = 'rcent',
                                  max_iter = 10, center_size = 20, min_samples = min_cluster_size,
                                  thresh_merge = thresh_merge, thresh_split = thresh_split, N_rsamp = N_rsamp,
                                  workspace = cluster_folder, log_file = log_file)
        
        # save data
        bpy.append_to_log('Writing clustered sequences to '+cout_file, log_file)
        data.to_csv(cout_file, index=False, compression='infer')
        bpy.append_to_log('clustering done', log_file)
        
if __name__ == "__main__":
    main(sys.argv)
