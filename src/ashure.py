#!/usr/bin/env python
# Consensus error correction for nanopore sequencing
# bilgenurb@gmail.com

# Loading libraries
import numpy as np
import scipy as sp
import pandas as pd

# For interfacing with the file system
import glob
import subprocess
import os
import time
import sys
import logging
import argparse
import json

# Loading the custom libraries
import bilge_pype as bpy

#########################################################################################################

def generate_pseudo_refdb(primers, reads, block=10000, fs='500-1000', workspace='./pseudo_refdb/', config='-x map-ont'):
    '''
    Build reference database using paired end search of some primer sequences in the uncorrected reads
    primers   = dataframe containing primer information with columns:
                [fwd_id, fwd_seq, rev_id, rev_seq]
    reads     = dataframe containing read information with columns:
                [id, sequence] or [id, sequence, quality]
    block     = block size of reads to work on each time. This is for saving on RAM.
    fs        = size range of fastq reads to search primers
    workspace = workspace for the aligner
    config    = configs to pass to minimap2 aligner
    returns dataframe of pseudo reference sequences with columns
    '''
    def parse_fsize(text):
        data = [i.split('-') for i in text.split(',')]
        return np.array(data).astype(int)

    def split_sequence(df):
        '''
        Function to get sequence info for each frag
        df = pandas dataframe with columns:
             [id, orientation, start, end, sequence, quality]
        '''
        data = []
        for r_id, ot, t1, t2, seq, qual in df.values:
            seq = seq[t1:t2]
            qual = qual[t1:t2]
            if ot == '-':
                seq = bpy.dna_revcomp(seq)
                qual = qual[::-1]
            data.append([r_id+'_rs'+str(t1)+'_'+str(t2), seq, qual])
        return pd.DataFrame(data, columns=['id', 'sequence', 'quality'])

    workspace = bpy.check_dir(workspace)
    # look for the primers
    data = []
    for s1,s2 in parse_fsize(fs):
        # filter for some read sizes to reduce search time
        s = (reads['length'] > s1) & (reads['length'] < s2)
        # chunk up the query so we can see progress
        N_chunks = np.ceil(len(reads[s])/block)
        for i, df_i in enumerate(np.array_split(reads[s], N_chunks)):
            logging.info('Working on '+str(i)+'/'+str(N_chunks)+' block size='+str(block)+' read size='+str(s1)+' to '+str(s2)+' bp')
            out = match_primer_to_read(primers, df_i, workspace=workspace, config=config)            
            if len(out) > 0:
                data.append(out)
    if len(data)==0:
        logging.error('generate_pseudo_refdb: no primers found in reads.')
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
    data = pd.DataFrame(x, columns=['id','orientation','start','end'])

    # merge and get pseudo reference sequences
    reads = bpy.add_seq(reads[reads['id'].isin(data['id'])]) # low mem optimization
    data = data.merge(reads[['id','sequence','quality']], on = 'id', how = 'left')
    df_d = split_sequence(data[['id','orientation','start','end','sequence','quality']]).drop_duplicates()
    logging.info('Reference sequences found '+str(len(df_d)))
    return df_d
    
def find_RCA_frags(ref, reads, block=10000, output='./frags/', config=''):
    '''
    Function to search sequencing primers and identify RCA reads from fastq info
    ref            = similarily related sequences to search repeated RCA frags with columns [id, sequence]
    reads          = fastq info
    block          = number of fastq reads to process each time
    output         = output where data is saved
    '''
    output = bpy.check_dir(output)
    workspace = output+'aligner'

    N_chunks = np.ceil(len(reads)/block)
    build_index = True
    for i, reads in enumerate(np.array_split(reads, N_chunks)):
        logging.info('Working on '+str(i)+'/'+str(N_chunks)+' block size='+str(block))
        frags = bpy.run_minimap2(reads, ref, workspace=workspace, config=config, cigar=True, build_index=build_index, use_index=True)
        build_index = False
        logging.info('Saving information')
        if len(frags) > 0:
            frags.to_csv(output+'fraginfo_'+str(i)+'.csv.gz', index=False, compression='infer')
        logging.info('frags found in this block = '+str(len(frags)))
    # clear the work folder
    subprocess.run(['rm','-r',workspace]) # remove the aligner folder after work is done

def annotate_reads(df):
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
    data = df[['query_id','q_len','q_start','q_end','t_len','t_start','t_end','orientation','AS','match_score']].values
    
    # storage variables
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
            N_1D2+=1
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
            df_frag.append([data[i,0],N_frags,N_1D2,data[i,1],q_unmapped,db_fwd_cover,db_rev_cover,match_size,np.mean(match_score),np.mean(AS),np.std(AS)])
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
    col = ['query_id','N frags','1D2','q_len','q_unmapped','db_fwd_cover','db_rev_cover','match_size','avg_match','avg_AS','std_AS']
    return pd.DataFrame(df_frag, columns=col)

def perform_MSA(df, frags, batch_size=100, folder='./msa/', thread_lock=True, config='-l 0 -r 2'):
    '''
    This function performs multisequence alignment on RCA reads
    df     = dataframe containing the original read fastq information
    frags  = RCA fragments with start and stop locations found by the previous aligner
    output = folder containing file of aligned RCA reads
    '''
    # make workspace folder if it does not exist
    folder = bpy.check_dir(folder)
    files = glob.glob(folder+'*')
    bpy.batch_file_remove(files)
    # batch submit the jobs
    frags = frags.rename(columns={'query_id':'id'})
    cols = ['q_len','q_start','q_end']
    for c in cols:
        frags[c] = frags[c].astype(int)
    # low mem check
    if 'filename' in df.columns:
        rlist = df[df['id'].isin(frags['id'])].sort_values(by=['filename'])['id'].drop_duplicates().values
    else:
        rlist = np.unique(frags['id'])
    N_chunks = np.ceil(len(rlist)/batch_size)
    for j, rid in enumerate(np.array_split(rlist, N_chunks)):
        # merge the input data
        data = frags[frags['id'].isin(rid)][['id','q_len','q_start','q_end','orientation']]
        # check if we are doing low_mem
        reads = bpy.add_seq(df[df['id'].isin(rid)])
        data = data.merge(reads[['id','sequence','quality']], on='id', how='left')
        data = data.sort_values(by=['id','q_start']).values
        MSA_infile = []
        seqlist = []
        k = 0
        for i in range(0,len(data)):
            s1 = data[i,2]
            s2 = data[i,3]
            # if read orientations match and gaps are not big, try merging them
            if i > 0 and data[i,4]==data[i-1,4]:
                g1 = data[i,2] - data[i-1,3]
                if g1*4 < data[i,3]-data[i,2]:
                    s1 = int((data[i,2]+data[i-1,3])/2)
            if i+1 < len(data) and data[i,4]==data[i+1,4]:
                g2 = data[i+1,2] - data[i,3]
                if g2*4 < data[i,3]-data[i,2]:
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
        if len(MSA_infile) > 0:
            logging.info('MSA progress = '+str(j)+'/'+str(N_chunks))
            bpy.run_msa(MSA_infile, aligner='spoa', config=config, thread_lock=thread_lock)
            MSA_infile = []
 
def get_spoa_consensus(MSA_outfile):
    '''
    Reads spoa outfile for the consensus from multi-sequence alignment
    MSA_outfile = multi-sequence alignment output from spoa

    output = consensus sequence information
    '''
    logging.info('Reading consensus info from spoa output files')
    data = []
    for i in range(0,len(MSA_outfile)):
        fname = MSA_outfile[i]
        f = open(fname,'r')
        text = f.read().split('\n')
        f.close()
        rcaID = fname.split('/')[-1].split('_Nfrags')[0]
        nfrags = int(fname.split('/')[-1].split('_Nfrags')[1].split('.')[0])
        if len(text) < 2:
            logging.info('file read error on '+fname)
        else:
            data.append([rcaID, fname.split('.out')[0], nfrags, text[1]])
        if i%20000 == 0:
            logging.info('Consensus progress = '+str(i)+'/'+str(len(MSA_outfile)))
    return pd.DataFrame(data, columns=['id','msa_input_file','N frags','consensus'])

def match_primer_to_read(primers, reads, thresh=10, config='-x map-ont', workspace='./match_pmr_to_rd/', compact=True, verbose=True):
    '''
    Takes a list of primer pairs and reads and matches primer pairs to the reads
    primers   = pandas dataframe with columns = [fwd_id, fwd_sequence, rev_id, rev_sequence]
    reads     = pandas dataframe with columns = [id, sequence]
    config    = config passed to minimap2
    workspace = workspace used by the aligner
    output    = pandas dataframe with primer mapping info and sequence id as the first column
    '''
    # search for forward primers
    logging.info('searching forward primers')
    database = primers.rename(columns={'fwd_id':'id','fwd_seq':'sequence'})
    df1 = bpy.run_minimap2(reads, database[['id','sequence']], workspace, config, build_index=True, use_index=True)
    # search for reverse primers
    logging.info('searching reverse primers')
    database = primers.rename(columns={'rev_id':'id','rev_seq':'sequence'})
    df2 = bpy.run_minimap2(reads, database[['id','sequence']], workspace, config, build_index=True, use_index=True)
    # remove the working directory
    subprocess.run(['rm','-r',workspace])

    # merge the primers that were found
    if len(df1) == 0 and len(df2)==0:
        logging.info('warning no fwd and rev primers found')
        return []
    elif len(df1) == 0:
        logging.info('warning no fwd primers found')
        df2['pdir'] = 'fwd'
        df = df2
    elif len(df2) == 0:
        logging.info('warning no rev primers found')
        df1['pdir'] = 'rev'
        df = df1
    else:
        logging.info('fwd and rev primers found')
        df1['pdir'] = 'fwd'
        df2['pdir'] = 'rev'
        df = pd.concat([df1,df2]) # concatenate
  
    # finding best matches
    df = bpy.remove_overlaps(df, 'AS', thresh)
    
    logging.info('finding best primer match for each read')
    # do a sort by query and position of the primers
    df = df.sort_values(by=['query_id','q_start']).reset_index(drop=True)
    data = df[['query_id','database_id','pdir','orientation','q_start','q_end','AS']].values
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
        logging.warning('warning no primers found')
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
    
    # some info about primer matching
    if verbose:
        logging.info(str(len(df_f1))+'/'+str(len(reads))+' reads with fwd and rev primer')
        logging.info(str(len(df_f2))+'/'+str(len(reads))+' reads with only fwd primer found')
        logging.info(str(len(df_r3))+'/'+str(len(reads))+' reads with only rev primer found')
        logging.info(str(len(reads) - len(df_f1) - len(df_f2) - len(df_r3))+'/'+str(len(reads))+' reads with no primer matches')
    return df

def trim_consensus(df_pmatch, df_cons):
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
    
    df = pd.DataFrame(out, columns=['id','fwd_primer_seq','consensus','rev_primer_seq'])
    return df.merge(df_cons[['id','N frags','msa_input_file']])

def perform_cluster(df_q, df_d=[], max_iter=1, csize=20, N=2000, th_s=2, th_m=0.9, pw_config='', msa_config='', workspace='./cluster/', track_file='', timestamp=True):
    '''
    Refines the cluster centers for a set of sequence data
    df_q = dataframe of query sequences to search for cluster centers.
    df_d = dataframe of cluster centers with you want to append to
           df_q and df_d must both have at least columns [id, sequence]
    csize = number of sequences from cluster center to use for recomputation of cluster center sequence
    max_iter = max iterations to run
    th_s = threshold lower quartile of accuracy to sample for rare sequences
    th_m = threshold eps for clustering related cluster centers
    N = size of the random sample to compute pairwise distance data with
    workspace = workspace folder for the aligner
    '''
    # initialize cluster centers
    if len(df_d) > 0:
        df_c = df_d
    else:
        df_c = bpy.cluster_sample(df_q, N=5, th=0.8, rsize=1, config=pw_config, workspace=workspace)
    df_c['split'] = True # refine all new centers by default
    track_file_list = []
    for k in range(0, max_iter):
        logging.info('perform_cluster: iter = '+str(k)+'/'+str(max_iter))
        df_c = cluster_split(df_q, df_c, N, csize, pw_config, msa_config, workspace)
        df_c = cluster_sweep(df_q, df_c, th_s, N, csize, pw_config, msa_config, workspace)
        df_c = cluster_merge(df_q, df_c, th_m, 100, csize, pw_config, msa_config, workspace)
        # track progress
        if track_file!= '':
            import datetime
            ts = ''
            if timestamp:
                ts = '_'+datetime.datetime.fromtimestamp(time.time()).strftime('%Y_%m_%d_%H_%M_%S')
            fname = track_file+ts+'_i'+str(k)+'_c'+str(len(df_c))+'.csv.gz'
            df_c.to_csv(fname, index=False, compression='infer')
            track_file_list.append(fname)
    # final refinement
    df_c['split'] = True
    df_c = cluster_split(df_q, df_c, 200, csize, pw_config, msa_config, workspace)
    df_c = cluster_merge(df_q, df_c, th_m, 100, csize, pw_config, msa_config, workspace)
    logging.info('perform_cluster: complete')
    return df_c, track_file_list

def cluster_eval(df_q, df_c, pw_config, workspace, cigar=True):
    config = pw_config+' --for-only'
    df_align = bpy.run_minimap2(df_q, df_c, config=config, cigar=cigar, workspace=workspace)
    cols = np.unique(df_align['database_id'])
    df_c = df_c[df_c['id'].isin(cols)]
    logging.info('cluster_eval: number of clusters = '+str(len(df_c)))
    return df_c, df_align

def cluster_sweep(df_q, df_c, th_s, N, csize, pw_config, msa_config, workspace):
    df_c, df_align = cluster_eval(df_q, df_c, pw_config, workspace, cigar=False)
    # normalize s1
    df_a = bpy.get_best(df_align, ['database_id'], metric='s1', stat='idxmax').rename(columns={'s1':'b1'})
    df_align = df_align.merge(df_a[['database_id','b1']],on='database_id',how='left')
    df_align['m1'] = df_align['s1']/df_align['b1']
    # partition and get lower quartile # debug
    df_align = bpy.get_best(df_align,['query_id'],metric='m1',stat='idxmax')
    df_align = df_align.rename(columns={'query_id':'id'})
    x = bpy.cluster_Kmeans(df_align[['id','m1']], n_clusters=th_s, n_init=10, n_iter=100, ordered=True)
    c = x['cluster_id'] == 0 # get bad alignments
    uncover = set(df_q['id'].astype(str)) - set(x[~c]['id'].astype(str))
    uncover = df_q[df_q['id'].isin([i for i in uncover])]
    logging.info('cluster_sweep: uncovered '+str(len(uncover))+'/'+str(len(df_q)))
    # subsample and run optics to get centers
    ridx = np.random.permutation(len(uncover)) 
    x = uncover.iloc[ridx[:N]]
    cout = cluster_compute(x, csize, outliers=False, pw_config=pw_config, msa_config=msa_config, workspace=workspace)
    cout['split'] = True
    df_c = df_c.append(cout[['id','sequence','split']])
    df_c['id'] = ['cluster'+str(i) for i in range(0,len(df_c))]
    return df_c 

def cluster_split(df_q, df_c, N, csize, pw_config, msa_config, workspace):
    df_c, df_align = cluster_eval(df_q, df_c, pw_config, workspace, cigar=False)
    squeue = df_c[df_c['split']]['id'].values 
    # mark clusters for splitting
    keep = set(df_c['id'].astype(str)) - set(squeue)
    df_c = df_c[df_c['id'].isin([i for i in keep])]
    # add weights for semi random sampling
    df_align['weight'] = cluster_weight(df_align, metric='s1')
    # do the split
    for i in range(0,len(squeue)):
        qout = squeue[i]
        logging.info('cluster_split: splitting on cid='+qout+' '+str(i)+'/'+str(len(squeue)))
        # subsample
        qlist = cluster_subsample(df_align, [qout], N=N, mode='semi', metric='s1')
        qlist = df_q[df_q['id'].isin(qlist)]
        cout = cluster_compute(qlist, csize, outliers=False, pw_config=pw_config, msa_config=msa_config, workspace=workspace)
        # record cluster centers
        cout['split'] = False
        if len(cout) > 1:
            cout['split'] = True
        df_c = df_c.append(cout[['id','sequence','split']])
    df_c['id'] = ['cluster'+str(i) for i in range(0,len(df_c))]
    return df_c 

def cluster_merge(df_q, df_c, th_m, N, csize, pw_config, msa_config, workspace):
    df_c, df_align = cluster_eval(df_q, df_c, pw_config, workspace, cigar=True)
    # get the cosimilarity matrix
    vec = bpy.get_feature_vector(df_align[['query_id','database_id','AS']])
    cols = vec.columns[vec.columns!='id']
    x = vec[cols].values
    cosim = bpy.dist_cosine(x.T, x.T)
    cosim = pd.DataFrame(cosim, columns=cols)
    cosim['id'] = cols
    df_clst = bpy.cluster_OPTICS(cosim, metric='precomputed', min_samples=2, max_eps=th_m)
    # clusters not to perform merge
    c = df_clst['cluster_id']== -1
    df_c = df_c[df_c['id'].isin(df_clst[c]['id'])]
    mqueue = []
    # extract clusters to merge
    for cid in np.unique(df_clst[~c]['cluster_id']):
        x = df_clst[df_clst['cluster_id'] == cid]['id'].values
        mqueue.append([j for j in x])
    df_align['weight'] = cluster_weight(df_align, metric='AS')
    # merge clusters
    logging.info('cluster_merge: '+str(np.sum(~c))+'/'+str(len(df_c))+' clusters to merge')
    for i in range(0,len(mqueue)):
        mout = mqueue[i]
        logging.info('cluster_merge: doing merging on '+str(len(mout))+' clusters, '+str(i)+'/'+str(len(mqueue)))
        qlist = cluster_subsample(df_align, mout, N=N, mode='sorted',metric='AS')
        qlist = df_q[df_q['id'].isin(qlist)]
        cout = cluster_compute(qlist, csize, outliers=False, pw_config=pw_config, msa_config=msa_config, workspace=workspace)
        cout['split']=True
        df_c = df_c.append(cout[['id','sequence','split']])
    df_c['id'] = ['cluster'+str(i) for i in range(0,len(df_c))]
    return df_c

def cluster_weight(df_align, metric='s1'):
    df_a = bpy.get_best(df_align, ['query_id'], metric=metric, stat='idxmax').rename(columns={metric:'b1'})
    df_align = df_align.merge(df_a[['query_id','b1']],on='query_id',how='left')
    w = df_align[metric]/(df_align[metric] + df_align['b1'])
    return w

def cluster_subsample(df_align, dlist, N=100, mode='sorted', metric='AS'):
    df_align = df_align.sort_values(by=['database_id',metric])[::-1] # set to high value first
    x = []
    for d in dlist:
        if mode=='rand':
            v = df_align[df_align['database_id']==d]['query_id'].values
            ridx = np.random.permutation(len(v))
            qlist = v[ridx[:N]]
        elif mode=='semi':
            v = df_align[df_align['database_id']==d][['query_id','weight']].values
            v[v[:,1] >= 0.5, 1] = 1
            v = v[v[:,1] > np.random.rand(len(v))] # weight probability on alignment
            ridx = np.random.permutation(len(v))
            qlist = v[ridx[:N],0]
        else:
            v = df_align[df_align['database_id']==d][['query_id','weight']].values
            v[v[:,1] >= 0.5, 1] = 1
            v = v[v[:,1] > np.random.rand(len(v))] # weight probability on alignment
            qlist = v[:N,0]
        x+= [q for q in qlist]
    return np.unique(x)

def cluster_compute(df_q, csize=20, metric='similarity', outliers=False, pw_config ='-k15 -w10 -D', msa_config='-l 0 -r 2', workspace='./clust_run/'):
    '''
    Find cluster centers for a set of sequence data
    qry = dataframe of query sequences to search for cluster centers.
           Must have at least the following columns:
           [id, sequence]
           [id, sequence, quality]
           [id, sequence, msa_input_file]
           msa_input_files should be in fasta or fastq format
           with filename as .fasta, .fa, .fastq, or .fq 
    csize = number of sequences to take from cluster center for multi-seq alignment
    workspace = workspace folder for the aligner
    '''
    # make sure names are strings
    df_q['id'] = df_q['id'].astype(str)
    if len(df_q) < 2:
        logging.warning('sample size=1, too small to do anything')
        if outliers:
            df_q['outlier'] = True
            return df_q[['id','sequence','outlier']]
        return pd.DataFrame([],columns=['id','sequence'])
    logging.info('cluster_compute: computing pairwise distance matrix')
    config = pw_config+' --dual=no --for-only'
    df = bpy.run_minimap2(df_q, df_q, config=config, workspace=workspace, cleanup=True)

    # catch error in case pairwise alignment fails
    if len(df) == 0:
        logging.warning('cluster_compute: not enough samples to cluster')
        if outliers:
            df_q['outlier'] = True
            return df_q[['id','sequence','outlier']]
        return pd.DataFrame([],columns=['id','sequence'])
    df = bpy.get_best(df, ['query_id','database_id'], 'AS')
    df = bpy.get_feature_vector(df[['query_id','database_id',metric]], symmetric=True)
    df = bpy.get_symmetric_matrix(df, sym_larger=False)
    # convert to distance
    logging.info('preparing precomputed data')
    for i in range(0,len(df)):
        df.iloc[:,i] = 1-df.iloc[:,i]
        df.iat[i,i] = 0
    # do optics on raw data
    logging.info('cluster_compute: running optics')
    df_clst = bpy.cluster_OPTICS(df, metric='precomputed', alt_label=False)
    df_clst = df_clst.merge(df_q, on='id', how='left')
    df_o = df_clst[df_clst['cluster_id'] == -1]
    df_o['outlier'] = True
    df_i = df_clst[df_clst['cluster_id'] > -1]
    # do merge via msa if needed
    if len(df_i) > 0:
        din = []
        for cid in np.unique(df_i['cluster_id']):
            df = df_i[df_i['cluster_id'] == cid]
            din.append(df.iloc[:csize])
        df_i = bpy.cluster_spoa_merge(pd.concat(din), config=msa_config, workspace=workspace, cleanup=True)
        df_i['outlier'] = False
    # return results
    if outliers:
        col = ['id','sequence','outlier']
        if len(df_i) > 0:
            return pd.concat([df_i,df_o])[col]
        else:
            return df_o[col]
    else:
        col = ['id','sequence']
        if len(df_i) == 0: 
            logging.warning('cluster_compute: no clusters found')
        return  df_i[col]

def update_config(default, new_config):
    '''
    Updates default variables with user specified variables
    default = dict contain default variables
    new_config = dict contain new variables
    '''
    # load settings of old config file
    if type(new_config['config_file'])==str and os.path.exists(new_config['config_file']):
        with open(new_config['config_file'], 'r') as f:
            default = json.load(f)
    folders = ['workspace','frag_folder','msa_folder','clst_folder']
    # update settings
    for k in new_config.keys():
        if new_config[k]!=None:
            default[k] = new_config[k]
            if k in folders:
                default[k] = bpy.check_dir(default[k], fname_only=True)
    # apply new settings to old config file
    if type(new_config['config_file'])==str:
        with open(new_config['config_file'], 'w') as f:
            logging.info('writing new settings to '+new_config['config_file'])
            json.dump(default, f, indent=2)
    return default

#########################################################################################################
# Actual sequence of code in pipeline begins here
def main():
    parser = argparse.ArgumentParser(description='aSHuRE: a consensus error correction pipeline for nanopore sequencing',
                        formatter_class=argparse.RawTextHelpFormatter)
    # variable to store default configs
    config = {}

    # add subcommands for modifying the run configuration 
    subparser = parser.add_subparsers(title='subcommands', dest='subcommand')
    run_parser = subparser.add_parser('run', formatter_class=argparse.RawTextHelpFormatter, help='suboptions for running the pipeline')
    config['fastq']=[]
    run_parser.add_argument('-fq', dest='fastq', type=str, nargs='+', help='fastq files from the basecaller')
    config['exclude']=[]
    run_parser.add_argument('-e', dest='exclude', type=str, nargs='+', default=[], help='frags files containing reads which are excluded from search')
    config['primer_file']='primers.csv'
    run_parser.add_argument('-p', dest='primer_file', type=str,
            help='''csv file containing forward and reverse primers used.
    This must have at least columns [fwd_id, fwd_seq, rev_id, rev_seq]''')
    config['db_file']='pseudodb.csv.gz'
    run_parser.add_argument('-db', dest='db_file', type=str,
            help='''reference sequences to search fragments with (optional). If this is not provided prfg is run.
    This file can be a csv with columns [id, sequence] or fastq or fastq format
    If prfg is run, then the new pseudo reference sequences are written to this file''')
    config['workspace']='./workspace/'
    run_parser.add_argument('-w', dest='workspace', type=str, help='workspace folder where frags, msa, clusters, and other output are stored')
    config['cons_file']='consensus.csv.gz'
    run_parser.add_argument('-o1', dest='cons_file', type=str, help='output csv file for consensus sequence information')
    config['pmatch_file']='pmatch.csv.gz'
    run_parser.add_argument('-o2', dest='pmatch_file', type=str, help='output csv file for primer matches to each fastq read')
    config['cin_file']='trimmed_'+config['cons_file']
    run_parser.add_argument('-o3', dest='cin_file', type=str, help='input file for clustering and output for trimmed consensus')
    config['cout_file']='clusters.csv.gz'
    run_parser.add_argument('-o4', dest='cout_file', type=str, help='output csv file for sequence of cluster centers')
    config['log_file']='ashure.log'
    config['low_mem']=False
    run_parser.add_argument('--low_mem', dest='low_mem', action='store_true', help='enable optimizations that reduce RAM used')
    run_parser.add_argument('-log', dest='log_file', type=str, help='log file where pipeline progress is logged')
    config['config_file']=''
    run_parser.add_argument('-c', dest='config_file', type=str,
            help='''Input json file where run configuration is stored.
    If none is provided, the default config is used.''')
    run_parser.add_argument('-r', dest='run', type=str,
            help='''Specify what subsets of this pipeline to run
    prfg = generate pseudo reference database using primer information
    fgs  = find repeated fragments in fastq reads
    msa  = run multi-sequence alignment on repeated reads
    cons = read consensus output after multi-sequence alignment
    fpmr = map primers to each consensus fastq read
    trmc = trim primers from consensus reads
    clst = running clustering on trimmed reads
    
    -r msa,cons runs only multi-seq alignment and consensus readout''')
 
    # config for pseudo reference generator
    prfg_parser = subparser.add_parser('prfg', formatter_class=argparse.RawTextHelpFormatter, help='suboptions for pseudo reference generator')
    config['prfg_fs']='500-1200'
    prfg_parser.add_argument('-fs', dest='prfg_fs', type=str,
            help='''fastq read size ranges to search for primers
    -fs 100-200             searches for fragments in reads of length 100-200bp
    -fs 100-200,500-900     searches for fragments in reads of length 100-200bp and 500-900bp''')
    config['prfg_config']='-k5 -w1 -s 20 -P'
    prfg_parser.add_argument('-c', dest='prfg_config', type=str, help='config passed to minimap2')
    prfg_parser.add_argument('-s', dest='config_file', type=str, help='write settings to configuration file')
    prfg_parser.add_argument('-r', dest='run_prfg', action='store_true', help='generate pseudo reference database now with the current configuration')
    prfg_parser.add_argument('-fq', dest='fastq', type=str, nargs='+', default=[], help='fastq reads to search')
    prfg_parser.add_argument('-e', dest='exclude', type=str, nargs='+', default=[], help='frags files containing reads which are excluded from search')
    prfg_parser.add_argument('-p', dest='primer_file', type=str,
            help='''csv file containing forward and reverse primers used.
    This must have at least columns [fwd_id, fwd_seq, rev_id, rev_seq]''')
    prfg_parser.add_argument('-o', dest='db_file', type=str, help='output csv file of pseudo reference sequences')
    prfg_parser.add_argument('--low_mem', dest='low_mem', action='store_true', help='enable optimizations that reduce RAM used')
    
    # config for repeat frag finder
    fgs_parser = subparser.add_parser('fgs', help='suboptions for repeat fragment finder')
    config['fgs_config']='-k10 -w1'
    fgs_parser.add_argument('-c', dest='fgs_config', type=str, help='config options passed to minimap2')
    fgs_parser.add_argument('-fq', dest='fastq', type=str, nargs='+', default=[], help='fastq reads to search')
    fgs_parser.add_argument('-e', dest='exclude', type=str, nargs='+', default=[], help='frags files containing reads which are excluded from search')
    fgs_parser.add_argument('-db', dest='db_file', type=str, help='file containing reference sequences used in fragment search')
    config['frag_folder']='./frags/'
    fgs_parser.add_argument('-o', dest='frag_folder', type=str, help='folder where frags csv files are stored')
    fgs_parser.add_argument('-r', dest='run_fgs', action='store_true', help='run repeat fragment finder')
    fgs_parser.add_argument('-s', dest='config_file', type=str, help='write settings to configuration file')
    fgs_parser.add_argument('--low_mem', dest='low_mem', action='store_true', help='enable optimizations that reduce RAM used')

    # config for multi-sequence alignment
    msa_parser = subparser.add_parser('msa', help='suboptions for multi-sequence alignment')
    config['msa_config']='-n -15 -g -10 -l 0 -r 2'
    msa_parser.add_argument('-c', dest='msa_config', type=str, help='config options passed to spoa multi-sequence aligner')
    msa_parser.add_argument('-fq', dest='fastq', nargs='+', type=str, help='fastq files')
    msa_parser.add_argument('-i', dest='frag_folder', type=str, help='folder where frags csv files are stored')
    config['msa_folder']='./msa/'
    msa_parser.add_argument('-o1', dest='msa_folder', type=str, help='folder where msa input and output files are saved')
    msa_parser.add_argument('-o2', dest='cons_file', type=str, help='output csv file for consensus sequences')
    msa_parser.add_argument('-r1', dest='run_msa', action='store_true', help='run multi-sequence alignment')
    msa_parser.add_argument('-r2', dest='run_cons', action='store_true', help='read consensus after multi-sequence alignment')
    msa_parser.add_argument('-s', dest='config_file', type=str, help='write settings to configuration file')
    msa_parser.add_argument('--low_mem', dest='low_mem', action='store_true', help='enable optimizations that reduce RAM used')

    # config for matching primers to reads
    fpmr_parser = subparser.add_parser('fpmr', help='suboptions for matching primers to consensus reads')
    config['fpmr_config']='-k5 -w1 -s 20 -P'
    fpmr_parser.add_argument('-c', dest='fpmr_config', help='config options passed to minimap2 aligner')
    fpmr_parser.add_argument('-i', dest='cons_file', help='input untrimmed consensus csv files')
    fpmr_parser.add_argument('-p', dest='primer_file', type=str,
            help='''csv file containing forward and reverse primers used.
    This must have at least columns [fwd_id, fwd_seq, rev_id, rev_seq]''')
    fpmr_parser.add_argument('-o', dest='pmatch_file', help='csv file containing primer match information')
    fpmr_parser.add_argument('-s', dest='config_file', type=str, help='write settings to configuration file')
    fpmr_parser.add_argument('-r1', dest='run_fpmr', action='store_true', help='run primer matching')
    fpmr_parser.add_argument('-r2', dest='run_trmc', action='store_true', help='trim primers from reads')

    # config for clustering
    clst_parser = subparser.add_parser('clst', help='suboptions for clustering')
    clst_parser.add_argument('-i', dest='cin_file', type=str, help='input csv, fasta, or fastq file of sequences to cluster')
    clst_parser.add_argument('-o', dest='cout_file', type=str, help='output csv file for cluster center sequences')
    config['clst_init']=''
    clst_parser.add_argument('-init', dest='clst_init', type=str, help='csv, fasta, or fastq file containing initial cluster centers.')
    config['clst_folder']='./clusters/'
    clst_parser.add_argument('-f', dest='clst_folder', type=int, help='folder where clustering work is saved')
    config['clst_min_k']=5
    clst_parser.add_argument('-k', dest='clst_min_k', type=int, help='min_cluster_size')
    config['clst_csize']=20
    clst_parser.add_argument('-cs', dest='clst_csize', type=int, help='number of sequences from cluster center to multi-align into center sequence')
    config['clst_th_m']=0.9
    clst_parser.add_argument('-tm', dest='clst_th_m', type=float, help='threshold for marking clusters to merge')
    config['clst_th_s']=4
    clst_parser.add_argument('-ts', dest='clst_th_s', type=float, help='how many partition to split the sequences for sweep')
    config['clst_N']=2000
    clst_parser.add_argument('-N', dest='clst_N', type=int, help='size of sequence subsample')
    config['clst_N_iter']=10
    clst_parser.add_argument('-iter', dest='clst_N_iter', type=int, help='number of iterations to run clustering')
    config['clst_iter_out'] = 'clst'
    clst_parser.add_argument('-tout', dest='clst_iter_out', type=str, help='prefix of file to save cluster centers after each iteration')
    config['clst_track'] = False
    clst_parser.add_argument('--track', dest='clst_track', action='store_true', help='keep cluster centers after each iteration')
    config['clst_pw_config']='-k15 -w10 -p 0.9 -D'
    clst_parser.add_argument('-pw_config', dest='clst_pw_config', type=str, help='config passed to minimap2')
    clst_parser.add_argument('-s', dest='config_file', type=str, help='write settings to configuration file')
    clst_parser.add_argument('-r', dest='run_clst', action='store_true', help='run clustering')
    args = parser.parse_args()
    bpy.init_log(config['log_file'])
    if args.subcommand==None:
        logging.error('please select subcommand to run')
        sys.exit(1)

    # apply new settings
    config = update_config(config, vars(args))

    # default runs
    keys = ['prfg','fgs','msa','cons','fpmr','trmc','clst']
    for k in keys:
        config['run_'+k]=False # set everything to false

    # if nothing is selected, then default is to run everything
    if args.subcommand=='run':
        if type(args.run)==type(None):
            if args.db_file==None:
                config['run_prfg']=True
            else:
                config['run_prfg']=False
            keys = ['fgs','msa','cons','fpmr','trmc','clst']
            for k in keys:
                config['run_'+k]=True
        elif type(args.run)==str:
            args.run = args.run.split(',')
            keys = ['prfg','fgs','msa','cons','fpmr','trmc','clst']
            for k in keys:
                config['run_'+k] = k in args.run
    elif args.subcommand in keys:
        user_config = vars(args)
        for k in user_config.keys():
            if 'run_' in k:
                config[k] = user_config[k]

    # Print run settings
    logging.info('Running the pipeline')
    logging.info('config_file = '+config['config_file'])
    logging.info('primer_file = '+config['primer_file'])
    logging.info('pseudo_fs   = '+config['prfg_fs'])
    logging.info('prfg_config = '+config['prfg_config'])
    logging.info('db_file     = '+config['db_file'])
    logging.info('N fastq     = '+str(len(config['fastq'])))
    logging.info('frag_folder = '+config['frag_folder'])
    logging.info('fgs_config  = '+config['fgs_config'])
    logging.info('msa_folder  = '+config['msa_folder'])
    logging.info('msa_config  = '+config['msa_config'])
    logging.info('cons_file        = '+config['cons_file'])
    logging.info('fpmr_config      = '+config['fpmr_config'])
    logging.info('pmatch_file      = '+config['pmatch_file'])
    logging.info('cin_file         = '+config['cin_file'])
    logging.info('cout_file        = '+config['cout_file'])
    logging.info('clst_init        = '+config['clst_init'])
    logging.info('clst_folder      = '+config['clst_folder'])
    logging.info('clst_pw_config   = '+config['clst_pw_config'])
    logging.info('clst_N      = '+str(config['clst_N']))
    logging.info('clst_th_s   = '+str(config['clst_th_s']))
    logging.info('clst_th_m   = '+str(config['clst_th_m']))
    logging.info('clst_min_k  = '+str(config['clst_min_k']))
    logging.info('clst_csize  = '+str(config['clst_csize']))
    logging.info('clst_N_iter = '+str(config['clst_N_iter']))
    logging.info('workspace        = '+config['workspace'])
    logging.info('Run Pseudoref generator = '+str(config['run_prfg']))
    logging.info('Run Fragment Search     = '+str(config['run_fgs']))
    logging.info('Run Multi-Seq Alignment = '+str(config['run_msa']))
    logging.info('Run Consensus           = '+str(config['run_cons']))
    logging.info('Run primer match        = '+str(config['run_fpmr']))
    logging.info('Run trim primers        = '+str(config['run_trmc']))
    logging.info('Run clustering          = '+str(config['run_clst']))
    logging.info('low_mem                 = '+str(config['low_mem']))

    # Load the fastq data
    if config['run_fgs'] or config['run_msa'] or config['run_prfg']:
        logging.info('loading fastq data')
        if len(config['fastq'])==0:
            logging.error('no fastq files provided')
            sys.exit(1)
        reads = bpy.load_ONT_fastq(config['fastq'], low_mem=config['low_mem'])
        # exclude reads which already have hits
        if len(config['exclude']) > 0:
            L1 = len(reads)
            exclude = pd.concat([pd.read_csv(f) for f in config['exclude']])
            reads = reads[reads['id'].isin(exclude['query_id'])==False]
            del exclude
            logging.info('searching subset '+str(len(reads))+'/'+str(L1))
    
    # run pseudo reference database generator if selected
    if config['run_prfg']:
        # load primers
        primers = bpy.load_primers(config['primer_file'])
        # generate pseudo reference database using primer info
        logging.info('Generating pseudo reference database')
        ref = generate_pseudo_refdb(primers, reads, block=10000, fs=config['prfg_fs'], config=config['prfg_config'])
        ref.to_csv(config['db_file'], index=False, compression='infer')
    
    # load reference database
    if config['run_fgs']:
        ref = bpy.load_file(config['db_file'])

    if config['run_fgs']:
        # search or RCA fragments in reads
        logging.info('Searching for RCA frags in reads')
        find_RCA_frags(ref, reads, block=10000, output=config['frag_folder'], config=config['fgs_config'])
        logging.info('Aligner done')

    if config['run_msa']:
        # Loading the data for multi-sequence alignment
        logging.info('Loading frags from '+config['frag_folder'])
        files = glob.glob(config['frag_folder']+'fraginfo_*.csv.gz')
        frags = pd.concat([pd.read_csv(f) for f in files])
        # filtering out overlapping frags
        frags = bpy.remove_overlaps(frags, metric='AS', thresh=50)
        
        # annotating information about the reads
        logging.info('annotating reads')
        df_frags = annotate_reads(frags)
        frags = frags.merge(df_frags[['query_id','N frags']], on='query_id', how='left')
        logging.info('annotation complete')
        
        # Run multi-sequence alignment on fragments
        logging.info('Running multi-sequence aligner')
        perform_MSA(reads, frags[frags['N frags'] > 1], batch_size=100, config=config['msa_config'], thread_lock=False, folder=config['msa_folder'])
        logging.info('multi-sequence alignment done')
    
    # Generate consensus sequence and save the info
    if config['run_cons']:
        MSA_outfile = glob.glob(config['msa_folder']+'*.out')
        logging.info('Number of MSA_outfiles='+str(len(MSA_outfile)))        
        df = get_spoa_consensus(MSA_outfile)

        # Write the data to file
        logging.info('Writing consensus output to: '+config['cons_file'])
        df.to_csv(config['cons_file'], index=False, compression='infer')
        logging.info('Consensus generation done')
        
    # match primers to read consensus
    if config['run_fpmr']:
        # load primer info
        primers = bpy.load_primers(config['primer_file'])

        # load consensus info
        logging.info('Loading cons_file: '+config['cons_file'])
        df = pd.read_csv(config['cons_file'])
        
        # match the primers to consensus reads        
        logging.info('Matching primers to reads')
        df = df.rename(columns={'consensus':'sequence'})
        df = match_primer_to_read(primers, df[['id','sequence']], config=config['fpmr_config'])
        if len(df) == 0:
            logging.error('Fwd and Rev primers were not found in reads')
            sys.exit(1)

        # Write the data to file
        logging.info('Writing primer_match to '+config['pmatch_file'])
        df.to_csv(config['pmatch_file'], index=False, compression='infer')
        logging.info('primer matching done')

    if config['run_trmc']:
        # load data
        logging.info('Loading pmatch_file: '+config['pmatch_file'])
        df_pmatch = pd.read_csv(config['pmatch_file'])
        logging.info('Loading cons_file: '+config['cons_file'])
        df_cons = pd.read_csv(config['cons_file'])
        df_cons = df_cons.rename(columns={'consensus':'sequence'})

        # run the trimming operation
        logging.info('Trimming the reads')
        df = trim_consensus(df_pmatch, df_cons)

        # Write the data to file
        logging.info('Writing trimmed consensus sequence')
        df.to_csv(config['cin_file'], index=False, compression='infer')
        logging.info('primer matching done')
    
    if config['run_clst']:
        # load data
        logging.info('Loading cin_file: '+config['cin_file'])
        df = bpy.load_file(config['cin_file'])
        df = df.rename(columns={'consensus':'sequence'})
        df['id'] = ['read'+str(i) for i in range(0,len(df))]
        # only work on sequences which have both forward and reverse primers
        if 'fwd_primer_seq' in df.columns and 'rev_primer_seq' in df.columns:
            df = df[df['fwd_primer_seq'].notna() & df['rev_primer_seq'].notna()]
        # load initial cluster centers if given
        df_d = []
        if config['clst_init']!='':
            df_d = bpy.load_file(config['clst_init'])

        # refine the clusters
        data, flist = perform_cluster(df, df_d, max_iter=config['clst_N_iter'], csize=config['clst_csize'], N=config['clst_N'], th_s=config['clst_th_s'], th_m=config['clst_th_m'], pw_config=config['clst_pw_config'], msa_config=config['msa_config'], workspace=config['clst_folder'], track_file=config['clst_iter_out'], timestamp=True)

        # clean up if needed
        if config['clst_track']==False:
            bpy.batch_file_remove(flist)

        # save data
        logging.info('Writing clustered sequences to '+config['cout_file'])
        data.to_csv(config['cout_file'], index=False, compression='infer')
        logging.info('clustering done')
        
if __name__ == "__main__":
    main()
