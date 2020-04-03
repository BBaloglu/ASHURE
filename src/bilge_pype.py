#!/usr/bin/env python
# Helper tools for sequence data
# bilgenurb@gmail.com

# Loading the libraries
import numpy as np
import pandas as pd

# signal processing libraries
import scipy as sp
import scipy.signal

# clustering libraries
import sklearn
import sklearn.cluster
import sklearn.decomposition
import hdbscan

# For interfacing with the file system
import subprocess
import os
import time
import datetime
import glob
import sys
import logging

def append_to_log(text, logfile='log.txt'):
    '''
    write text to a log file and also print out what was written
    '''
    timestamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S.%f ')
    
    f = open(logfile,'a')
    f.write(timestamp + ' pid['+str(os.getpid())+'] ' + text+'\n')
    f.close()
    print(text)

def init_log(fname='ashure.log', level='DEBUG'):
    fmt = '%(asctime)s.%(msecs)03d %(levelname)s: %(message)s'
    dfmt = '%Y-%m-%d %H:%M:%S' 
    if level=='WARNING':
        logging.basicConfig(filename=fname, format=fmt, datefmt=dfmt, level=logging.WARNING)
    elif level=='INFO':
        logging.basicConfig(filename=fname, format=fmt, datefmt=dfmt, level=logging.INFO)
    elif level=='DEBUG':
        logging.basicConfig(filename=fname, format=fmt, datefmt=dfmt, level=logging.DEBUG)

def get_ONT_header(text):
    text = text.split(' ')
    readid = text[0][1:]
    for i in range(0,len(text)):
        line = text[i].split('=')
        if line[0]=='read':
            read = int(line[1])
            line = text[i+1].split('=')
            chan = int(line[1])
            line = text[i+2].split('=')
            datetime = line[1].split('T')
            day = datetime[0]
            time = datetime[1].split('Z')[0]
            return [readid, read, chan, day+' '+time]
    return -1 # catch errors
    
def read_ONT_fastq(file):
    '''
    Function to parse fastq file and generate a pandas dataframe of the following info
    
    file = fastq file to parse
    output = pandas dataframe with columns = [file, id, read, channel, datetime, length, sequence, quality]
    '''    
    
    # open file and read it
    f = open(file,'r')
    text = f.read()
    f.close()
    
    # split lines and get relevant info
    text = text.split('\n')
    header = text[0::4] # header is the first line
    seq = text[1::4] # sequence is second line
    quality = text[3::4] # quality score is third lines
    header = header[:len(seq)]
    data = []
    for i in range(0,len(header)):
        h = get_ONT_header(header[i])
        data.append([file,h[0],h[1],h[2],h[3],len(seq[i]),seq[i],quality[i]])
    # return as pandas dataframe table
    col = ['file','id','read','chan','datetime','length','sequence','quality']
    return pd.DataFrame(data, columns=col)

# Function to load fastq data and primer info
def load_basecalled_data(files, log_file='log.txt'):
    '''
    Function to load fastq info from list of ONT fastq files
    files = list of fastq files
    output = pandas dataframe of ONT fastq info
    '''
    # Generate table of sequences via iteration
    start = time.time()
    df = []
    for i in range(0,len(files)):
        df.append(read_ONT_fastq(files[i])) # Save tables as a vector
        if i%20==0: # print out progress
            append_to_log('Processing '+str(i)+'/'+str(len(files))+', current file = '+files[i], log_file)
            append_to_log('elapse = {0:.2f}'.format(time.time()-start)+'s', log_file)
    df = pd.concat(df) # Merge all the tables
    df['datetime'] = pd.to_datetime(df['datetime']) # Format to datetime
    return df

def load_primers(primer_file):
    '''
    Reads primer information from a csv file and checks it has the expected structure
    primer csv must have the following columns:
    [fwd_id,fwd_sequence,rev_id,rev_sequence]

    fwd_id  = name of forward primer
    fwd_seq = sequence of forward primer
    rev_id = name of reverse primer
    rev_seq = sequence of reverse primer

    output = pandas dataframe of primer information
    '''
    print('Reading primers from ',primer_file)
    df = pd.read_csv(primer_file, delimiter = ',')
    # safety checks on the primer data
    L1 = len(df)
    df = df.dropna()
    L2 = len(df)
    if L1!=L2:
        print('Warning: nan info was dropped from primer data')
    
    # check basic info is there
    col = ['fwd_id','fwd_seq','rev_id','rev_seq']
    for c in col:
        if (c in df.columns) == False:
            print('column '+c+' not found in ',primer_file)
            sys.exit(1)
    return df
        
def read_fastq(file, low_mem=False):
    '''
    Reads list of sequences from a fastq file format.
    
    file = file to read the information
    returns list of names, sequence, quality
    '''
    # open file and read it
    f = open(file,'r')
    text = f.read()
    f.close()
    
    # split lines and get relevant info
    text = text.split('\n')
    header = text[0::4] # header is the first line
    sequence = text[1::4] # sequence is second line
    quality = text[3::4] # quality score is third lines
    header = header[:len(sequence)]
    data = []
    for i in range(0,len(header)):
        data.append([header[i],sequence[i],quality[i]])
    return pd.DataFrame(data, columns=['id','sequence','quality'])

def write_fastq(file, data):
    '''
    Write list of sequences into a fastq file format.
    
    data =  array with [id, sequence, quality] as the columns
    file = file to write the information
    '''
    f = open(file,'w')
    if len(data.shape) > 1: # case for table
        for i in data:
            text = '@' + str(i[0]) + '\n'
            text+= str(i[1]) + '\n'
            text+= '+\n'
            text+= str(i[2]) + '\n'
            f.write(text)
    else: # case for single row
        text = '@' + data[0] + '\n'
        text+= str(data[1]) + '\n'
        text+= '+\n'
        text+= str(data[2]) + '\n'
        f.write(text)
    f.close()

def read_fasta(file):
    '''
    Reads list of sequences from a fasta file format.
    
    file = file to read the information
    returns list of sequence and names
    '''
    f = open(file,'r')
    text = f.read().split('>')[1:]
    f.close()
    data = []
    for i in text:
        i = i.split('\n')
        data.append([i[0],''.join(i[1:])])
    out = pd.DataFrame(data, columns=['id','sequence'])
    return out

def write_fasta(file, data):
    '''
    Write list of sequences into a fasta file format.
    
    data =  array with [id, sequence] as the columns
    file = file to write the information
    '''
    f = open(file,'w')
    if len(data.shape) > 1: # case for table
        for i in data:
            text = '>' + str(i[0]) + '\n'
            text+= str(i[1]) + '\n'
            f.write(text)
    else: # case for single row
        text = '>' + str(data[0]) + '\n'
        text+= str(data[1]) + '\n'
        f.write(text)
    f.close()

def get_SAMinfo(read, key): # debug rewrite sam file parser to be more efficient
    '''
    Function to parse useful info from a line of the SAM file
    '''
    read = read.split('\t')
    qname = read[0]
    flag = '{0:10b}'.format(int(read[1])) # decode the sam flag via binary to see if reverse complement is aligned
    rname = read[2]
    if flag[-5] == '1':
        orientation = '-'
        pos2 = int(read[3])
        pos1 = pos2 + int(read[8])
    else:
        orientation = '+'
        pos1 = int(read[3])
        pos2 = pos1 + int(read[8])
    mapq = int(read[4])
    cigar = read[5]
    # get another sam info
    info = [0]*len(key)
    if rname != '*':
        for i in range(11,len(read)):
            x = read[i].split(':') # parse it backwards
            if x[0] in key.keys():
                info[key[x[0]]] = x[-1]
    data = [qname, rname, orientation, flag, pos1, pos2, mapq, cigar]
    return data, info

def parse_SAMfile(file):
    '''
    This function parses the information from a SAM file and returns a dataframe   
    '''
    f = open(file)
    text = f.read()
    f.close()

    # Removing the header
    text = text.split('\n@PG')[-1].split('\n')[1:-1]
    data = []
    i = 0
    col1 = ['query_id','database_id','orientation','flag','t_start','t_end','mapq','CIGAR']
    col2 = ['AS','XS','XN','XM','XO','XG','NM']
    key = {col2[i]:i for i in range(0,len(col2))}
    while i < len(text):
        out, info = get_SAMinfo(text[i], key)
        if out[3][-2] == '1':
            # if it is a paired read, look for its pair on the next line
            i = i+1
            out2, info2 = get_SAMinfo(text[i], key)
            out[0]+= ' '+out2[0]
            out[7]+= ' '+out2[7]
            for j in range(0,len(info)):
                info[j]+=info2[j]
        data.append(out+info)
        # increment the counter
        i = i+1
    return pd.DataFrame(data, columns=col1+col2)

def run_bowtie2(query, database, workspace='./bowtie2/',
                config='-a --very-sensitive-local --threads 1 --quiet',
                build_index=True, cleanup=False, log_file='log.txt'):
    '''
    This is a wrapper for the bowtie2 aligner.

    query = list of sequences to search for in database
    query must be in pandas dataframe format with columns:
    [id, sequence] or [id, sequence, quality]

    database = database of sequences being search through
    database must be in pandas dataframe format with columns:
    [id, sequence] or [id, sequence, quality]
    [fwd_id, fwd_seq, rev_id, rev_seq] or
    [fwd_id, fwd_seq, fwd_quality, rev_id, rev_seq, rev_quality]

    build_index = determines if a burrow-wheeler transform index should be build for the database
    workspace = defines where input and output files for bowtie2 resides
    config = defines the config options to pass to bowtie2
    '''
    # Define locations of input and output files for bowtie
    qfile = workspace+'fwdp.fa'
    mfile = workspace+'revp.fa'
    dbfile = workspace+'database.fa'
    outfile = workspace+'results.sam'
    btfile = workspace+'index'

    check_dir(workspace)

    # Only build index when we need to. Index building takes a lot of time
    if build_index == True:
        value = check_seqdf(database, log_file)
        if value == 2:
            write_fasta(dbfile, database[['id','sequence','quality']].values)
            cmd = 'bowtie2-build --quiet -q '+dbfile+' '+btfile+' --threads 2'
        elif value == 1:
            write_fasta(dbfile, database[['id','sequence']].values)
            cmd = 'bowtie2-build --quiet -f '+dbfile+' '+btfile+' --threads 2'
        else:
            return -1
        subprocess.run(cmd.split(' '))

    # Write input files
    value = check_seqdf(database, log_file)
    if value == 4:
        write_fastq(qfile, query[['fwd_id','fwd_seq','fwd_quality']].values)
        write_fastq(mfile, query[['rev_id','rev_seq','fwd_quality']].values)
        cmd = 'bowtie2 '+config+' -x '+btfile+' -q -1 '+qfile+' -2 '+mfile+' -S '+outfile
    elif value == 3:
        write_fasta(qfile, query[['fwd_id','fwd_seq']].values)
        write_fasta(mfile, query[['rev_id','rev_seq']].values)
        cmd = 'bowtie2 '+config+' -x '+btfile+' -f -1 '+qfile+' -2 '+mfile+' -S '+outfile
    elif value == 2: # write fastq
        write_fastq(qfile, query[['id','sequence','quality']].values)
        cmd = 'bowtie2 '+config+' -x '+btfile+' -q -U '+qfile+' -S '+outfile
    elif value == 1: # write fasta
        write_fasta(qfile, query[['id','sequence']].values)
        cmd = 'bowtie2 '+config+' -x '+btfile+' -f -U '+qfile+' -S '+outfile
    else:
        return -1

    # call bowtie2
    subprocess.run(cmd.split(' '))
    data = parse_SAMfile(outfile)

    # add q_len and t_len info
    q_len = np.transpose([query['id'].values, [len(i) for i in query['sequence']]])
    q_len = pd.DataFrame(q_len, columns = ['query_id','q_len'])
    data = data.merge(q_len, on='query_id', how='left')

    t_len = np.transpose([database['id'].values, [len(i) for i in database['sequence']]])
    t_len = pd.DataFrame(t_len, columns = ['database_id','t_len'])
    data = data.merge(t_len, on='database_id', how='left')
    data = data.dropna() # drop shit that doesn't have sequence length

    # Format the fields to integer
    x = ['q_len','t_len','t_start','t_end','mapq','AS','XS','XN','XM','XO','XG','NM']
    for i in x:
        data[i] = data[i].astype(int)

    # compute similarity
    v1 = np.min([data['q_len'].values, data['t_len'].values], axis = 0)
    data['similarity'] = 1-data['NM']/v1
        
    # remove work folder
    if cleanup:
        subprocess.run(['rm', '-r', workspace])
    return data

def run_bwa(query, database, workspace='./bwa/', config=' mem -a ',
            build_index=True, log_file='log.txt', cleanup=False):
    '''
    This is a wrapper for the bwa aligner.

    query = list of sequences to search for in database
    query must be in pandas dataframe format with columns:
    [id, sequence] or [id, sequence, quality]
    [fwd_id, fwd_seq, rev_id, rev_seq] or
    [fwd_id, fwd_seq, fwd_quality, rev_id, rev_seq, rev_quality]

    database = database of sequences being search through
    database must be in pandas dataframe format with columns:
    [id, sequence] or [id, sequence, quality]

    build_index = determines if a burrow-wheeler transform index should be build for the database
    workspace = defines where input and output files for bowtie2 resides
    config = defines the config options to pass to bowtie2
    '''
    # Define locations of input and output files for bowtie
    qfile = workspace+'read1.fq'
    mfile = workspace+'read2.fq'
    dbfile = workspace+'database.fa'
    outfile = workspace+'results.sam'

    # Make the working workspace if it does not exist
    check_dir(workspace)

    # Only build index when we need to. Index building takes a lot of time
    if build_index:
        value = check_seqdf(database, log_file)
        if value == 2:
            write_fasta(dbfile, database[['id','sequence','quality']].values)
        elif value == 1:
            write_fasta(dbfile, database[['id','sequence']].values)
        else:
            return -1
        cmd = 'bwa index '+dbfile
        subprocess.run(cmd.split(' '))

    # Check if we are aligning paired reads
    value = check_seqdf(query, log_file)
    
    paired = False
    if value == 4:
        write_fastq(qfile, query[['fwd_id','fwd_seq','fwd_quality']].values)
        write_fastq(mfile, query[['rev_id','rev_seq','rev_quality']].values)
        paired = True
    elif value == 3:
        write_fasta(qfile, query[['fwd_id','fwd_seq']].values)
        write_fasta(mfile, query[['rev_id','rev_seq']].values)
        paired = True
    elif value == 2:
        write_fastq(qfile, query[['id','sequence','quality']].values)
    elif value == 1:
        write_fasta(qfile, query[['id','sequence']].values)
    else:
        return -1
    
    if paired:
        cmd = 'bwa '+config+' '+dbfile+' '+qfile+' '+mfile+' -o '+outfile
    else:
        cmd = 'bwa '+config+' '+dbfile+' '+qfile+' -o '+outfile

    # call the aligner
    subprocess.run(cmd.split(' '))

    # extract the SAM file information
    data = parse_SAMfile(outfile)

    # add q_len and t_len info
    q_len = np.transpose([query['id'].values, [len(i) for i in query['sequence']]])
    q_len = pd.DataFrame(q_len, columns = ['query_id','q_len'])
    data = data.merge(q_len, on='query_id', how='left')

    t_len = np.transpose([database['id'].values, [len(i) for i in database['sequence']]])
    t_len = pd.DataFrame(t_len, columns = ['database_id','t_len'])
    data = data.merge(t_len, on='database_id', how='left')
    data = data.dropna() # drop shit that doesn't have sequence length

    # Format the fields to integer
    x = ['q_len','t_len','t_start','t_end','mapq','AS','XS','XN','XM','XO','XG','NM']
    for i in x:
        data[i] = data[i].astype(int)

    # compute similarity
    v1 = np.min([data['q_len'].values, data['t_len'].values], axis = 0)
    data['similarity'] = 1-data['NM']/v1

    # remove work folder
    if cleanup:
        subprocess.run(['rm','-r',workspace])
    return data

def parse_PAFfile(file):
    '''
    This function parses the information from a PAF file and returns a dataframe
    '''
    # parse the sam file
    f = open(file)
    text = f.read().split('\n')
    f.close()
    
    if len(text[-1])==0:
        text = text[:-1] # removes the last line that is empty    
    data = []
    
    # information we are parsing
    col1 = ['query_id','q_len','q_start','q_end','orientation',
            'database_id','t_len','t_start','t_end','match','tot','mapq']
    col2 = ['tp','cm','s1','s2','NM','AS','ms','nn','rl','cg']
    key = {col2[i]:i for i in range(0,len(col2))}
    for line in text:
        out = line.split('\t')
        # parses for other info
        info = [0]*len(key)
        for i in range(len(col1),len(out)):
            x = out[i].split(':')
            if x[0] in key:
                info[key[x[0]]] = x[-1]
        # save the data
        data.append(out[:len(col1)]+info)
    data = pd.DataFrame(data, columns=col1+col2)
    data = data.rename(columns = {'cg':'CIGAR'})
    return data

def run_minimap2(query, database, workspace = './minimap2/', config = '-x map-ont', cigar = True,
                 build_index = True, use_index = False, cleanup = False, log_file = 'log.txt'):
    '''
    This is a wrapper for the minimap2 aligner. This aligner is better suited for long reads (>100bp)

    query = list of sequences to search for in database
    query must be in pandas dataframe format with columns:
    [id, sequence] or [id, sequence, quality]
    [fwd_id, fwd_seq, rev_id, rev_seq] or
    [fwd_id, fwd_seq, fwd_quality, rev_id, rev_seq, rev_quality]

    database = database of sequences being search through
    database must be in pandas dataframe format with columns:
    [id, sequence] or [id, sequence, quality]

    build_index = determines if a burrow-wheeler transform index should be build for the database
    use_index = determines if index will be used or just read from database fasta file
    workspace = defines where input and output files for bowtie2 resides
    config = defines the config options to pass to bowtie2
    '''
    # Define locations of input and output files for bowtie
    qfile = workspace+'read1.fq'
    mfile = workspace+'read2.fq'
    dbfile = workspace+'database.fq'
    outfile = workspace+'results.paf'
    btfile = workspace+'index.mmi'

    check_dir(workspace)

    # Write database to fasta or fastq
    value = check_seqdf(database, log_file)
    if value == 2:
        write_fastq(dbfile, database[['id','sequence','quality']].values)
    elif value == 1:
        write_fasta(dbfile, database[['id','sequence']].values)
    else:
        return -1

    # Only build index when we need to. Index building takes a lot of time
    if build_index:
        cmd = 'minimap2 '+config+' -d '+btfile+' '+dbfile
        subprocess.run(cmd.split(' '))
        
    # Write query to fasta or fastq
    value = check_seqdf(query, log_file)
    paired = False
    if value == 4:
        write_fastq(qfile, query[['fwd_id','fwd_seq','fwd_quality']].values)
        write_fastq(mfile, query[['rev_id','rev_seq','rev_quality']].values)
        paired = True
    elif value == 3:
        write_fasta(qfile, query[['fwd_id','fwd_seq']].values)
        write_fasta(mfile, query[['rev_id','rev_seq']].values)
        paired = True
    elif value == 2:
        write_fastq(qfile, query[['id','sequence','quality']].values)
    elif value == 1:
        write_fasta(qfile, query[['id','sequence']].values)
    else:
        return -1
    
    if cigar: # adds flag to generate cigar string
        config = '-c '+config

    if use_index:
        dbfile = btfile
    if paired:
        cmd = 'minimap2 '+config+' '+dbfile+' '+qfile+' '+mfile+' -o '+outfile
    else:
        cmd = 'minimap2 '+config+' '+dbfile+' '+qfile+' -o '+outfile

    # call the aligner
    subprocess.run(cmd.split(' '))
    
    # extract the PAF file information
    data = parse_PAFfile(outfile)
    
    # Format the fields to integer
    x = ['q_len','q_start','q_end','t_len','t_start','t_end','match','tot','mapq',
         'cm','s1','s2','NM','AS','ms','nn','rl']
    for i in x:
        data[i] = data[i].astype(int)
    data['match_score'] = data['match']/data['tot']*1.0 # compute blast like score
    
    # compute similarity
    v1 = np.min([data['q_len'].values, data['t_len'].values], axis = 0)
    data['similarity'] = data['match']/v1
    
    # drop cigar if it is not present
    if cigar == False:
        data.drop(columns = ['CIGAR'], inplace = True)
    
    # remove work folder
    if cleanup and ~use_index and ~build_index:
        subprocess.run(['rm','-r',workspace])
    return data

def check_seqdf(df, log_file = 'log.txt'):
    '''
    This function looks at the header of the dataframe and checks if it is fasta, fastq, or paired fasta or fastq
    return 1 if fasta
    return 2 if fastq
    return 3 if paired fasta
    return 4 if paired fastq
    return -1 if nothing matches
    '''
    cols = df.columns
    
    s1 = 'id' in cols
    s2 = 'sequence' in cols
    s3 = 'quality' in cols
    s4 = 'fwd_id' in cols
    s5 = 'fwd_seq' in cols
    s6 = 'fwd_quality' in cols
    s7 = 'rev_id' in cols
    s8 = 'rev_seq' in cols
    s9 = 'rev_quality' in cols
    
    # if it is fasta dataframe with columns
    if s1 & s2 & ~s3:
        return 1
    # if it is a fastq dataframe with columns
    elif s1 & s2 & s3:
        return 2
    # if it is a paired fasta dataframe
    elif s4 & s5 & ~s6 & s7 & s8 & ~s9 :
        return 3
    # if it is a paired fastq dataframe
    elif s4 & s5 & s6 & s7 & s8 & s9 :
        return 4
    # if nothing matches, then print error message
    else:
        append_to_log('columns of data = '+''.join(cols), log_file)
        append_to_log('''
        Input dataframe does not have the proper header.
        
        The following are valid columns
        
        fasta data   = [id, sequence]
        fastq data   = [id, sequence, quality]
        paired fasta = [fwd_id, fwd_seq, rev_id, rev_seq]
        paired fastq = [fwd_id, fwd_seq, fwd_quality, rev_id, rev_seq, rev_quality]
        ''', log_file)
        return -1

def check_dir(folder, overwrite=True):
    '''
    This function checks if the directory exists. If it does not exist, make the directory
    '''
    if folder[-1]!='/': folder+='/'

    if os.path.exists(folder) == False and overwrite:
        print('Making directory ',folder)
        cmd = 'mkdir ' + folder
        subprocess.run(cmd.split(' '))
        return folder
    elif overwrite==False:
        i=0
        while os.path.exists(folder[:-1]+str(i)+'/'):
            i+=1
        return check_dir(folder[:-1]+str(i)+'/')

def find_overlaps(df, thresh = 1):
    '''
    Function that iterates through list and clusters together overlapping id
    df = dataframe containing columns [id, start, stop]
    thresh = allow some threshold nucleotide overlap before calling it an overlap    
    Returns dataframe with addition column = 'clusterN'    
    '''
    data = df.values
    s = np.lexsort((data[:,1],data[:,0]))
    data = data[s]
    c = np.zeros(len(data))
    
    # set initial values
    cur_id = ''
    cur_start = 0
    cur_stop = 0
    cluster = 0
    clusterN = []
    N=0
    
    # iterate through the data
    for i in range(0,len(data)):
        # edge of id reached or new start exceeds prev stop
        # might be buggy here --> look at this later #debug
        s1 = (cur_id!=data[i,0])
        s2 = (data[i,1] > cur_stop)
        s3 = (data[i,1] > cur_stop - thresh and thresh*2 < data[i][2] - data[i][1])
        if s1 or s2 or s3:
            # only run if edge of cluster is reached
            [cur_id, cur_start, cur_stop] = data[i][0:3]
            clusterN.append(N)
            cluster+=1
            N=0
        elif cur_stop < data[i,2]: # keep moving boundary of cluster forward
            cur_stop = data[i,2]
        
        c[i] = cluster
        N+=1
    clusterN.append(N) # end of list reached, append final value
    
    # add cluster num label to data
    data = df.iloc[s].copy()
    data['cluster_num'] = c

    # generate reference list of total counts for each cluster num label
    clusterN = np.transpose([range(0,len(clusterN)),clusterN])
    clusterN = pd.DataFrame(clusterN, columns = ['cluster_num','cluster_item_count'])
    
    # merge the data
    return data.merge(clusterN, on = 'cluster_num', how = 'left').set_index(df.index[s]) # merge and reset index

def split_overlaps(df, thresh = 1):
    '''
    Function that iterates through list and reduces it to best scoring non-overlapping id
    df = dataframe with columns [id, start, stop, score]
    thresh = allow some threshold nucleotides of overlap before calling it an overlap
    Returns cleaned items
    '''
    data = df.values # work on numpy array because its faster
    s = np.lexsort((data[:,3],data[:,0]))[::-1]
    data = data[s]
        
    # initial value for keep --> keep everything by default
    keep = np.array([True]*len(data))
    for i in range(0,len(data)):
        if keep[i]:
            for j in range(i+1,len(data)):
                # exit loop if end of read is reached
                if data[i,0]!=data[j,0]:
                    break
                else:
                    # q_stop < cur_start or q_start > cur_stop --> then frags are not overlapping
                    s1 = data[i,2] < data[j,1]+thresh
                    s2 = data[i,1] > data[j,2]-thresh
                    # f1 in f2
                    # f2 overlap f2 beyond thresh
                    # thresh not triggered, but frags too small
                    c = (s1 & s2) | ~(s1 | s2) | (thresh*2 > data[j,2] - data[j,1])
                    keep[j] = bool(~c)
    # decode the original ordering and return data
    data = df.iloc[s[keep]]
    return data

def remove_overlaps(frags, metric='AS', thresh=0, log_file='log.txt'):
    '''
    Function to filter out bad overlapping matches in sequences
    
    frags = output from find_RCA_frags()
            must be dataframe with columns [query_id, q_start, q_end, <alignment score>]
    thresh = amount of overlap allowed
    metric = <alignment score> defined by user, default is AS from bwa or minimap2
    
    returns same dataframe with overlapping fragments removed
    '''
    df = frags.reset_index(drop=True)
    
    append_to_log('filtering to best aligned fragments', log_file)
    f1 = find_overlaps(df[['query_id','q_start','q_end']], thresh = thresh)
    s = f1[f1['cluster_item_count']>1]
    # do work only on overlapping info
    if len(s) > 0:
        f2 = split_overlaps(df.iloc[s.index][['query_id','q_start','q_end',metric]], thresh = thresh)
        f2 = df.iloc[f2.index] # get the split up frags
        f1 = df.iloc[f1[f1['cluster_item_count']==1].index] # get single frags
        df = pd.concat([f1,f2]) # add them together
    # report results
    append_to_log('Cleaned out '+str(len(frags) - len(df))+' of '+str(len(frags)), log_file)
    return df

def get_best(df, col, metric='AS', stat='idxmax'):
    df=df.reset_index(drop=True)
    idx = df.groupby(by = col).agg({metric:stat}).reset_index()
    return df.iloc[idx[metric].values]

def run_msa(MSA_infiles, aligner='spoa', config='-l 0 -r 2', thread_lock=True):
    '''
    Wrapper for spoa or mafft multi-sequence aligner
    MSA_infile = list of files for input with .fq or .fa at end of filename
    aligner = aligner to use. Default = spoa
    config = configs to pass to aligner
    output is saved in same directory as input files except the file extension is now .out
    thread_lock = makes sure all spoa calls run to completion
                  switching to false and using system pipe can be faster,
                  but danger for process collisions
    '''
    processes = []
    if thread_lock:
        for msa_in in MSA_infiles:
            MSA_outfile = msa_in + '.out'
            f = open(MSA_outfile, 'w') # open file for writing std out
            cmd = aligner+' '+config+' '+msa_in
            processes.append([subprocess.Popen(cmd.split(' '), stdout = f), f]) # submit the job

        # wait until all jobs finish. Jobs have max time of 120s
        for j, f in processes:
            retval = j.wait(120)
            f.close() # close the file
            if retval!=0:
                print('syntax error on spoa call: ',j)
                return -1
    else:
        cmd = ''
        for msa_in in MSA_infiles:
            MSA_outfile = msa_in + '.out'
            f = open(MSA_outfile, 'w') # open file for writing std out
            cmd+= aligner+' '+config+' '+msa_in+' > '+MSA_outfile+' & '
        subprocess.call(cmd[:-3], shell = True) # submit the job

def read_spoa(fname):
    '''
    Parses output from spoa multi-sequence aligner
    '''
    f = open(fname,'r')
    text = f.read().split('\n')[:-1]
    f.close()

    consensus = ''
    msa = []
    if text[0].split(' ')[0] == 'Consensus':
        consensus = text[1]
        
    if text[2] == 'Multiple sequence alignment':
        msa = text[3:]
    elif text[0] == 'Multiple sequence alignment':
        msa = text[1:]
    return [consensus, msa]
    
def batch_file_remove(files):
    '''
    Function to call rm and batch remove a list of files
    '''
    # split file remove to batches of 100 files
    if len(files) > 0:
        split_files = np.array_split(files, np.ceil(len(files)/100))
        for files in split_files:
            cmd =['rm']
            for fname in files:
                cmd.append(fname)
            subprocess.run(cmd)
    
def get_qual_consensus(msa_out):
    # loads sequence info
    msa_in = msa_out.split('.out')[0]
    [cons, seq] = read_spoa(msa_out)
    
    # loads phred info and adds gaps
    x = read_fastq(msa_in)
    phred = x['quality'].values
    pgap = []
    for i in range(0,len(seq)):
        temp = phred_gap(seq[i],phred[i]) # add gaps to phred string 
        pgap.append(phred_ascii_to_int(temp)) # translate phred ascii to integer quality score
    pgap = np.transpose(pgap)
    
    lets = ['A','T','G','C','-']
    p1 = np.zeros([len(seq[0]),5])
    p2 = np.zeros([len(seq[0]),5])
    p3 = np.zeros([len(seq[0]),5])
    p4 = np.zeros([len(seq[0]),5])
    
    for j in range(0,len(seq)):
        s = np.array([k for k in seq[j]])
        for i in range(0,len(lets)):
            # compute prod of P(not A | A_i)
            p1[:,i]+=pgap[:,j]*(s==lets[i])
            # compute prod of P(A | A_i) = 1 - P(not A | A_i )
            p2[:,i]-=np.log10(1-10**(-pgap[:,j]/10))/10*(s==lets[i])
            # compute prod of P(not A | C_i) = 1 - P(not C | C_i)/4
            p3[:,i]-=np.log10(1-0.25*10**(-pgap[:,j]/10))/10*(s!=lets[i])
            # compute prod of P(A | C_i) = P(not C | C_i)/4
            p4[:,i]+=(pgap[:,j]-10*np.log10(0.25))*(s!=lets[i])

    # get P(A | AAACC)
    x = 10**(-(p1+p3)/10)
    y = 10**(-(p2+p4)/10)
    z = y/(x+y)
    
    # P(gap | AAACC) = 1 - P(A or C or G or T | AAACC)
    
    
    # find max and set as value
    idx = np.argmax(z, axis=1)
    out = ''.join([lets[i] for i in idx])
    phred = np.zeros(len(idx))
    for i in range(0,len(idx)):
        phred[i] = z[i,idx[i]]

    # count blank space
    blk = np.zeros(len(seq[0]))
    for i in seq:
        s = np.array([k for k in i])
        blk+= (s=='-')

    # count agreement
    agree = np.zeros(len(seq[0]))
    x = []
    for s in seq:
        x.append([i for i in s])
    x = np.array(x)
    
    for i in range(0,len(idx)):
        agree[i] = np.sum(x[:,i]==lets[idx[i]])
    
    avg_agree = np.mean(agree)
    std_agree = np.std(agree)
    avg_blk = np.mean(blk)
    std_blk = np.std(blk)
    
    phred = np.array(phred)
    error = np.sum(1-phred)
    avg = np.mean(1-phred)
    stdev = np.std(1-phred)

    phred = np.round(np.log10(1-phred)*-10)    
    phred[phred>93] = 93

    phred = [chr(int(i)+33) for i in phred]
    return out, ''.join(phred), avg, stdev, error, avg_agree, std_agree, avg_blk, std_blk

def phred_ascii_to_int(x):
    # Translate phred ASCII to quality integers
    return [ord(i)-33 for i in x]

def phred_gap(seq, x):
    # insert gaps into phred ASCII
    seq = seq.split('-')
    s1 = 0
    s2 = 0
    out = ''
    for i in seq:
        s2+=len(i)
        out+=x[s1:s2]+'"' # add a space character which evaluates to -1 quality
        s1=s2
    return out[:-1]    

def find_endgap(seq):
    s1 = -1
    s2 = -1
    for i in range(0,len(seq)):
        if seq[i]!='-' and s1 == -1:
            s1 = i
        if seq[len(seq)-i-1]!='-' and s2 == -1:
            s2 = len(seq)-i
        if s1 != -1 and s2 != -1:
            break
    return [s1,s2]

def get_pairwise_distance(seq, block=500, output_folder='./pairwise_dist/', workspace='./pw_aligner/',
                          config='-k15 -w10 -PD', log_file='log.txt', symmetric=False):
    '''
    Function to compute pairwise distance matrix of a list of sequences
    seq = pandas dataframe of a list of sequences with columns [id, sequence]
    block = number of sequences to compare in a block --> break up computation and save it
    output_folder = where output csv files are stored
    workspace = workspace folder for the aligner
    configs = alignment config options for minimap2
    log_file = log file to record progress
    symmetric = evaluates only pairwise alignments of a longer sequence 
                to a shorter one and returns data for a symmetric matrix
                symmetric evaluations take half time
    returns list of csv files where data is stored
    '''
    check_dir(output_folder)
    df = seq.copy()
    # sort it by longest sequence first
    if symmetric:
        L = np.array([len(s) for s in df['sequence']])
        s = np.argsort(L)[::-1]
        df = df.iloc[s]
    
    # run the aligner on blocks of the data
    df = np.array_split(df, np.ceil(len(df)/block))
    
    # figure out how many blocks we compute
    if symmetric:
        N = (len(df)+1)*len(df)/2
    else:
        N = len(df)*len(df)
    
    data = []
    count = 0
    for i in range(0,len(df)):
        js = 0
        if symmetric:
            js = i
        # only build index first time we use a new database
        build_index = True
        for j in range(js,len(df)):
            append_to_log('processing block '+str(i)+','+str(j)+' '+str(count)+'/'+str(N), log_file)
            count+=1
            out = run_minimap2(df[j], df[i], config = config, workspace = workspace, 
                               cigar = True, build_index = build_index, use_index = True)
            build_index = False
            fname = output_folder+'dst_'+str(i)+'_'+str(j)+'.csv.gz'
            out.to_csv(fname, index=False, compression='infer')
            data.append(fname)
    return data

def get_feature_vector(df, d_dct=None, q_dct=None, symmetric=False):
    '''
    Reorganizes array so that feature vectors are columns
    df = dataframe with columns [query_id, database_id, metric]
    d_dct = dictionary mapping database_id to column number
    q_dct = dictionary mapping query_id to row number
    symmetric = if d_dct and q_dct are not given, then make a symmetric matrix
                from unique values of query_id and database_id
    return dataframe with columns [d_i, ... , d_n, id]
    '''
    # create dict for features
    if symmetric and d_dct == None and q_dct == None:
        qlist = np.unique([i for i in df['query_id']]+[i for i in df['database_id']])
        q_dct = {qlist[i]:i for i in range(0,len(qlist))}
        d_dct = q_dct
    else:
        if d_dct==None:
            dlist = np.unique(df['database_id'])
            d_dct = {dlist[i]:i for i in range(0,len(dlist))}
        if q_dct==None:
            qlist = np.unique(df['query_id'])
            q_dct = {qlist[i]:i for i in range(0,len(qlist))}

    # initial variables for storage
    v = np.zeros([len(q_dct), len(d_dct)])
    for q, d, m in df.values:
        if (q in q_dct) and (d in d_dct):
            v[q_dct[q], d_dct[d]] = m
    data = pd.DataFrame(v, columns = [i for i in d_dct])
    data['id'] = [i for i in q_dct]
    return data

def get_symmetric_matrix(df, sym_larger=False):
    '''
    Makes a square distance matrix symmetric
    df = pandas dataframe with columns [id, d_i, ..., d_n]
         id columns should equal d_i, ..., d_n columns
    sym_larger = for symmetric matrices, favor the larger value
    '''
    cols = df.columns[df.columns!='id']
    dist = df[cols].values
    for i in range(0, len(dist)):
        v1 = dist[i,i:]
        v2 = dist[i:,i]
        c1 = (v1 == 0)
        c2 = (v2 == 0)
        c3 = (v2 > v1)
        if sym_larger:
            v1[c3] = v2[c3] # larger value is kept
            dist[i,i:] = v1
            dist[i:,i] = v1
        else: # favor smaller distance, 0 = invalid
            s1 = c1 & ~c2
            s2 = ~c2 & ~c3
            v1[s1] = v2[s1]
            v1[s2] = v2[s2]
            dist[i,i:] = v1
            dist[i:,i] = v1

    # make dist matrix into dataframe
    dist = pd.DataFrame(dist, columns = cols)
    dist['id'] = df['id'].values
    return dist

def run_PCA(df, n_comp=2, recenter=False, ofile='', timestamp=True, log_file='log.txt'):
    '''
    Runs principle component analysis to compression dimensions 
    df = pandas dataframe with columns [id, d_i, ..., d_n]
    n_comps = components to reduce to
    ofile = output file prefix to save pca data
    timestamp = add timestamp to output file
    log_file = log file to write to
    '''
    
    # initialize PCA settings
    append_to_log('running pca with n_comp = '+str(n_comp), log_file)

    # get the vector
    col = df.columns[df.columns!='id']
    v = df[col].values
    
    # center the feature vectors
    if recenter:
        vT = v.T
        for i in range(0,len(vT)):
            vT[i] = vT[i] - np.mean(vT[i])
        v = vT.T
    
    # do PCA
    pca = sklearn.decomposition.SparsePCA(n_components = n_comp)
    v = pca.fit_transform(v)
    
    col = ['f_'+str(i) for i in range(0, n_comp)]
    data = pd.DataFrame(v, columns = col)
    data['id'] = df['id'].values

    if ofile!='':
        ts = ''
        if timestamp:
            ts = '_'+datetime.datetime.fromtimestamp(time.time()).strftime('%Y_%m_%d_%H_%M_%S')
        fname = ofile+ts+'.csv.gz'
        data.to_csv(fname, index = False, compression='infer')
    return data

def run_TSNE(df, n_comp=2, method='sklearn', metric='euclidean', perplexity=30, learning_rate=200, n_iter=1000,
             ofile='', timestamp=True, log_file='log.txt', verbose=0):
    '''
    Run tsne on input data
    df = pandas dataframe with columns [id, v_0, ..., v_i]
    '''
    append_to_log('running sklearn tsne with n_comp = '+str(n_comp), log_file)
    tsne = sklearn.manifold.TSNE(n_components=n_comp, perplexity=perplexity, metric=metric,
                         early_exaggeration=12.0, learning_rate=learning_rate,
                         n_iter=n_iter, n_iter_without_progress=300, verbose=verbose)
    # get the vector
    col = df.columns[df.columns!='id']
    v = tsne.fit_transform(df[col].values)
    cols = ['f_'+str(i) for i in range(0, n_comp)]
    data = pd.DataFrame(v, columns = cols)
    data['id'] = df['id'].values

    if ofile!='':
        ts = ''
        if timestamp:
            ts = '_'+datetime.datetime.fromtimestamp(time.time()).strftime('%Y_%m_%d_%H_%M_%S')
        fname = ofile+ts+'.csv.gz'
        data.to_csv(fname, index=False, compression='infer')
    return data

def cluster_sample(df_q, df_d=[], max_iter=100, block_thresh=0.9, N_rsamp=1, optics_thresh=5000,
                   config ='-k15 -w10 -p 0.9 -D', workspace = './clust_sample/', log_file = 'log.txt'):
    '''
    Find get rough cluster centers for a set of sequence data
    df_q = dataframe of query sequences to search for cluster centers.
    df_d = dataframe of cluster centers with you want to append to
           df_q and df_d must both have at least columns [id, sequence]
    max_iter = max iterations to run
    block_thresh = threshold match_score in order to keep the alignment
                   set higher if you want finer spacing between clusters
    N_rsamp = numbers of random draws for new cluster centers to do
    optics_thresh = threshold number of query sequences between sampling defaults to optics pairwise clustering
    workspace = workspace folder for the aligner
    config = alignment config for minimap2
    log_file = log file to record progress
    '''
    qry = df_q['id'].astype(str).drop_duplicates().values
    if len(df_d) > 0:
        df_centers = df_d.copy()
    else:
        df_centers = []
    
    for j in range(0, max_iter):
        append_to_log('cluster sample, i = '+str(j)+'/'+str(max_iter), log_file)
        append_to_log('query length = '+str(len(qry))+' total clusters = '+str(len(df_centers)), log_file)
        
        # stop if there is nothing more to search
        if len(qry) == 0: 
            append_to_log('query length = 0 terminating loop', log_file)
            break
        
        # run optics if the query size is small enough
        elif len(qry) < optics_thresh:
            seq = df_q[df_q['id'].isin(qry)]
            df_c, df_clst = cluster_iter_optics(seq, center_size = 5, min_samples = 20, config = config,
                                                workspace = workspace, log_file = log_file)
            if len(df_centers) > 0 and len(df_c) > 0:
                df_centers = pd.concat([df_centers[['id','sequence']], df_c])
            elif len(df_centers) == 0 and len(df_c) > 0:
                df_centers = df_c
            break
        
        # pick random sequence to partition the data
        else:
            # draw new centers and record them
            ridx = np.random.permutation(len(qry))
            db = df_q[df_q['id'].isin(qry[ridx[:N_rsamp]])][['id','sequence']]
            if len(df_centers) > 0:
                df_centers = pd.concat([df_centers[['id','sequence']], db])
            else:
                df_centers = db

            seq = df_q[df_q['id'].isin(qry)]
            out = run_minimap2(seq, db, config = config, cigar = True, workspace = workspace)
            out = out[out['match_score'] > block_thresh]
            
            # make the new block to sequence in
            qry = set(qry) - set(out['query_id'].astype(str)) - set(df_centers['id'].astype(str))
            qry = np.array([i for i in qry])
            
    # return the data
    if len(df_centers) > 0:
        df_centers['id'] = ['cluster'+str(i) for i in range(0,len(df_centers))]
    return df_centers

def cluster_refine(df_q, df_d=[], max_iter=1, center_size=20, min_samples=5,
                   N_rsamp=2000, thresh_split=0.9, thresh_merge=0.1, ofile='',
                   timestamp=True, config ='-k15 -w10 -p 0.9 -D',
                   workspace='./clust_refine/', log_file='log.txt'):
    '''
    Refines the cluster centers for a set of sequence data
    df_q = dataframe of query sequences to search for cluster centers.
    df_d = dataframe of cluster centers with you want to append to
           df_q and df_d must both have at least columns [id, sequence]
    center_size = number of sequences from cluster center to use for recomputation of cluster center sequence
    max_iter = max iterations to run
    thresh_split = 
    thresh_merge = 
    N_rsamp = numbers of random draws for new cluster centers to do
    workspace = workspace folder for the aligner
    ofile = output file prefix to use
    config = alignment config for minimap2
    log_file = log file to record progress
    '''
    
    if len(df_d) > 0:
        df_centers = df_d.copy()
    else:
        # use pseudo reference database file as initial cluster centers
        append_to_log('guessing cluster centers via random sampling', log_file)
        df_centers = cluster_sample(df_q, max_iter = 20, block_thresh = 0.8, N_rsamp = 1,
                                    optics_thresh = N_rsamp, config = config,
                                    workspace = workspace, log_file = log_file)
    df_centers['split'] = False

    for k in range(0, max_iter):
        append_to_log('cluster_refine: iter = '+str(k)+'/'+str(max_iter), log_file)
        append_to_log('cluster_refine: sweeping for rare sequences', log_file)
        df_align = run_minimap2(df_q, df_centers, config = config, cigar = True,
                                workspace = workspace, log_file = log_file)
        #df_align = df_align[df_align['tp'] == 'P'] # keep only primary alignments
        df_align = get_best(df_align, ['query_id'], 'match_score')
        
        cols = df_align['database_id'].drop_duplicates().values
        df_centers = df_centers[df_centers['id'].isin(cols)]
        append_to_log('cluster_refine: number of clusters = '+str(len(df_centers)), log_file)
        
        # do the sweep for rare sequences
        p = 5 # get the lower 5 percent of bad alignments
        p = np.percentile(df_align['match_score'].values, p)
        c = df_align['match_score'] > p
        uncover = set(df_q['id'].astype(str)) - set(df_align[c]['query_id'].astype(str))
        append_to_log('cluster_refine: looking at match_score < '+str(p), log_file)
        append_to_log('cluster_refine: query not covered '+str(len(uncover))+'/'+str(len(df_q)), log_file)
        if len(uncover) > 0:
            append_to_log('cluster_refine: looking for rare sequences in '+str(len(uncover))+' reads', log_file)
            uncover = df_q[df_q['id'].isin([i for i in uncover])]
            cout = cluster_sample(uncover, max_iter = 100, block_thresh = 0.8, N_rsamp = 1, 
                                  optics_thresh = N_rsamp, config = config,
                                  workspace = workspace, log_file = log_file)
            if len(cout) > 0:
                cout['split'] = True
                df_centers = pd.concat([df_centers, cout])
                df_centers['id'] = ['cluster'+str(i) for i in range(0,len(df_centers))]
        
        append_to_log('cluster_refine: splitting clusters based on threshold = '+str(thresh_split), log_file)
        df_align = run_minimap2(df_q, df_centers, config = config, cigar = True,
                                workspace = workspace, log_file = log_file)
        
        #df_align = df_align[df_align['tp'] == 'P'] # bin the cluster boundaries
        df_align = get_best(df_align, ['query_id'], 'match_score')
        
        cols = df_align['database_id'].drop_duplicates().astype(str).values
        df_centers = df_centers[df_centers['id'].isin(cols)]
        append_to_log('cluster_refine: number of clusters = '+str(len(df_centers)), log_file)

        d_dct = {cols[i]:i for i in range(0,len(cols))}
        vec = get_feature_vector(df_align[['query_id','database_id','match_score']], d_dct = d_dct)
        x = vec[cols].values.T
        squeue = []
        for i in range(0,len(cols)):
            v = x[i]
            n = v > 0
            c1 = np.percentile(v[n], 25)
            c2 = df_centers[df_centers['id']==cols[i]]['split'].values
            if c1 < thresh_split or c2:
                squeue.append(cols[i])
        
        clst = []
        for i in range(0,len(squeue)):
            qout = squeue[i]
            append_to_log('cluster_refine: doing splitting on '+qout+' '+str(i)+'/'+str(len(squeue)), log_file)
            cout = cluster_recenter([qout], df_q, vec, config = config,
                                    N_rsamp = N_rsamp, thresh_split = 0,
                                    center_size = center_size, min_samples = min_samples,
                                    workspace = workspace, log_file = log_file)
            if len(cout) > 0:
                clst.append(cout)
        
        keep = set(df_centers['id'].astype(str)) - set(squeue)
        df_centers = df_centers[df_centers['id'].isin([i for i in keep])]
        clst.append(df_centers)
        df_centers = pd.concat(clst)
        df_centers['id'] = ['cluster'+str(i) for i in range(0,len(df_centers))]
        
        append_to_log('cluster_refine: merging clusters', log_file)
        df_align = run_minimap2(df_q, df_centers, config = config, cigar = True,
                                workspace = workspace, log_file = log_file)
        cols = df_align['database_id'].drop_duplicates().astype(str).values
        df_centers = df_centers[df_centers['id'].isin(cols)]
        append_to_log('cluster_refine: number of clusters = '+str(len(df_centers)), log_file)
        
        # get the cosimilarity matrix
        d_dct = {cols[i]:i for i in range(0,len(cols))}
        vec = get_feature_vector(df_align[['query_id','database_id','AS']], d_dct = d_dct)
        x = vec[cols].values
        cosim = dist_cosine(x.T, x.T)

        # run hierarchical clustering to find items to merge
        cosim = pd.DataFrame(cosim, columns = cols)
        cosim['id'] = cols
        df_clst = cluster_hierarchical(cosim, metric = 'precomputed', linkage = 'complete',
                                       thresh = 0.8, log_file = log_file)
        # extract clusters to merge
        mqueue = []
        for i in range(0, np.max(df_clst['cluster_id'])+1):
            x = df_clst[df_clst['cluster_id'] == i]['id'].values
            if len(x) > 1:
                mqueue.append([j for j in x])
        n_merge = np.sum(df_clst['cluster_id'] > -1)
        
        append_to_log('cluster_refine: '+str(n_merge)+'/'+str(len(df_centers))+' clusters to merge', log_file)
        
        #c = df_align['tp'] == 'P'
        #vec = get_feature_vector(df_align[c][['query_id','database_id','match_score']], d_dct = d_dct)
        df_align = get_best(df_align, ['query_id'], 'match_score')
        vec = get_feature_vector(df_align[['query_id','database_id','match_score']], d_dct = d_dct)
        
        keep = set(df_centers['id'].astype(str))
        clst = []
        for i in range(0,len(mqueue)):
            mout = mqueue[i]
            append_to_log('cluster_refine: doing merging on '+str(len(mout))+' clusters, '+str(i)+'/'+str(len(mqueue)), log_file)
            cout = cluster_recenter(mout, df_q, vec, config = config,
                                    N_rsamp = N_rsamp, thresh_split = 500,
                                    center_size = center_size, min_samples = min_samples,
                                    workspace = workspace, log_file = log_file)
            if len(cout) > 0:
                clst.append(cout)
            keep = keep - set(mout)
        
        df_centers = df_centers[df_centers['id'].isin([i for i in keep])]
        clst.append(df_centers)
        df_centers = pd.concat(clst)
        df_centers['id'] = ['cluster'+str(i) for i in range(0,len(df_centers))]
        
        # track progress
        if ofile!= '':
            ts = ''
            if timestamp:
                ts = '_'+datetime.datetime.fromtimestamp(time.time()).strftime('%Y_%m_%d_%H_%M_%S')
            fname = ofile+'_'+str(k)+ts+'_c'+str(len(df_centers))+'.csv.gz'
            df_centers.to_csv(fname, index=False, compression='infer')
        
    append_to_log('cluster_refine: complete', log_file)
    return df_centers

def cluster_metric(df_q, df_d, config='-k15 -w10 -p 0.9', 
                   workspace='./clust_metric/', log_file='log.txt'):
    '''
    Outputs some stats on cluster quality and coverage
    '''
    # global alignment
    append_to_log('cluster_metric: getting global alignment on cluster centers', log_file)
    df = run_minimap2(df_q, df_d, config = config, cigar = True, workspace = workspace, log_file = log_file)

    # bin the cluster boundaries
    c = df['tp'] == 'P'
    centers = np.unique(df[c]['database_id'])
    d_dct = {centers[i]:i for i in range(0,len(centers))}
    
    append_to_log('cluster_metric: getting cluster center similarity', log_file)
    df_vec = get_feature_vector(df[['query_id','database_id','AS']], d_dct = d_dct)
    cols = df_vec.columns[df_vec.columns!='id']
    A = df_vec[cols].values
    cosim = dist_cosine(A.T, A.T)
    
    append_to_log('cluster_metric: getting quartiles', log_file)
    vecP = get_feature_vector(df[c][['query_id','database_id','AS']], d_dct = d_dct)
    A = vecP[cols].values
    data = []
    for i in range(0,len(centers)):
        v1 = A.T[i]
        n1 = v1 > 0
        a = np.percentile(v1[n1], 25)
        d = cosim[i,:]
        data.append([centers[i], np.sum(n1), a, np.min(d[d>0])])
    
    out = np.concatenate((data, cosim), axis=1)
    col1 = ['id','n','ci_25','tot_cosim']
    col2 = ['cosim_'+i for i in centers]
    out = pd.DataFrame(out, columns = col1 + col2)
    out = out.merge(df_d[['id','split']], on = 'id', how = 'left')
    
    # formating the data types for output
    for i in col2 + col1[2:]:
        out[i] = out[i].astype(float)
    out['id'] = out['id'].astype(str)
    out['n'] = out['n'].astype(int)
    return out, vecP, df[df['tp']=='P']

def insert_to_queue(items, queue):
    items = {i for i in items}
    for i in range(0,len(queue)):
        if len(items.intersection(queue[i])) > 0:
            queue[i] = queue[i].union(items)
            return queue
    queue.append(items)
    return queue

def dist_cosine(v1, v2):
    '''
    Use numpy to compute cosine similarity. This runs faster than scipy if you have openblas or mkl intalled
    '''
    B = np.matmul(v1, v2.T)
    C = np.reshape([np.sqrt(B[i,i]) for i in range(0,len(B))], (len(B),1))
    C = np.matmul(C, C.T)
    D = B/C
    D[D>1] = 1 # handling numerical imprecision
    D = 1 - D
    return D

def dist_jensenshannon(v1, v2 = []):
    N1 = v1.shape[0]
    p = []
    for i in range(0,N1):
        if i%1000 == 0:
            print('computing p = ',i,'/',N1)
        x = v1[i][v1[i]>0]
        p.append(np.sum(x*np.log(x)))
    
    if len(v2)==0:
        v2 = v1
        q = p
        N2 = N1
        square = True
    else:
        q = []
        N2 = v2.shape[0]
        for i in range(0,N2):
            if i%1000 == 0:
                print('computing q = ',i,'/',N1)
            x = v2[i][v2[i]>0]
            q.append(np.sum(x*np.log(x)))
        square = False
    
    out = np.zeros((N1,N2))
    for i in range(0,N1):
        js = 0
        if square:
            js = i
        if i%50 == 0:
            print('computing m = ',i,'/',N1)
        for j in range(js,N2):
            m = (v1[i]+v2[j])/2
            m = m[m>0]
            out[i,j] = (p[i] + q[j])/2 - np.sum(m*np.log(m))
            out[j,i] = out[i,j]
    return out

def cluster_recenter(queue, df_q, df_vec, config='-k15 -w10 -p 0.9',
                     center_size=20, min_samples=20, N_rsamp=2000, thresh_split=0.9,
                     workspace='./clust_recenter/', log_file='log.txt'):
    '''
    Reclusters some given cluster centers
    
    '''
    # get it via number of sequences closest from the center
    x = []
    for col in queue:
        v = df_vec[['id',col]].values
        v = v[v[:,1]>0]
        if thresh_split == 0:
            s = v
        else:
            s = np.argsort(v[:,1])[::-1]
            v = v[s]
            s = v[:thresh_split]
        x.append(s[:,0])
    qlist = np.unique(np.concatenate((x)))
    
    # check if min_samples is satisfied, otherwise its outlier shit
    if (len(qlist) > min_samples) == False:
        append_to_log('cluster_recenter: not enough samples to cluster', log_file)
        return []

    # subsample if its too big to hold square float array in memory
    append_to_log('clustering '+str(len(qlist))+' sequences', log_file)

    # default sampling
    mem_limit = N_rsamp
    if len(qlist) > mem_limit and len(queue) == 1:
        append_to_log('subsampling because memlimit = '+str(mem_limit)+' sequences reached', log_file)
        ridx = np.random.permutation(len(qlist))
        qlist = qlist[ridx[:N_rsamp]]
    qlist = df_q[df_q['id'].isin(qlist)]
    
    # defaults for cluster recentering
    method = 'pairwise'
    # sample in more clever way if query is large
    if len(qlist) > N_rsamp and len(queue) == 1:
        df_d = cluster_sample(qlist, max_iter = 500, block_thresh = 0.8, N_rsamp = 20,
                              config = config, workspace = workspace, log_file = log_file)
        method = 'pairwise'
    
    # run clustering
    df_c, df_clst = cluster_iter_optics(qlist, center_size = center_size, min_samples = min_samples,
                        config = config, workspace = workspace, log_file = log_file, method = method)
    if len(df_c) > 0:
        df_c['split'] = False
    return df_c

def cluster_iter_optics(qry, center_size=20, min_samples=20,
                        config ='-k15 -w10 -p 0.9 -D', workspace='./clust_iter_optics/',
                        log_file='log.txt', method=''):
    '''
    Find cluster centers for a set of sequence data
    qry = dataframe of query sequences to search for cluster centers.
           Must have at least the following columns:
           [id, sequence]
           [id, sequence, quality]
           [id, sequence, msa_input_file]
           msa_input_files should be in fasta or fastq format
           with filename as .fasta, .fa, .fastq, or .fq 
    df_d = 
    center_size = number of sequences to take from cluster center for multi-seq alignment
    min_samples = min_samples parameter for optics
    workspace = workspace folder for the aligner
    log_file = log file to record progress
    '''
    # make sure names are strings
    df_q = qry.copy()
    df_q['id'] = df_q['id'].astype(str)
    
    # compute pairwise similarity --> slower but more accurate
    append_to_log('cluster_optics_iter: computing pairwise distance matrix', log_file)
    files = get_pairwise_distance(df_q, block=500, output_folder = workspace, workspace = workspace,
                                  config = config, log_file = log_file, symmetric = True)
    df = pd.concat([pd.read_csv(f) for f in files])
    batch_file_remove(files) # clean up the files
    append_to_log('cluster_optics_iter: getting the distance matrix', log_file)
    
    metric = 'match_score'
    df = get_best(df, ['query_id','database_id'], metric)
    df = get_feature_vector(df[['query_id','database_id',metric]], symmetric = True)
    df = get_symmetric_matrix(df, sym_larger = False)
    y = df[df.columns[df.columns!='id']].values
    for i in range(0,len(y)):
        y[i,i] = 1
    y = pd.DataFrame(1-y, columns = ['d'+str(i) for i in range(0,len(y))])
    y['id'] = df['id'].values
    df = y # debug
    
    # safety check that we have enough data to cluster
    if len(df) < min_samples+1:
        append_to_log('cluster_optics_iter: not enough samples to cluster', log_file)
        return [], []
    
    ofile = 'tsne'
    # run clustering
    if method == 'tsne + hdbscan':
        # dimensional reduction with tsne for large datasets
        append_to_log('cluster_optics_iter: running tsne', log_file)
        df = run_TSNE(df, n_comp = 2, metric = 'precomputed', perplexity = min_samples, 
                      ofile = ofile, verbose = 1)
        
        # do importance weighting of the samples
        if weights:
            append_to_log('cluster_optics_iter: applying rca weights to samples', log_file)
            w = df.merge(df_q[['id','N frags']], on = 'id', how = 'left')
            w = np.ceil(np.log(w['N frags'].values)/np.log(3))
            df = cluster_weighting(df, w)
        
        # do clustering
        if method == 'tsne + hdbscan':
            df_clst = cluster_HDBSCAN(df, metric = 'euclidean', min_cluster_size = min_samples,
                                      min_samples = min_samples*2, log_file = log_file)
        df_clst = df_clst[['id','cluster_id']].drop_duplicates()
    else:
        # do optics on raw data
        append_to_log('cluster_optics_iter: running optics', log_file)
        df_clst = cluster_OPTICS(df, metric='precomputed', alt_label = True, log_file=log_file)

    # get cluster centers
    if np.max(df_clst['cluster_id']) == -1:
        append_to_log('No clusters found', log_file)
        return [], df_clst
    else:
        append_to_log('Getting cluster centers', log_file)
        din = []
        df_clst = df_clst.merge(df_q, on='id', how = 'left')
        for i in range(0, np.max(df_clst['cluster_id'])+1):
            c = df_clst['cluster_id'] == i
            din.append(df_clst[c].iloc[:center_size])
        din = pd.concat(din)
        # compute the center sequence
        df_centers = cluster_spoa_merge(din, workspace = workspace, log_file = log_file)
        return df_centers, df_clst

def cluster_weighting(df, weight):
    '''
    Adds more weight to sample datapoint by adding samples
    df = pandas dataframe containing at least the columns ['id']
    weights = integer weights on that data point
    '''
    out = []
    rid = df['id'].values
    for i in range(0, len(weight)):
        out = out + [rid[i]]*int(weight[i])
    out = pd.DataFrame(out, columns = ['id'])
    return out.merge(df, on = 'id', how = 'left')

def cluster_hierarchical(df, metric='precomputed', linkage='single', thresh=0.5, 
                         log_file='log.txt', n_clusters=None):
    '''
    Perform hierarchical clustering on a distance matrix
    '''
    # running the clustering
    append_to_log('Running hierarchical clustering', log_file)

    # initialize settings
    clust = sklearn.cluster.AgglomerativeClustering(n_clusters = n_clusters, distance_threshold = thresh,
                                                    affinity = metric, linkage = linkage)
    cols = df.columns[df.columns!='id'] # get only non read id values
    clust.fit(df[cols].values)
    
    L = clust.labels_
    for i in np.unique(L):
        n = np.sum(L==i)
        if n == 1:
            L[L==i] = -1 # set cluster size = 1 as outliers
    return pd.DataFrame(np.transpose([df['id'].values, L]), columns = ['id','cluster_id'])

def cluster_HDBSCAN(df, metric='euclidean', min_cluster_size=100, min_samples=None,
                    log_file='log.txt', n_jobs=4):
    '''
    Do clustering with HDBSCAN
    df = data frame containing columns [id, d_i..., d_n]
    metric = distance metric to use
    min_samples = min_sample parameter to hdbscan
    min_cluster_size = min_cluster_size parameter to hdbscan
    n_jobs = number of parallel processes to run if supported by algorithm, -1 means use all processors
    log_file = file to log steps of the algorithm
    '''
    # initialize the 
    clust = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples,
                            core_dist_n_jobs=n_jobs, metric=metric)
    
    append_to_log('Running HDBSCAN', log_file)
    cols = df.columns[df.columns!='id'] # get only non read id values
    clust.fit(df[cols].values)
    
    append_to_log('getting and ordering the results', log_file)
    outlier = clust.outlier_scores_
    prob = clust.probabilities_
    labels = clust.labels_
    r_id = df['id'].values
    
    append_to_log('number of clusters = '+str(np.max(labels)+1), log_file)
    append_to_log('unclustered = '+str(np.sum(labels==-1)), log_file)
    
    s = np.lexsort((outlier, labels))
    df = np.transpose([r_id[s], labels[s], range(0,len(r_id)), prob[s], outlier[s]])
    df = pd.DataFrame(df, columns = ['id', 'cluster_id', 'ordering', 'probability', 'outlier_score'])
    df['cluster_id'] = df['cluster_id'].astype(int)
    df['ordering'] = df['ordering'].astype(int)
    df['probability'] = df['probability'].astype(float)
    df['outlier_score'] = df['outlier_score'].astype(float)
    return df

def cluster_OPTICS(df, metric='euclidean', min_samples=None, min_cluster_size=None, xi=0.05,
                   max_eps=None, cluster_method='dbscan', n_jobs=None, alt_label=False,
                   log_file='log.txt'):
    '''
    Function to perform optics via sklearn library implementation
    df = data frame containing columns [id, d_i..., d_n]    
    metric = distance metric to use
    min_samples = number of samples in a neighborhood for a point to be considered a core point
            smaller value results in smaller clusters, but longer run time
            if min_samples = None, the function will search for the best value
    xi = minimum steepness before calling something a cluster boundary
    cluster_method = clustering method to pass to optics, using either dbscan or xi
    max_eps = max distance threshold to check if two points are adjacent to each other
              smaller value results in faster run time at the cost of accuracy
    alt_label = use custom reachability plot labeling scheme that is not in sklearn
    n_jobs = number of parallel processes to run
    log_file = file to log steps of the algorithm
    '''
    # running the clustering
    append_to_log('Running OPTICS', log_file)
    cols = df.columns[df.columns!='id'] # get only non read id values
    if max_eps == None:
        z = np.abs(df[cols].values)
        max_eps = (np.max(z)+np.min(z))/2
    append_to_log('max_eps = '+str(max_eps), log_file)
    
    # do optimization on min_samples to get clustering with the least outliers
    if min_samples == None or min_samples > len(df):
        min_samples = np.ceil(len(df)/2)
        max_iter = 100
    else: # if min_samples is provided, dont optimize
        max_iter = 1
    
    # initialize settings
    prev_nout = len(df)
    prev_nclust = 0
    prev_min = min_samples*2
    best_clust = None
    best_nclust = 0
    best_min = prev_min
    best_nout = prev_nout
    delta = min_samples/2
    for i in range(0, max_iter):
        min_samples = int(min_samples)
        append_to_log('clust_OPTICS: iter='+str(i)+' using min_samples='+str(min_samples), log_file)
        clust = sklearn.cluster.OPTICS(min_samples=min_samples, min_cluster_size=min_cluster_size, xi=xi,
                            max_eps=max_eps, cluster_method=cluster_method, metric=metric, n_jobs=n_jobs)
        clust.fit(df[cols].values)
        # check labels for cluster number and coverage
        nout = np.sum(clust.labels_ == -1)
        nclust = np.max(clust.labels_)+1
        append_to_log('clust_OPTICS: clusters='+str(nclust)+' outliers='+str(nout)+' delta='+str(delta), log_file)
        # always keep the change if we get more clusters
        if nclust > prev_nclust or (nclust == prev_nclust and nout <= prev_nout):
            tmp = min_samples
            delta = int((prev_min - min_samples)/2)
            min_samples-= delta
        # if cluster number goes down as we increase min_samples
        if nclust < prev_nclust:
            tmp = min_samples
            delta = int((prev_min + min_samples)/2) - min_samples
            min_samples = prev_min + delta
        # record the best clustering
        if nclust > best_nclust or (nclust == best_nclust and nout <= best_nout):
            if np.min(clust.reachability_) >= 0:
                best_clust = clust
                best_min = min_samples
                best_nclust = nclust
                best_nout = nout
        # if shift in min_samples is zero, exit loop
        if delta==0 or min_samples <= 2:
            break
        prev_nclust = nclust
        prev_nout = nout
        prev_min = tmp

    # get the rest of the values
    s = best_clust.ordering_
    reach = best_clust.reachability_[s]
    labels = best_clust.labels_[s]
    r_id = df['id'].values[s]
    
    # cap the reachability so fft can be computed
    y = np.copy(reach)
    y[reach >= np.inf] = -1
    reach[y == -1] = np.max(y)
    
    if np.min(reach) >= 0:
        # using custom relabeling on the reachability plot
        if alt_label:
            append_to_log('using alt labeling', log_file)
            labels = reachability_alt_label(reach)
        append_to_log('n_clusters='+str(np.max(labels)+1)+' n_unclustered='+str(np.sum(labels==-1))+' N='+str(len(df)), log_file)
    else:
        append_to_log('clust_OPTICS: reachability < 0, please change the min_samples parameter you are using')
    
    # format and return the data
    df = np.transpose([r_id, labels, range(0,len(r_id)), reach])
    df = pd.DataFrame(df, columns = ['id', 'cluster_id', 'ordering', 'reachability'])
    df['cluster_id'] = df['cluster_id'].astype(int)
    df['ordering'] = df['ordering'].astype(int)
    df['reachability'] = df['reachability'].astype(float)
    return df

def reachability_alt_label(reach):
    # get the slope change
    d1 = sp.signal.fftconvolve(reach, [1,-1], mode='same')
    d1[:-1] = d1[1:]
    '''
    d2 = sp.signal.fftconvolve(d1,[1,-1], mode='same')
    # compute local median and stdev
    L = 100
    f = np.ones(L)
    med = sp.ndimage.median_filter(d1, size=L)
    # use mean shift to reduce numerical imprecision for variance calc
    K = np.median(d1)
    m = sp.signal.fftconvolve(d1-K, f, mode='same')
    v = sp.signal.fftconvolve((d1-K)**2, f, mode='same')
    s = v/(L-1) - m**2/L/(L-1) # get the local standard deviation
    '''
    labels = np.zeros(len(reach))
    # flag steep drops as outliers
    km = sklearn.cluster.AgglomerativeClustering(n_clusters=5, affinity='euclidean', linkage='ward')
    km.fit(np.transpose([d1]))
    c = []
    for j in np.unique(km.labels_):
        x = d1[km.labels_==j]
        c.append(np.mean(x))
    s = np.argsort(c)
    labels[(km.labels_==s[0])] = -1

    # flag high reachability as outliers
    y = np.copy(reach)
    y[(km.labels_==s[0])] = np.max(y)
    km = sklearn.cluster.AgglomerativeClustering(n_clusters=3, affinity='euclidean', linkage='ward')
    km.fit(np.transpose([y]))
    c = []
    for j in np.unique(km.labels_):
        x = y[km.labels_==j]
        c.append(np.mean(x))
    s = np.argsort(c)
    labels[(km.labels_==s[-1])] = -1

    # assign new labels
    th = labels!=-1
    c = th[1:]!=th[:-1]
    cuts = np.arange(1,len(c)+1)[c]
    k = 0
    for i in range(1,len(cuts),2):
        s1 = cuts[i-1]
        s2 = cuts[i]
        labels[s1:s2] = k
        k+=1
    return labels

def cluster_Kmeans(df, n_clusters, n_init=10, n_iter=100, log_file='log.txt'):
    '''
    Function to perform clustering via k means algorithm
    '''
    # running the clustering
    append_to_log('Running kmeans with n_clusters = '+str(n_clusters), log_file)
    cols = df.columns[df.columns!='id'] # get only non read id values
    x = df[cols.values]
    
    # initialize settings
    clust = sklearn.cluster.MiniBatchKMeans(n_clusters = n_clusters, n_init = n_init, max_iter = n_iter)
    clust.fit(x)
    
    append_to_log('Getting results', log_file)
    labels = clust.labels_
    inertia = clust.inertia_
    centers = clust.cluster_centers_
    r_id = df['id'].values
    
    # build center distance matrix
    c = np.zeros(x.shape)
    for i in range(0,len(x)):
        c[i] = centers[labels[i]]
    d = np.std(x-c, axis = 1) # get squared distance
    s = np.lexsort((d,labels))
    
    df = np.transpose([r_id[s], labels[s], range(0,len(r_id)), d[s]])
    df = pd.DataFrame(df, columns = ['id', 'cluster_id', 'ordering', 'inertia'])
    df['cluster_id'] = df['cluster_id'].astype(int)
    df['ordering'] = df['ordering'].astype(int)
    return df

def cluster_spoa_merge(df, config = '-n -15 -g -10 -l 0 -r 2', workspace = './clust_spoa_merge/',
                       batch_size = 100, cleanup = True, log_file = 'log.txt'):
    '''
    Function to return list of cluster centers via multi-sequence alignment
    df = dataframe of query sequences to search for cluster centers.
        Must have at least the following columns:
        [id, cluster_id, sequence]
        [id, cluster_id, sequence, quality]
        [id, cluster_id, msa_input_file]
        msa_input_files should be in fasta or fastq format with filename:
        .fasta, .fa, .fastq, or .fq 
    config = config passed to multi-sequence aligner
    workspace = folder containing data from multi-sequence alignments
    batch_size = batch size of files to align
    cleanup = delete workspace folder after the run
    '''
    # make directory if it does not exist
    check_dir(workspace)
    files = glob.glob(workspace+'*')
    batch_file_remove(files)
    
    # sort the dataframe and prepare to iterate through it
    df = df.sort_values(by = ['cluster_id', 'id'])
    # check if sequences are consensus and if original msa files are present
    MSA_infile = []
    # deal with merging consensus from multi-seq alignment
    if 'msa_input_file' in df.columns:
        df = df[['id','cluster_id','msa_input_file']]
        df_fa = []
        df_fq = []
        for i in range(0, len(df)):
            [r_id, c_id, file] = df.values[i]
            ftype = file.split('.')
            if 'fa' in ftype or 'fasta' in ftype:
                df_fa.append(read_fasta(file))
            elif 'fq' in ftype or 'fastq' in ftype:
                df_fq.append(read_fastq(file))

            # write file if boundary of new cluster found
            if i+1 >= len(df) or c_id!= df.values[i+1, 1]:
                c1 = len(df_fa) > 0
                c2 = len(df_fq) > 0
                if c1:
                    df_fa = pd.concat(df_fa)
                if c2:
                    df_fq = pd.concat(df_fq)
                if c1 and c2:
                    fname = workspace + str(c_id) + '.fa'
                    seqlist = pd.concat([df_fa[['id','sequence']], df_fq[['id','sequence']]])
                    write_fasta(fname, seqlist[['id','sequence']].values)
                elif c1:
                    fname = workspace + str(c_id) + '.fa'
                    write_fasta(fname, df_fa[['id','sequence']].values)
                elif c2:
                    fname = workspace + str(c_id) + '.fq'
                    write_fastq(fname, df_fq[['id','sequence','quality']].values)
                df_fa = []
                df_fq = []
                MSA_infile.append(fname)
    
    # deal with sequences that have quality info
    elif 'quality' in df.columns:
        df = df[['id','cluster_id','sequence','quality']]
        seqlist = []
        for i in range(0, len(df)):
            [r_id, c_id, seq, qual] = df.values[i]
            seqlist.append([r_id, seq, qual])
            # write file if boundary of new cluster found
            if i+1 >= len(df) or c_id!= df.values[i+1, 1]:
                fname = workspace + str(c_id) + '.fq'
                write_fastq(fname, np.array(seqlist))
                seqlist = []
                MSA_infile.append(fname)
    
    # deal with sequence only
    elif 'sequence' in df.columns:
        df = df[['id','cluster_id','sequence']]
        seqlist = []
        for i in range(0, len(df)):
            [r_id, c_id, seq] = df.values[i]
            seqlist.append([r_id, seq])
            # write file if boundary of new cluster found
            if i+1 >= len(df) or c_id!= df.values[i+1, 1]:
                fname = workspace + str(c_id) + '.fa'
                write_fasta(fname, np.array(seqlist))
                seqlist = []
                MSA_infile.append(fname)
    
    # run multi-sequence alignment when batch_size limit is hit
    N_chunks = np.ceil(len(MSA_infile)/batch_size)
    for i, infile in enumerate(np.array_split(MSA_infile, N_chunks)):
        run_msa(infile, 'spoa', config)
        append_to_log('cluster_spoa_merge: spoa on '+str(i)+'/'+str(N_chunks), log_file)
    
    # get consensus
    append_to_log('cluster_spoa_merge: reading consensus', log_file)
    MSA_outfile = [i+'.out' for i in MSA_infile]
    data = []
    for i in range(0,len(MSA_outfile)):
        out = read_spoa(MSA_outfile[i])
        c_id = 'cluster' + MSA_outfile[i].split('/')[-1].split('.')[0]
        data.append([c_id, out[0]])

    # remove work folder after completion
    if cleanup:
        subprocess.run(['rm', '-r', workspace])
    return pd.DataFrame(data, columns=['id', 'sequence'])

def cigar_align(query, ref, cigar):
    '''
    Function to align query and reference strings based on CIGAR string
    query = query sequence --> this must be relative to how it was used in aligner
    ref = reference sequence
    cigar = CIGAR string from aligner
    
    return [aligned query, aligned reference, alignment ticks]
    '''
    pattern = re.compile('([0-9]*)([DIMSX])')
    query = query.upper()
    ref = ref.upper()
    q_aligned = ''
    ref_aligned = ''
    align = ''
    i1 = 0
    i2 = 0
    for n, c in pattern.findall(cigar):
        span = int(n)
        if c == 'M':
            q_aligned+=query[i1:i1+span]
            ref_aligned+=ref[i2:i2+span]
            align+=''.join(['|']*span)
            i1+=span
            i2+=span
        elif c=='X':
            q_aligned+=query[i1:i1+span]
            ref_aligned+=ref[i2:i2+span]
            align+=''.join([' ']*span)
            i1+=span
            i2+=span
        elif c=='S' or c=='H':
            q_aligned+=query[i1:i1+span].lower() # query contains soft clipped sequence
            ref_aligned+=ref[i2:i2+span]
            align+=''.join([' ']*span)
            i1+=span
            i2+=span
        elif c=='I':
            q_aligned+=query[i1:i1+span] # query contains an insert
            ref_aligned+=''.join(['-']*span) # adds gap to reference
            align+=''.join([' ']*span)
            i1+=span
        elif c=='D':
            q_aligned+=''.join(['-']*span) # query contains a deletion --> add a gap
            ref_aligned+=ref[i2:i2+span]
            align+=''.join([' ']*span)
            i2+=span
    return [q_aligned, ref_aligned, align]

def print_pairwise_alignment(data, chars=100):
    '''
    Formatted printing of a set of aligned strings
    data = [query, reference, alignment]
    chars = chars per line
    '''
    [query, ref, aligned] = data
    out = ''
    x = np.arange(0,len(query),chars)
    if x[-1] < len(query):
        x = np.concatenate((x,[len(query)]))
    for i in range(0,len(x)-1):
        out+=query[x[i]:x[i+1]]+  ' qry \n'
        out+=aligned[x[i]:x[i+1]] + ' '+str(x[i+1])+'\n'
        out+=ref[x[i]:x[i+1]]+ ' ref \n\n'
    print(out)

def print_multi_alignment(data, chars=100):
    '''
    Formatted printing of a set of aligned strings
    data = pandas array of sequences with the columns ['id','sequence']
    chars = chars per line
    '''
    df = data[['id','sequence']].values
    out = ''
    x = np.arange(0,len(df[0,1]),chars)
    if x[-1] < len(df[0,1]):
        x = np.concatenate((x,[len(df[0,1])]))
    for i in range(0,len(x)-1):
        for j in range(0,len(df)):
            out+=df[j,1][x[i]:x[i+1]]+' '+df[j,0]+' '+str(x[i+1])+'\n'
        out+='\n'   
    print(out)
    
def dna_revcomp(seq):
    '''
    Return the reverse complement of a DNA string
    '''
    key = {'A':'T','T':'A','G':'C','C':'G',
           'R':'Y','Y':'R',
           'S':'S','W':'W',
           'K':'M','M':'K',
           'B':'V','V':'B','D':'H','H':'D',
           '-':'-','N':'N',' ':' '}
    seq = seq.upper()
    return ''.join([key[seq[-i]] for i in range(1,len(seq)+1)])

def get_nondegen_sequence(seq):
    '''
    Convert a sequence with degenerate letters to all possible combination of degenerate letters.
    seq = <string>
    return <array> of strings
    '''
    base = ['AG','CT','GC',
            'AT','GT','AC',
            'CGT','AGT','ACT','ACG','ATGC']
    base_dict = {'R':0,'Y':1,'S':2,'W':3,'K':4,'M':5,'B':6,'D':7,'H':8,'V':9,'N':10}
    
    data = ['']
    for i in range(0,len(seq)):
        if seq[i] in base_dict:
            data = add_letters(data, base[base_dict[seq[i]]])
        else:
            data = add_letters(data, seq[i])
    return data

def add_letters(data, lets):
    out = []
    for i in range(0,len(data)):
        for j in range(0,len(lets)):
            out.append(data[i]+lets[j])
    return out

def stats_subsample(df, col, N=100):
    '''
    Subsamples the group
    df = pandas dataframe with atleast columns [group_id]
    col = column to subsample by
    N = number of samples for each group
    
    returns dataframe which is a subsample of original group
    '''
    groups = np.unique(df[col])
    out = []
    for g in groups:
        x = df[df[col] == g]
        ridx = np.random.permutation(len(x))
        out.append(x.iloc[ridx[:N]])
    return pd.concat(out)

def stats_count_xygroup(df, xgroup, ygroup):
    '''
    df = pandas dataframe with columns = [id, xgroup, ygroup]
    xgroup could be primers or method
    ygroup can be haplotypes
    '''
    # get counts
    df = df.values
    data = np.zeros((len(xgroup),len(ygroup)))
    
    for i in range(0,len(xgroup)):
        for j in range(0,len(ygroup)):
            data[i,j] = np.sum((df[:,1] == xgroup[i])&(df[:,2] == ygroup[j]))
    return data

def stats_group_data(df, groups = None):
    '''
    Reformats dataframe into array for plotting jitter or violin plots
    df = pandas dataframe with columns = [group_id, val]
    groups = array of group_id to plot. If None is given, the list is derived from given data
    
    returns array with each row containing [[group_id],[val_0 ... val_n]]
    '''
    df = df.values

    # get group labels
    if type(groups) == type(None):
        groups = np.unique(df[:,0])
    
    # populates the values for each group
    data = []
    for g in groups:
        x = df[df[:,0] == g]
        if len(x) > 0:
            data.append([g, x[:,1]])
        else:
            data.append([g, [-1]])
    return data

def mpl_violin(ax, data, groups=None, vert=False, alpha=0.5,
               face_color='gray', edge_color='black'):
    '''
    Use matplotlib to make violin plots
    axis = axis handle obtained via plt.gca()
    data = dataframe containing columns with [group_id, values]
    groups = array of group_id to plot
    vert = make vertical box whiskers or horizonal
    
    face_color = color of each face of the violin
    edge_color = color of the edge of the violin
    '''
    # get the data
    d = stats_group_data(data, groups)

    # make the violin plot
    labels = []
    y = []
    for i in range(0,len(d)):
        if len(d[i][1]) > 1:
            x = np.array(d[i][1]).astype(float)
            parts = ax.violinplot([x], positions = [i], vert=vert,
                                  showmeans=False, showmedians=False, showextrema=False)
            # coloring for the violin
            for pc in parts['bodies']:
                pc.set_facecolor(face_color)
                pc.set_edgecolor(edge_color)
                pc.set_alpha(alpha)

def mpl_box_whisker(ax, data, groups=None, vert=False, alpha=1, median_color='white',
                    cmap=['gray'], text_x_offset=0, text_y=None, 
                    counts=False, labels=False):
    '''
    Use matplotlib to make box whisker plot
    axis = axis handle obtained via plt.gca()
    data = dataframe containing columns with [group_id, values]
    groups = array of group_id to plot
    vert = make vertical box whiskers or horizonal
    median_color = color for median data point
    alpha = alpha of each box whisker
    cmap = array of colors for each box whiskers

    count = plot the counts for each group_id
    text_x_offset = how much the count annotate along each x index
    text_y = how much to offset the count index along each y index
    '''
    # get the data
    data = stats_group_data(data, groups)
    
    # process the percentiles
    inds = []
    x = []
    L = []
    for i in range(0,len(data)):
        inds.append(i)
        L.append(data[i][0])
        a = np.percentile(data[i][1], [10, 25, 50, 75, 90])
        n = len(data[i][1])
        x.append([a[0], a[1], a[2], a[3], a[4], n])
    x = np.array(x)
    L = np.array(L)
    
    # positioning for counts information
    if type(text_y) == type(None):
        text_y = np.max(x[:,4])
    
    # make the box whiskers
    if vert:
        for i in range(0,len(x)):
            ax.vlines(inds[i], x[i,0], x[i,4], color = cmap[i%len(cmap)], linestyle='-', lw=1, alpha = alpha)
            ax.vlines(inds[i], x[i,1], x[i,3], color = cmap[i%len(cmap)], linestyle='-', lw=5, alpha = alpha)
            ax.scatter(inds[i], x[i,2], marker='o', color=median_color, s=30, alpha = alpha, zorder = 3)
            if counts:
                ax.text(inds[i]+text_x_offset, text_y, str(int(x[i,5])), rotation=90)
    else:
        for i in range(0,len(x)):
            ax.hlines(inds[i], x[i,0], x[i,4], color = cmap[i%len(cmap)], linestyle='-', lw=1, alpha = alpha)
            ax.hlines(inds[i], x[i,1], x[i,3], color = cmap[i%len(cmap)], linestyle='-', lw=5, alpha = alpha)
            ax.scatter(x[i,2], inds[i], marker='o', color=median_color, s=30, alpha = alpha, zorder = 3)
            if counts:
                ax.text(text_y, inds[i]+text_x_offset, str(int(x[i,5])))
                
    if vert and labels:
        mpl_set_xaxis(ax, inds, L)
    elif ~vert and labels:
        mpl_set_yaxis(ax, inds, L)
        
def mpl_set_xaxis(ax, inds, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.set_xticks(inds)
    ax.set_xticklabels(labels)
    ax.set_xlim(inds[0]-1,inds[-1]+1)
    # rotate the tick labels
    for tick in ax.get_xticklabels():
        tick.set_rotation(90)

def mpl_set_yaxis(ax, inds, labels):
    ax.get_yaxis().set_tick_params(direction='out')
    ax.set_yticks(inds)
    ax.set_yticklabels(labels)
    ax.set_ylim(inds[0]-1,inds[-1]+1)
