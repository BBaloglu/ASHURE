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

def init_log(fname=None, level='DEBUG'):
    fmt = 'pid['+str(os.getpid())+'] %(asctime)s.%(msecs)03d %(levelname)s: %(message)s'
    dfmt = '%Y-%m-%d %H:%M:%S'
    formatter = logging.Formatter(fmt, datefmt=dfmt)
    # print log to a file
    N=1
    if type(fname)==str:
        file_handler = logging.FileHandler(filename=fname)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        N=2 # adds file_handler only once
        if len(logging.getLogger().handlers) < N:
            logging.getLogger().addHandler(file_handler)
    # print to stdout
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(level)
    stdout_handler.setFormatter(formatter)
    # add stream handler only once
    if len(logging.getLogger().handlers) < N:
        logging.getLogger().addHandler(stdout_handler)
    logging.getLogger().setLevel(level)

def parse_ONT_header(text):
    text = text.split(' ')
    key = {'read':1,'ch':2,'start_time':3}
    out = [None]*(len(key)+1)
    out[0] = text[0]
    for line in text[1:]:
        line = line.split('=')
        if line[0] in key:
            out[key[line[0]]] = line[1] 
    return out

def add_ONT_header(df):
    df = df.copy()
    data = []
    for h in df['id'].values:
        data.append(parse_ONT_header(h))
    data = pd.DataFrame(data, columns=['id','read','ch','start_time'])
    for c in data.columns:
        df[c] = data[c].values
    col = ['read','ch']
    for c in col:
        df[c] = df[c].astype(int)
    return df

# Function to load fastq data and primer info
def load_ONT_fastq(files, low_mem=False):
    '''
    Function to load fastq info from list of ONT fastq files
    files = list of fastq files
    output = pandas dataframe of ONT fastq info
    '''
    # Generate table of sequences via iteration
    start = time.time()
    df = []
    for i in range(0,len(files)):
        df.append(read_fastq(files[i], low_mem=low_mem)) # Save tables as a vector
        if i%20==0: # print out progress
            logging.info('Processing '+str(i)+'/'+str(len(files))+', current file = '+files[i])
            logging.info('elapse = {0:.2f}'.format(time.time()-start)+'s')
    df = pd.concat(df) # Merge all the tables
    logging.info('adding ONT header info')
    df = add_ONT_header(df)
    # add sequence length info
    if 'rlen' in df.columns:
        df['length'] = df['rlen']
    elif 'sequence' in df.columns:
        df['length'] = [len(i) for i in df['sequence']]
    # format datetime
    df['start_time'] = pd.to_datetime(df['start_time'])
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
    logging.info('Reading primers from '+primer_file)
    df = pd.read_csv(primer_file, delimiter = ',')
    # safety checks on the primer data
    L1 = len(df)
    df = df.dropna()
    L2 = len(df)
    if L1!=L2:
        logging.warning('Warning: nan info was dropped from primer data')
    
    # check basic info is there
    col = ['fwd_id','fwd_seq','rev_id','rev_seq']
    for c in col:
        if (c in df.columns) == False:
            logging.error('column '+c+' not found in '+primer_file)
            sys.exit(1)
    return df

def add_seq(df):
    '''
    Adds sequence and quality information back to dataframe
    df = dataframe containing at least the columns:
        [filename, id, seek1, rlen] for fasta files
        [filename, id, seek1, seek2, rlen] for fastq files
    return dataframe with columns [id, sequence, quality] or [id, sequence]
    '''
    # if data is already there, dont do anything
    if 'sequence' in df.keys():
        return df
    # check if relevant columns are there
    cols = ['filename','id','seek1','rlen']
    for c in cols:
        if (c in df.keys())==False:
            logging.warning('columns '+c+' not found in dataframe')
            return df
    # sort and iterate
    df = df.sort_values(by=['filename','seek1']).copy()
    data = []
    if 'seek2' in df.keys():
        # work on fastq
        x = df[['filename','id','seek1','seek2','rlen']].values
        files = np.unique(x[:,0])
        for fname in files:
            y = x[x[:,0]==fname]
            with open(fname,'r') as f:
                for i in range(0,len(y)):
                    f.seek(y[i,2])
                    seq = f.read(y[i,4])
                    f.seek(y[i,3])
                    q = f.read(y[i,4])
                    data.append([seq, q])
        data = np.array(data)
        df['sequence'] = data[:,0]
        df['quality'] = data[:,1]
        cols = ['seek1','seek2','rlen','filename']
        df = df.drop(columns=cols)
    else:
        # work on fasta
        x = df[['filename','id','seek1','rlen']].values
        files = np.unique(x[:,0])
        for fname in files:
            y = x[x[:,0]==fname]
            with open(fname,'r') as f:
                for i in range(0,len(y)):
                    f.seek(y[i,2])
                    seq = f.read(y[i,3])
                    seq = seq.replace('\n','')
                    data.append(seq)
        data = np.array(data)
        df['sequence'] = data
        cols = ['seek1','rlen','filename']
        df = df.drop(columns=cols)
    return df

def read_fastq(fname, low_mem=False):
    '''
    Reads list of sequences from a fastq file format.
    
    fname = fname to read the information
    low_mem = saves memory by not storing sequence and quality strings.
            Instead this function returns columns [filename, id, seek1, seek2, rlen]
            Use add_seq function to add back sequence and quality to the dataframe
    returns list of names, sequence, quality
    '''
    data = []
    with open(fname,'r') as f:
        fsize = os.path.getsize(fname)
        while fsize > f.tell():
            rid = f.readline()[:-1]
            s1 = f.tell()
            seq = f.readline()[:-1]
            f.seek(f.tell()+2)
            q = f.readline()[:-1]
            L = len(seq)
            if low_mem:
                data.append([rid[1:], s1, L])
            else:
                data.append([seq, q, rid[1:]])
    col = ['sequence','quality','id','seek1','rlen']
    if low_mem:
        data = pd.DataFrame(data, columns=col[2:])
        data['seek2'] = data['seek1']+data['rlen']+3
        data['filename'] = fname
    else:
        data = pd.DataFrame(data, columns=col[:3])
    return data

def write_fastq(fname, data):
    '''
    Write list of sequences into a fastq file format.
    
    data =  array with [id, sequence, quality] as the columns
    fname = file to write the information
    '''
    with open(fname,'w') as f:
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

def read_fasta(fname, low_mem=False):
    '''
    Reads list of sequences from a fasta file format.
    
    fname = file to read the information
    returns list of sequence and names
    '''
    data = []
    with open(fname,'r') as f:
        seq = ''
        fsize = os.path.getsize(fname)
        while fsize > f.tell():
            line = f.readline()
            if line[0]=='>':
                # record it once new line is hit
                if seq!='':
                    L = len(seq)
                    seq = seq.replace('\n','')
                    if low_mem:
                        data.append([rid, s1, L])
                    else:
                        data.append([seq, rid])
                # start the new entry 
                rid = line[1:-1]
                s1 = f.tell()
                seq = ''
            else:
                seq+= line
    # end of file reached, run last loop
    L = len(seq)
    seq = seq.replace('\n','')
    # export the data
    col = ['sequence','id','seek1','rlen']
    if low_mem:
        data.append([rid, s1, L])
        data = pd.DataFrame(data, columns=col[1:])
        data['filename'] = fname
    else:
        data.append([seq, rid])
        data = pd.DataFrame(data, columns=col[:2])
    return data

def write_fasta(fname, data):
    '''
    Write list of sequences into a fasta file format.
    
    data =  array with [id, sequence] as the columns
    fname = file to write the information
    '''
    with open(fname,'w') as f:
        if len(data.shape) > 1: # case for table
            for i in data:
                text = '>' + str(i[0]) + '\n'
                text+= str(i[1]) + '\n'
                f.write(text)
        else: # case for single row
            text = '>' + str(data[0]) + '\n'
            text+= str(data[1]) + '\n'
            f.write(text)

def get_SAM_info(read, key):
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
    info = np.array([0]*len(key))
    if rname != '*':
        for i in range(11,len(read)):
            x = read[i].split(':') # parse it backwards
            if x[0] in key.keys():
                info[key[x[0]]] = x[-1]
    data = [qname, rname, orientation, flag, pos1, pos2, mapq, cigar]
    return data, info

def parse_SAM(fname):
    '''
    This function parses the information from a SAM file and returns a dataframe   
    '''
    fsize = os.path.getsize(fname)
    record = False
    col1 = ['query_id','database_id','orientation','flag','t_start','t_end','mapq','CIGAR']
    col2 = ['AS','XS','XN','XM','XO','XG','NM']
    key = {col2[i]:i for i in range(0,len(col2))}
    data = []
    with open(fname,'r') as f:
        while fsize > f.tell():
            text = f.readline()
            # dont record anything until we read past the header:w
            if ~record and text[:3]=='@PG':
                record = True
                text = f.readline()
            # start recording data
            if record:
                out, info = get_SAM_info(text, key)
                if out[3][-2] == '1':
                    text = f.readline()
                    out2, info2 = get_SAM_info(text, key)
                    out[0]+= ' '+out2[0]
                    out[7]+= ' '+out2[7]
                    info = info+info2
                data.append(out+[j for j in info])
    return pd.DataFrame(data, columns=col1+col2)

def run_bowtie2(query, database, workspace='./bowtie2/', config='-a --very-sensitive-local --threads 1 --quiet', build_index=True, cleanup=False):
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
    workspace = check_dir(workspace)
    # load sequence incase low mem option is used 
    query = add_seq(query) 
    database = add_seq(database)
    # Define locations of input and output files for bowtie
    qfile = workspace+'read1.fa'
    mfile = workspace+'read2.fa'
    dbfile = workspace+'database.fa'
    outfile = workspace+'results.sam'
    btfile = workspace+'index'

    # Only build index when we need to. Index building takes a lot of time
    if build_index == True:
        value = check_seqdf(database)
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
    value = check_seqdf(database)
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
    data = parse_SAM(outfile)

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

def run_bwa(query, database, workspace='./bwa/', config=' mem -a ', build_index=True, cleanup=False):
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
    workspace = check_dir(workspace)
    # load sequence incase low mem option is used 
    query = add_seq(query) 
    database = add_seq(database)
    # Define locations of input and output files for bowtie
    qfile = workspace+'read1.fq'
    mfile = workspace+'read2.fq'
    dbfile = workspace+'database.fa'
    outfile = workspace+'results.sam'

    # Only build index when we need to. Index building takes a lot of time
    if build_index:
        value = check_seqdf(database)
        if value == 2:
            write_fasta(dbfile, database[['id','sequence','quality']].values)
        elif value == 1:
            write_fasta(dbfile, database[['id','sequence']].values)
        else:
            return -1
        cmd = 'bwa index '+dbfile
        subprocess.run(cmd.split(' '))

    # Check if we are aligning paired reads
    value = check_seqdf(query)
    
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
    data = parse_SAM(outfile)

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

def parse_PAF(fname):
    '''
    This function parses the information from a PAF file and returns a dataframe
    '''
    with open(fname,'r') as f:
        text = f.read().split('\n')
    
    if len(text[-1])==0:
        text = text[:-1] # removes the last line that is empty    
    data = []
    
    # information we are parsing
    col1 = ['query_id','q_len','q_start','q_end','orientation','database_id','t_len','t_start','t_end','match','tot','mapq']
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
    data = data.rename(columns={'cg':'CIGAR'})
    return data

def run_minimap2(query, database, workspace='./minimap2/', config='-x map-ont', cigar=True, build_index=True, use_index=False, cleanup=False):
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
    workspace = check_dir(workspace)
    # load sequence incase low mem option is used 
    query = add_seq(query) 
    database = add_seq(database)

    # Define locations of input and output files for bowtie
    qfile = workspace+'read1.fq'
    mfile = workspace+'read2.fq'
    dbfile = workspace+'database.fq'
    outfile = workspace+'results.paf'
    btfile = workspace+'index.mmi'
    
    # Write database to fasta or fastq
    value = check_seqdf(database)
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
    value = check_seqdf(query)
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
    data = parse_PAF(outfile)
    
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
    if cleanup and ~use_index:
        subprocess.run(['rm','-r',workspace])
    return data

def check_seqdf(df):
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
        logging.error('columns of data = '+''.join(cols))
        logging.error('''
        Input dataframe does not have the proper header.
        
        The following are valid columns
        
        fasta data   = [id, sequence]
        fastq data   = [id, sequence, quality]
        paired fasta = [fwd_id, fwd_seq, rev_id, rev_seq]
        paired fastq = [fwd_id, fwd_seq, fwd_quality, rev_id, rev_seq, rev_quality]
        ''')
        return -1

def check_dir(folder, overwrite=True, fname_only=False):
    '''
    This function checks if the directory exists. If it does not exist, make the directory
    '''
    # make folder names consistent
    folder = folder.split('/')
    tmp = ''
    for j in folder:
        if j!='':
            tmp+=j+'/'
    folder = tmp
    # don't make directory if not asked to
    if fname_only:
        return folder
    # make directory if it does not exist
    if os.path.exists(folder) == False and overwrite:
        logging.info('Making directory '+folder)
        cmd = 'mkdir ' + folder
        subprocess.run(cmd.split(' '))
        return folder
    elif overwrite:
        return folder
    else:
        # make find a new name if it does exist
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
    return data.merge(clusterN, on='cluster_num', how='left').set_index(df.index[s]) # merge and reset index

def split_overlaps(df, thresh=1):
    '''
    Function that iterates through list and reduces it to best scoring non-overlapping id
    df = dataframe with columns [id, start, stop, score]
    thresh = allow some threshold nucleotides of overlap before calling it an overlap
    Returns cleaned items
    '''
    data = df.values # work on numpy array because its faster
    s = np.lexsort((data[:,3],data[:,0]))[::-1]
    data = data[s]
    # initial value for keep --> keep everything that is big enough
    keep = (thresh*2 < data[:,2] - data[:,1])
    for i in range(0,len(data)):
        if keep[i]:
            for j in range(i+1,len(data)):
                # exit loop if end of read is reached
                if data[i,0]!=data[j,0]:
                    break
                elif keep[j]:
                    # q_stop < cur_start or q_start > cur_stop --> then frags are not overlapping
                    s1 = data[i,2] < data[j,1]+thresh
                    s2 = data[i,1] > data[j,2]-thresh
                    # f1 in f2
                    # f1 overlap f2 beyond thresh
                    keep[j] = (s1 | s2)
    # decode the original ordering and return data
    data = df.iloc[s[keep]]
    return data

def remove_overlaps(frags, metric='AS', thresh=0):
    '''
    Function to filter out bad overlapping matches in sequences
    
    frags = output from find_RCA_frags()
            must be dataframe with columns [query_id, q_start, q_end, <alignment score>]
    thresh = amount of overlap allowed
    metric = <alignment score> defined by user, default is AS from bwa or minimap2
    
    returns same dataframe with overlapping fragments removed
    '''
    df = frags.reset_index(drop=True)
    
    logging.info('filtering to best aligned fragments')
    f1 = find_overlaps(df[['query_id','q_start','q_end']], thresh=thresh)
    s = f1[f1['cluster_item_count']>1]
    # do work only on overlapping info
    if len(s) > 0:
        f2 = split_overlaps(df.iloc[s.index][['query_id','q_start','q_end',metric]], thresh=thresh)
        f2 = df.iloc[f2.index] # get the split up frags
        f1 = df.iloc[f1[f1['cluster_item_count']==1].index] # get single frags
        df = pd.concat([f1,f2]) # add them together
    # report results
    logging.info('Cleaned out '+str(len(frags) - len(df))+' of '+str(len(frags)))
    return df

def get_best(df, col, metric='AS', stat='idxmax'):
    df=df.reset_index(drop=True)
    idx = df.groupby(by=col).agg({metric:stat}).reset_index()
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
                logging.error('syntax error on spoa call: '+j)
                return -1
    else:
        cmd = ''
        for msa_in in MSA_infiles:
            MSA_outfile = msa_in + '.out'
            f = open(MSA_outfile, 'w') # open file for writing std out
            cmd+= aligner+' '+config+' '+msa_in+' > '+MSA_outfile+' & '
        subprocess.run(cmd[:-3], shell = True) # submit the job

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

def get_pairwise_distance(seq, block=500, output_folder='./pairwise_dist/', workspace='./pw_aligner/', config='-k15 -w10 -PD', symmetric=False, cleanup=False):
    '''
    Function to compute pairwise distance matrix of a list of sequences
    seq = pandas dataframe of a list of sequences with columns [id, sequence]
    block = number of sequences to compare in a block --> break up computation and save it
    output_folder = where output csv files are stored
    workspace = workspace folder for the aligner
    configs = alignment config options for minimap2
    symmetric = evaluates only pairwise alignments of a longer sequence 
                to a shorter one and returns data for a symmetric matrix
                symmetric evaluations take half time
    returns list of csv files where data is stored
    '''
    output_folder = check_dir(output_folder)
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
            logging.info('processing block '+str(i)+','+str(j)+' '+str(count)+'/'+str(N))
            count+=1
            out = run_minimap2(df[j], df[i], config = config, workspace = workspace, 
                               cigar = True, build_index = build_index, use_index = True, cleanup=False)
            build_index = False
            fname = output_folder+'dst_'+str(i)+'_'+str(j)+'.csv.gz'
            out.to_csv(fname, index=False, compression='infer')
            data.append(fname)
    # cleanup the workspace folder
    if cleanup:
        subprocess.run(['rm','-r',workspace])
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

def run_PCA(df, n_comp=2, recenter=False, ofile='', timestamp=True):
    '''
    Runs principle component analysis to compression dimensions 
    df = pandas dataframe with columns [id, d_i, ..., d_n]
    n_comps = components to reduce to
    ofile = output file prefix to save pca data
    timestamp = add timestamp to output file
    '''
    
    # initialize PCA settings
    logging.info('running pca with n_comp = '+str(n_comp))

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
             ofile='', timestamp=True, verbose=0):
    '''
    Run tsne on input data
    df = pandas dataframe with columns [id, v_0, ..., v_i]
    '''
    logging.info('running sklearn tsne with n_comp = '+str(n_comp))
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

def cluster_sample(df_q, df_d=[], N=100, th=0.9, rsize=1, config ='-k15 -w10 -p 0.9 -D', workspace = './clust_sample/'):
    '''
    Find get rough cluster centers for a set of sequence data
    df_q = dataframe of query sequences to search for cluster centers.
    df_d = dataframe of cluster centers with you want to append to
           df_q and df_d must both have at least columns [id, sequence]
    N = number of orthogonal sequences to draw
    th = threshold match_score for excluding sequences which are not similar enough to the sample
         set higher if you want finer spacing between clusters
    rsize = number of samples to draw on each iteration
    workspace = workspace folder for the aligner
    config = config passed to minimap2
    '''
    qry = np.unique(df_q['id'].astype(str))
    if len(qry)==0:
        logging.info('cluster_sample: nothing to cluster')
        return df_d
    # initialize cluster center
    df_c = pd.DataFrame([], columns=['id','sequence'])
    if len(df_d) > 0:
        df_c = df_c.append(df_d[['id','sequence']])
        N+= len(df_c)
    else:
        # pick random sequence to partition the data
        ridx = np.random.permutation(len(qry))
        df_c.append(df_q[df_q['id'].isin(qry[ridx[:rsize]])][['id','sequence']])
    db = df_c
    while N > len(df_c):
        logging.info('cluster_sample: qlen='+str(len(qry))+' progress='+str(len(df_c))+'/'+str(N))
        # align centers to data and look for unaligned samples
        seq = df_q[df_q['id'].isin(qry)]
        out = run_minimap2(seq, db, config=config, cigar=True, workspace=workspace)
        out = out[out['match_score'] > th]
        # make the new block to sequence in
        qry = set(qry) - set(out['query_id'].astype(str)) - set(df_c['id'].astype(str))
        if len(qry)==0:
            break
        # pick random sequence to partition the data
        qry = np.array([i for i in qry])
        ridx = np.random.permutation(len(qry))
        db = df_q[df_q['id'].isin(qry[ridx[:rsize]])][['id','sequence']]
        df_c = df_c.append(db)
    df_c['id'] = ['cluster'+str(i) for i in range(0,len(df_c))]
    return df_c

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

def add_sample_weight(df, weight):
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

def cluster_hierarchical(df, metric='precomputed', linkage='single', thresh=0.5, n_clusters=None):
    '''
    Perform hierarchical clustering on a distance matrix
    '''
    # running the clustering
    logging.info('Running hierarchical clustering')

    # initialize settings
    clust = sklearn.cluster.AgglomerativeClustering(n_clusters = n_clusters, distance_threshold = thresh,
                                                    affinity = metric, linkage = linkage)
    cols = df.columns[df.columns!='id'] # get only non read id values
    clust.fit(df[cols].values)
    
    L = clust.labels_
    for i in np.unique(L):
        if np.sum(L==i) == 1:
            L[L==i] = -1 # set cluster size = 1 as outliers
    return pd.DataFrame(np.transpose([df['id'].values, L]), columns = ['id','cluster_id'])

def cluster_HDBSCAN(df, metric='euclidean', min_cluster_size=100, min_samples=None, n_jobs=4):
    '''
    Do clustering with HDBSCAN
    df = data frame containing columns [id, d_i..., d_n]
    metric = distance metric to use
    min_samples = min_sample parameter to hdbscan
    min_cluster_size = min_cluster_size parameter to hdbscan
    n_jobs = number of parallel processes to run if supported by algorithm, -1 means use all processors
    '''
    # initialize the 
    clust = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples,
                            core_dist_n_jobs=n_jobs, metric=metric)
    
    logging.info('Running HDBSCAN')
    cols = df.columns[df.columns!='id'] # get only non read id values
    clust.fit(df[cols].values)
    
    logging.info('getting and ordering the results')
    outlier = clust.outlier_scores_
    prob = clust.probabilities_
    labels = clust.labels_
    r_id = df['id'].values
    
    logging.info('number of clusters = '+str(np.max(labels)+1))
    logging.info('unclustered = '+str(np.sum(labels==-1)))
    
    s = np.lexsort((outlier, labels))
    df = np.transpose([r_id[s], labels[s], range(0,len(r_id)), prob[s], outlier[s]])
    df = pd.DataFrame(df, columns = ['id', 'cluster_id', 'ordering', 'probability', 'outlier_score'])
    df['cluster_id'] = df['cluster_id'].astype(int)
    df['ordering'] = df['ordering'].astype(int)
    df['probability'] = df['probability'].astype(float)
    df['outlier_score'] = df['outlier_score'].astype(float)
    return df

def cluster_OPTICS(df, metric='euclidean', min_samples=None, min_cluster_size=None, xi=0.05,
                   max_eps=None, cluster_method='dbscan', n_jobs=None, alt_label=False):
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
    '''
    # running the clustering
    logging.info('Running OPTICS')

    # Safety check on number of samples given
    if len(df) < 3:
        logging.warning('Sample size too small to use with OPTICS')
        df['cluster_id'] = range(0,len(df))
        df['ordering'] = range(0,len(df))
        df['reachability'] = -1
        return df[['id','cluster_id','ordering','reachability']]

    cols = df.columns[df.columns!='id'] # get only non read id values
    if max_eps == None:
        z = np.abs(df[cols].values)
        max_eps = (np.max(z)+np.min(z))/2
    logging.info('max_eps = '+str(max_eps))
    
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
        logging.info('clust_OPTICS: iter='+str(i)+' using min_samples='+str(min_samples))
        clust = sklearn.cluster.OPTICS(min_samples=min_samples, min_cluster_size=min_cluster_size, xi=xi,
                            max_eps=max_eps, cluster_method=cluster_method, metric=metric, n_jobs=n_jobs)
        clust.fit(df[cols].values)
        # check labels for cluster number and coverage
        nout = np.sum(clust.labels_ == -1)
        nclust = np.max(clust.labels_)+1
        logging.info('clust_OPTICS: clusters='+str(nclust)+' outliers='+str(nout)+' delta='+str(delta))
        # always keep the change if we get more clusters
        if nclust > prev_nclust or (nclust == prev_nclust and nout <= prev_nout):
            tmp = min_samples
            delta = int((prev_min - min_samples)/2)
            min_samples-= delta
        # if cluster number goes down as we increase min_samples
        elif nclust < prev_nclust:
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
        if delta==0 or min_samples <= 2 or nout==0:
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
        if alt_label and len(reach) > 5:
            logging.info('using alt labeling')
            labels = reachability_alt_label(reach)
        logging.info('n_clusters='+str(np.max(labels)+1)+' n_unclustered='+str(np.sum(labels==-1))+' N='+str(len(df)))
    else:
        logging.info('clust_OPTICS: reachability < 0, please change the min_samples parameter you are using')
    
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

def cluster_Kmeans(df, n_clusters, n_init=10, n_iter=100):
    '''
    Function to perform clustering via k means algorithm
    '''
    # running the clustering
    logging.info('Running kmeans with n_clusters = '+str(n_clusters))
    cols = df.columns[df.columns!='id'] # get only non read id values
    x = df[cols.values]
    
    # initialize settings
    clust = sklearn.cluster.MiniBatchKMeans(n_clusters = n_clusters, n_init = n_init, max_iter = n_iter)
    clust.fit(x)
    
    logging.info('Getting results')
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

def cluster_spoa_merge(df, config='-l 0 -r 2', workspace='./clust_spoa_merge/', batch_size=100, cleanup=True):
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
    workspace = check_dir(workspace)
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
            [r_id, c_id, fname] = df.values[i]
            ftype = fname.split('.')
            if 'fa' in ftype or 'fasta' in ftype:
                df_fa.append(read_fasta(fname))
            elif 'fq' in ftype or 'fastq' in ftype:
                df_fq.append(read_fastq(fname))

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
        logging.info('cluster_spoa_merge: spoa on '+str(i)+'/'+str(N_chunks))
    
    # get consensus
    logging.info('cluster_spoa_merge: reading consensus')
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
