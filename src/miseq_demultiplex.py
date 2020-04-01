#!/usr/bin/env python
# Matches reads against the database of sequences
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

def write_fasta(file, data):
    '''
    Write list of sequences into a fasta file format.
    
    data =  array with [id, sequence] as the columns
    file = file to write the information
    '''
    f = open(file,'w')
    if len(data.shape) > 1: # case for table
        for i in data:
            text = '>' + i[0] + '\n'
            text+= i[1] + '\n'
            f.write(text)
    else: # case for single row
        text = '>' + data[0] + '\n'
        text+= data[1] + '\n'
        f.write(text)
    f.close()

def get_SAMinfo(read):
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
    
    # get other bowtie info: AS, XS, XN, XM, XO, XG, NM
    # use a dictionary to quickly save the variables instead of if statements
    case = {'AS':0, 'XS':1, 'XN':2, 'XM':3, 'XO':4, 'XG':5, 'NM':6}
    info = [0]*7
    if rname != '*':
        for i in range(11,len(read)):
            val = read[i].split(':') # parse it backwards
            if val[0] in case.keys():
                info[case[val[0]]] = int(val[-1])
    return [rname, qname, orientation, flag, int(pos1), int(pos2), int(mapq), info, cigar]

def parse_SAMfile(infile, outfile, max_RAM_size = 100e6):
    '''
    This function parses the information from a SAM file and saves the info into a csv
    '''
    print('Chunking and parsing ',infile)
    # initialize variables
    past_header = False
    offset = -10
    max_RAM_size = int(max_RAM_size)
    file_size = os.path.getsize(infile)

    # iterate through chunks of the file
    f = open(infile)
    text = f.read(max_RAM_size)
    while text!='':
        print('Reading ',infile,' ', f.tell()/1e6, '/', file_size/1e6)
        # check if sam header is passed
        if past_header == False:
            text = text.split('\n@PG')
            # execute special instruction when end of header is detected in the block
            if len(text) > 1:
                text = text[1]
                temp = text.split('\n')
                if len(temp) > 1:
                    past_header = True
                    text = text[len(temp[0])+1:]

        # only start saving data if header is past
        if past_header:
            # check if lines were not cleanly read for the block
            if text[-1]=='\n':
                text = text.split('\n')[:-1]
                offset = 0
            else:
                text = text.split('\n')
                offset = -len(text[-1])
                text = text[:-1]
            
            data = []
            # process the lines of text
            for line in text:
                out = get_SAMinfo(line)
                qname = out[1]
                info = out[7]
                cigar = out[8]
                data.append([qname, out[2], out[0], out[4], out[5], out[6],
                            info[0], info[1], info[2], info[3], info[4], info[5], info[6], cigar, out[3]])
        
            # toss queries which don't match anything in the database
            #data = np.array(data)
            #data = data[data[:,2]!='*']
            
            # dump the data into a file to save RAM
            out = pd.DataFrame(data, columns=['query_id','orientation','database_id','t_start','t_end','mapq',
                                              'AS','XS','XN','XM','XO','XG','NM','CIGAR','flag'])
            
            # export data
            timestamp = time.strftime('%Y%m%d_%H%M%S', time.localtime())
            fname = outfile.split('.csv')[0]+'_TS_'+timestamp+'.csv'
            print('writing file',fname)
            out.to_csv(fname, index = False)

        # read next line of text
        f.seek(f.tell() + offset)
        text = f.read(max_RAM_size)
        
    # close the file
    f.close()

def run_bwa_on_fq(query, database, folder = './bwa_work/', config = ' mem -a ', build_index = True):
    '''
    This is a wrapper for the bwa aligner.

    query = list of sequences to search for in database
    query must be in pandas dataframe format with columns = ['id','sequence']

    database = database of sequences being search through
    database must be in pandas dataframe format with columns = ['id','sequence']

    build_index = determines if a burrow-wheeler transform index should be built for the database
    folder = defines where input and output files for bowtie2 resides
    config = defines the config options to pass to bowtie2
    '''
    
    start = time.time()
    # Define locations of input and output files for bowtie
    qfile = query
    dbfile = folder+'database.fa'
    outfile = folder+'results.sam'

    # Make the working folder if it does not exists
    if os.path.exists(folder) == False:
        print('Making directory ',folder)
        subprocess.call('mkdir '+folder, shell=True)
    
    # Only build index when we need to. Index building takes a lot of time
    if build_index == True:
        if len(database.columns) == 3:
            write_fasta(dbfile, database[['id','sequence','quality']].values)
        else:
            write_fasta(dbfile, database[['id','sequence']].values)
        cmd = 'bwa index '+dbfile
        print(cmd)
        subprocess.call(cmd, shell = True)
        print('bwa elapse = ',time.time()-start,'s')
    
    # call the aligner
    cmd = 'bwa '+config+' '+dbfile+' '+qfile+' > '+outfile
    print(cmd)
    subprocess.call(cmd, shell = True)
    print('bwa elapse = ',time.time()-start,'s')
    return outfile

def merge_paired_reads(outheader, ofile1, ofile2):
    file1 = glob.glob(outheader+'_'+ofile1+'*.csv')
    file2 = glob.glob(outheader+'_'+ofile2+'*.csv')
    # load files
    df1 = []
    for f in file1:
        print('loading ', f)
        print('finding best matches for primer tags')
        metric = 'mapq'
        out = pd.read_csv(f)
        #x = out.groupby('query_id',as_index=False).agg({metric:'idxmax'})
        #out = out.iloc[x[metric].values]
        df1.append(out)

    df2 = []
    for f in file2:
        print('loading ', f)
        print('finding best matches for primer tags')
        metric = 'mapq'
        out = pd.read_csv(f)
        #x = out.groupby('query_id',as_index=False).agg({metric:'idxmax'})
        #out = out.iloc[x[metric].values]
        df2.append(out)
    df1 = pd.concat(df1)
    df2 = pd.concat(df2)
    
    # rename column headers
    col1 = ['query_id']
    col2 = ['query_id']
    for i in range(1,len(df1.columns)):
        col1.append(df1.columns[i]+'1')
        col2.append(df2.columns[i]+'2')
    df1.columns = col1
    df2.columns = col2

    # merge the files
    print('merging paired end reads in csv files')
    df = df1.merge(df2, how = 'outer')
    
    # print some stats
    print('stats for',outheader)
    print('total reads ', len(df))
    unpaired1 = np.sum(df['orientation1'].isna())
    unpaired2 = np.sum(df['orientation2'].isna())
    print('unmatched r1 ', unpaired1)
    print('unmatched r2 ', unpaired2)
    print('paired reads ', len(df) - unpaired1 - unpaired2)
    
    opp = np.sum((df['orientation1']=='+') & (df['orientation2']=='+'))
    omm = np.sum((df['orientation1']=='-') & (df['orientation2']=='-'))
    opm = np.sum((df['orientation1']=='+') & (df['orientation2']=='-'))
    omp = np.sum((df['orientation1']=='-') & (df['orientation2']=='+'))
    
    print('orientation ++ reads        = ', opp)
    print('orientation -- reads        = ', omm)
    print('orientation +- reads        = ', opm)
    print('orientation -+ reads        = ', omp)
    print('orientation -+ and +- reads = ', opm + omp)
    
    # get frequency of read counts for each primer combination
    x1 = df['database_id1'].dropna().values
    x2 = df['database_id2'].dropna().values
    
    pmr = np.concatenate((x1, x2))
    pmr = np.unique(pmr)
    data = []
    for i in range(0, len(pmr)):
        for j in range(i, len(pmr)):
            print('Summing for ',pmr[i],pmr[j])
            s1 = (df['orientation1']=='+')&(df['orientation2']=='+')
            s2 = (df['orientation1']=='-')&(df['orientation2']=='-')
            s3 = (df['orientation1']=='+')&(df['orientation2']=='-')
            s4 = (df['orientation1']=='-')&(df['orientation2']=='+')
            s5 = (df['database_id1']==pmr[i])&(df['database_id2']==pmr[j])
            s6 = (df['database_id1']==pmr[j])&(df['database_id2']==pmr[i])
            x1 = np.sum(s1&(s5 | s6))
            x2 = np.sum(s2&(s5 | s6))
            x3 = np.sum(s3&(s5 | s6))
            x4 = np.sum(s4&(s5 | s6))
            x5 = np.sum(s5 | s6)
            if x5 > 0:
                data.append([pmr[i], pmr[j], x1, x2, x3, x4, x5])

    data = pd.DataFrame(data, columns = ['primer1','primer2',
                                         'orientation ++','orientation --',
                                         'orientation +-','orientation -+','total'])
    data.to_csv('stats.csv', index=False)
    
    # export the data
    print('writing to csv file')
    df.to_csv('merged_'+ofile1+'.csv', index=False)

# degenerate primer fix
def get_nondegen_primer(seq):
    '''
    Get a new list of primers without degenerate sequences
    
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

def print_help(argv):
    print('Usage: '+argv[0]+' -p <primer_file> -i <fastq_files> -o <out_dir>')
    print('-h              print help')
    print('-p              <primer_file> = name of primer csv file, such as ./primers.csv')
    print('-i              <fastq_files> = fastq files, such as ./data/*.fq.gz')
    print('-o              <out_file> = name of output file header for csv files, such as output')
    print('-of             <out_dir> = name of output folder for csv files, output_folder/')
    print('-mem            <int> = MB of RAM to use when reading sam file output from BWA, such as 100')
    sys.exit(1)

def main(argv):
    if len(argv) < 2:
        print_help(argv)
        sys.exit(1)
    print(argv)
    
    # Default settings
    primer_file = 'primers.csv'
    fq_file = './data/*.fastq*'
    out_file = 'output'
    out_folder = 'output_folder/'
    qfile = 'primers.fa'
    max_RAM_size = int(100e6) # use at most 100MB of RAM when reading the sam file
    
    if len(argv) < 2:
        print_help(argv)
    else:
        for i in range(1,len(argv)):
            if argv[i] == '-h' or argv[i] == '--help': print_help(argv)
            elif argv[i] == '-p': primer_file = argv[i+1]
            elif argv[i] == '-i':
                fq_file = []
                j = 1
                while i+j < len(argv):
                    if argv[i+j][0]=='-':
                        break
                    fq_file.append(argv[i+j])
                    j+=1
            elif argv[i] == '-o': out_file = argv[i+1]
            elif argv[i] == '-of': out_folder = argv[i+1]
            elif argv[i] == '-mem': max_RAM_size = int(argv[i+1])*1e6
    
    print('primer file   = ',primer_file)
    print('fastq files   = ',fq_file)
    print('output header = ',out_file)
    print('output file   = ',out_folder)
    print('RAM to use    = ',max_RAM_size/1e6,'MB')
    
    start = time.time()
    
    # Make the output folder if it does not exists
    if os.path.exists(out_folder) == False:
        print('Making directory ',out_folder)
        subprocess.call('mkdir '+out_folder, shell=True)
    
    # read primer info
    primers = pd.read_csv(primer_file)
    primers['sequence'] = primers['barcode']+primers['comment']
    primers['id'] = 'pID_'+primers['ID']+'_uID_'+primers['Unique_ID']
      
    # get non-degenerate primer combos
    plist = []
    for i in range(0,len(primers)):
        x = get_nondegen_primer(primers.iloc[i]['sequence'])
        x = pd.DataFrame(x, columns = ['sequence'])
        x['id'] = primers.iloc[i]['id']
        plist.append(x)
    plist = pd.concat(plist)
    
    # Search for the primers in miSeq fastq
    '''
    config = 'mem'
    for query in fq_file:
        # run bwa aligner
        sam_file = run_bwa_on_fq(query, plist[['id','sequence']], folder = './aligner_work/', config = config)
        # extract the SAM file information
        q_file = query.split('/')[-1].split('.')[0]
        parse_SAMfile(sam_file, out_folder+out_file+'_'+q_file, max_RAM_size)

    # clear the work folder
    subprocess.call('rm -r aligner_work', shell=True)
    '''
    
    # merge paired end info together
    for i in np.arange(0, len(fq_file), 2):
        q1_file = fq_file[i].split('/')[-1].split('.')[0]
        q2_file = fq_file[i+1].split('/')[-1].split('.')[0]
        merge_paired_reads(out_folder+out_file, q1_file, q2_file)

if __name__ == "__main__":
    main(sys.argv)