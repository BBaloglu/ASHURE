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

# import useful tools
import ONT_tools

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

def print_help(argv):
    print('Usage: '+argv[0]+' [options] -i miseq_pair1.fq miseq_pair2.fq -o output.csv')
    print('''
    Options:
    -h or --help        print help
    
    -i2 <file1>         miseq fastq pair1
        <file2>         miseq fastq pair2
    
    -p <file>           csv file of forward and reverse primer pairs used in the run
    
    -d <file>           file containing reference database sequences
                        file format can be .fasta, .fa, .fastq, .fq, .csv, or .csv.gz
    
    -o <file>           where OTU and count information is written
    
    -um <folder>        folder of fastq sequences which did not map to reference database
    
    -log <file>         log file
    ''')
    sys.exit(1)

def main(argv):
    if len(sys.argv)<2:
        print_help(argv)

    # Important file and folder defaults
    fastq_pair1 = ''
    fastq_pair2 = ''
    primer_file = ''
    db_file = ''
    out_file = 'output.csv'
    unmapped_folder = 'unmapped_sequences/'
    log_file = 'log.txt'
    
    # Other work folders
    minimap2_folder = './minimap2_work/'

    # parse some arguments from user
    print(argv)
    ONT_tools.append_to_log('Running the pipeline', log_file)
    
    for i in range(0,len(argv)):
        if argv[i] == '-h' or argv[i] == '--help': print_help(argv)
        elif argv[i] == '-i2':
            fastq_pair1 = argv[i+1]
            fastq_pair2 = argv[i+2]
        elif argv[i] == '-p': primer_file = argv[i+1]
        elif argv[i] == '-o': out_file = argv[i+1]
        elif argv[i] == '-d': db_file = argv[i+1]
        elif argv[i] == '-um': unmapped_folder = argv[i+1]

    # if nothing is select, then default is to run everything

    # Print run settings
    ONT_tools.append_to_log('fastq_folder     = '+fastq_folder, log_file)
    ONT_tools.append_to_log('primer_file      = '+primer_file, log_file)
    ONT_tools.append_to_log('db_file          = '+db_file, log_file)
    ONT_tools.append_to_log('log_file         ='+log_file, log_file)
    
    # load read info
    
    # map reads to existing database, generate list of mapped and unmapped
    if db_file=='': # only run if db_file was given
        ONT_tools.append_to_log('no reference database file given.', log_file)
    else:
        ONT_tools.append_to_log('mapping reads to reference database', log_file)
    
    # merge fwd and rev reads in miseq
    
    # cluster and map samples to unmerged files
    
    # map primers to reads
    
if __name__ == "__main__":
    main(sys.argv)
