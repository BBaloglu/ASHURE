#!/usr/bin/env python
# Matches reads against the database of sequences
# bilgenurb@gmail.com

# Loading libraries
import numpy as np
import pandas as pd
# For interfacing with the file system
import subprocess
import os
import time
import argparse
# Loading the custom libraries
import bilge_pype as bpy

def load_file(filename, run_info=None):
    '''
    Checks the file type and load the data
    '''
    file_format = filename.split('.')
    # check if it is a fasta file
    if ('fasta' in file_format) or ('fa' in file_format):
        A = bpy.read_fasta(filename)
    # check if it is a fastq file
    elif ('fastq' in file_format) or ('fq' in file_format):
        A = bpy.read_fastq(filename)
    # assume it is a csv file otherwise
    else:
        A = pd.read_csv(filename)
        if 'consensus' in A.columns:
            A = A.rename(columns={'consensus':'sequence'})
        elif run_info in A.columns:
            # filter for the run info
            A = A[A[run_info] == 'T']
            A = A.rename(columns={'HapID_Bilge':'id','Hap_Sequence':'sequence'})
    # replace ' ' with '_' in sequence id otherwise minimap2 or bowtie2 will complain
    A['id'] = [i.replace(' ','_') for i in A['id']]
    print(len(A),'sequences loaded from ',filename)
    return A

def get_best_matches(df, metric):
    '''
    Finds the best match and also annotates information about the second best match
    '''
    # label primary or secondary matches
    df = df.sort_values(by = ['query_id',metric])
    x = df[['query_id',metric]].values
    k = 0
    rank = -np.ones(len(df))
    for i in range(0,len(df)):
        rank[i] = k
        # if end of list is not reached
        if i+1 < len(df) and x[i,0]==x[i+1,0]:
            k+=1
        else:
            k=0
    df['rank'] = rank
    return df

def main():
    parser = argparse.ArgumentParser(description='Performs alignment of query sequences against database sequences',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', dest='qfile', type=str, required=True,
                        help='''Input query file. Valid formats include .fasta, .fa, .fastq, .fq, or .csv
                                If the file is a csv, it must have columns [id, sequence].
                                If it is from a sequence database, specify the run_info tag to filter for it''')
    parser.add_argument('-o', dest='acc_file', type=str, required=True,
                        help='''Output file for alignment information. This will be in csv format''')
    parser.add_argument('-d', dest='dfile', type=str, default='./species_info/Nanopore_ref_db_all_runs_info.csv',
                        help='''Database file that query sequences will be aligned against.
                                Valid formats include .fasta, .fa, .fastq, .fq or .csv''')
    parser.add_argument('-r', dest='run_info', type=str, default='Nanopore_run_02Feb19',
                        help='''run_info tag used to filter a custom reference database if you have it''')
    parser.add_argument('-a', dest='aligner', type=str, default='minimap2',
                        help='''Aligner to use. Valid options: bwa, bowtie2, or minimap2''')
    parser.add_argument('-all', action='store_true', help='''Report all alignments or not.''')
    parser.add_argument('-wf', dest='workspace', type=str, default='./aligner_folder/',
                        help='''Workspace folder for the aligner''')
    parser.add_argument('-bin', dest='bin_path', type=str, help='path to the aligner')

    args = parser.parse_args()  
    print('query file       = ', args.qfile)
    print('database file    = ', args.dfile)
    print('run_info         = ', args.run_info)
    print('output file      = ', args.acc_file)
    workspace = bpy.check_dir(args.workspace, overwrite=False)
    print('aligner folder   = ', workspace)
    print('all alignments   = ', args.all)
    start = time.time()
    # loading the data
    print('Loading query file:', args.qfile)
    A = load_file(args.qfile)
    print('Loading database file:', args.dfile)
    B = load_file(args.dfile, args.run_info)
    # Run the aligner
    if args.aligner=='bowtie2':
        if args.bin_path!=None:
            bpy._bowtie2 = args.bin_path
            bpy.check_toolchain()
        print('Running bowtie2')
        configs = '--very-sensitive-local'
        data = bpy.run_bowtie2(A[['id','sequence']], B[['id','sequence']], workspace, configs)
        metric = 'AS'
    elif args.aligner=='bwa':
        if args.bin_path!=None:
            bpy._bwa = args.bin_path
            bpy.check_toolchain()
        print('Running bwa')
        configs = 'mem -a -L 100'
        data = bpy.run_bwa(A[['id','sequence']], B[['id','sequence']], workspace, configs)
        metric = 'AS'
    elif args.aligner=='minimap2':
        if args.bin_path!=None:
            bpy._bowtie2 = args.bin_path
            bpy.check_toolchain()
        print('Running minimap2')
        configs = '-k8 -w1 --score-N 0'
        data=bpy.run_minimap2(A[['id','sequence']], B[['id','sequence']], workspace, configs, cigar=True, build_index=True, use_index=True)
        metric = 'AS'
    subprocess.run(['rm','-r',workspace])
    print('Found ',len(data[data['database_id']!='*']),' matching sequences')
    # find best matches
    print('Finding best matches')
    if args.all:
        data = get_best_matches(data, metric)
    elif args.aligner=='minimap2':
        data = bpy.remove_overlaps(data, metric='AS', thresh=20)
        x = data.groupby(by=['query_id']).agg({'query_id':'count'}).rename(columns={'query_id':'count'}).reset_index()
        # flag chimeras
        x['chimeric'] = x['count'] > 1
        print('Found',np.sum(x['chimeric']),'chimeric sequences')
        data = data.merge(x[['query_id','chimeric']], on='query_id', how='left')
    else:
        data = bpy.get_best(data,['query_id'],metric='AS',stat='idxmax')
    print('Found',len(data),'matches')
    # print results summary
    print(len(A['id']) - len(np.unique(data['query_id'])), 'unmatched')
    # save the data
    print('Saving data to',args.acc_file)
    data = data.sort_values(by=['database_id','query_id'])
    data.to_csv(args.acc_file, index=False, compression='infer')
    print('Elapse=',time.time()-start)

if __name__ == "__main__":
    main()
