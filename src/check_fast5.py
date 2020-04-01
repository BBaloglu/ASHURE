#!/usr/bin/env python
# bilgenurb@gmail.com

import argparse
import glob
import subprocess

import numpy as np
import pandas as pd
import h5py
import bilge_pype as bpy

def get_rid(fname):
    ftype = fname.split('.')
    if 'fast5' in ftype:
        f1 = h5py.File(fname1, 'r')
        s1 = set(f1.keys())
        f1.close()
        return s1
    elif 'fastq' in ftype or 'fq' in ftype:
        df = bpy.read_fastq(fname)
        return set(df['id'].astype(str).values)

parser = argparse.ArgumentParser(description='Checks contents of two folders for fast5 file duplication',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-q', dest='files1', required=True, nargs='+', help='only fast5 or fastq files in folder1')
parser.add_argument('-d', dest='files2', required=True, nargs='+', help='only fast5 or fastq files in folder2')
parser.add_argument('-o', dest='outfile', type=str, default='out.csv.gz', help='Output file')
parser.add_argument('-all', action='store_true', help='Compare all files to each other')
parser.add_argument('-remove', action='store_true', help='Remove files listed in -d')
parser.add_argument('-i', dest='input', help='Input csv file to use for removing files')
args = parser.parse_args()

# run fast5 file check
if args.input==None:
    data = []
    for fname1 in args.files1:
        print('examining fname1=',fname1)
        r1 = get_rid(fname1)
        for i in range(0, len(args.files2)):
            fname2 = args.files2[i]
            r2 = get_rid(fname2)
            isec = len(r1.intersection(r2))
            if isec > 0 and fname1!=fname2:
                print('r1 ',isec,'/',len(r1),' r2 ',isec,'/',len(r2))
                print('reads intersecting between', fname1, fname2)
                data.append([fname1, len(r1), len(r1-r2), fname2, len(r2), len(r2-r1), isec])
                if args.all==False:
                    del args.files2[i]
                    break
    df = pd.DataFrame(data, columns=['file1','len1','len1_only','file2','len2','len2_only','intersection'])
    df.to_csv(args.outfile, compression='infer', index=False)
else:
    df = pd.read_csv(args.input)

if args.remove:
    for f1,f2,isec,r1,r2 in df[['file1','file2','intersection','len1','len2']].values:
        print('r1 ',isec,'/',r1,' r2 ',isec,'/',r2)
        subprocess.call(['ls','-l',f1])
        subprocess.call(['ls','-l',f2])
        subprocess.call(['rm','-i',f1])

