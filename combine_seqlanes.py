#!/usr/bin/env python

'''
combine_seqlanes.py: Combines data from multiple lanes into a single .fq.gz file for further analysis.
This is often needed for sequencers like the NextSeq.
Example usage:
python combine_seqlanes.py indir outdir -common_suffix _S*_L00*_R*_001.fastq.gz
# python combine_seqlanes.py indir outdir -common_suffix '_S*_L00*_R*_001.fastq.gz'
Requires Python >= 3.5
'''

__author__ = 'Mary Thompson'

import sys
import os
import re
import logging
import shutil
import glob
import argparse
from collections import defaultdict

def concatenate_files(indir, outdir, r1_mark='R1', r2_mark='R2', common_suffix='_S*_L00*_R*_001.fastq.gz', 
                      lane_marker='L00*', nested_dirs=True):
    '''
    Find files matching patterns, sort lanes, and combine R1 and R2 files.
    All files must be in the same input directory.
    File R2 is optional and its absence indicates single-end epxeriments.
    nested_dirs = search subdirectories for files matching pattern
    '''
    os.makedirs(outdir, exist_ok=True)
    logging.basicConfig(level=logging.DEBUG, filename=os.path.join(outdir, 'combine_nextseq.log'),
    filemode='w', format='%(message)s')

    common_suffix_re = common_suffix.replace('*', '.')
    lane_marker_re = lane_marker.replace('*', '.')
    glob_pattern = f'{indir}/**/*{common_suffix}'
    files = glob.glob(glob_pattern, recursive=nested_dirs)

    # Make list of dict of sample: [R1 files, R2_files]
    file_info = defaultdict(lambda: [[], []])
    for file in files:
        sample_name = os.path.basename(re.split(common_suffix_re, file)[0])
        if re.search(r1_mark, file):
            file_info[sample_name][0].append(file)
        elif re.search(r2_mark, file):
            file_info[sample_name][1].append(file)

    # Sort the R1 and R2 files
    # Then copy file contents to the outfile
    outdict = {}
    for s in file_info:
        R2_exists = file_info[s][1] != []
        logging.info(f'processing sample {s}')
        for i in file_info[s]:
            i.sort(key = lambda x: int(re.search(lane_marker_re, x).group(0)[-1]))

        R1_outfile = os.path.join(outdir, f'{s}_R1.fq.gz')
        outdict[R1_outfile] = file_info[s][0]
        outnames = [R1_outfile]
        logging.info(f'combining files {file_info[s][0]} into {R1_outfile}')

        if R2_exists:
            R2_outfile = os.path.join(outdir, f'{s}_R2.fq.gz')
            outdict[R2_outfile] = file_info[s][1]
            outnames = [R1_outfile, R2_outfile]
            logging.info(f'combining files {file_info[s][1]} into {R2_outfile}')

        for outname, infiles in zip(outnames, file_info[s]):
            with open(outname, 'wb') as outfile:
                for lane in infiles:
                    with open(lane, 'rb') as f:
                        shutil.copyfileobj(f, outfile)
        logging.info(f'processed sample {s}')
        logging.info('\n')
    return outdict

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('indir')
    parser.add_argument('outdir')
    parser.add_argument('-r1_mark', default='R1')
    parser.add_argument('-r2_mark', default='R2')
    parser.add_argument('-common_suffix', default='_S*_L00*_R*_001.fastq.gz')
    parser.add_argument('-lane_marker', default='L00*')
    parser.add_argument('--no_nested_dirs', action='store_true', default=False)
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    nested_dirs = not args.no_nested_dirs
    outdict = concatenate_files(args.indir, args.outdir, r1_mark=args.r1_mark, r2_mark=args.r2_mark, 
                      common_suffix=args.common_suffix, lane_marker=args.lane_marker, 
                      nested_dirs=nested_dirs)