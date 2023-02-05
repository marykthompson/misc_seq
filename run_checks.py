'''
Test that files in one directory with the same names match those in another directory
with the same names.
'''

import argparse
import os
from test_combine_seqlanes import check_read_content, check_read_order

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('dir1')
    parser.add_argument('-dir2')
    parser.add_argument('-file_suffix', default='.fq.gz')
    parser.add_argument('--check_read_content', action='store_true', default=False)
    parser.add_argument('--check_read_order', action='store_true', default=False)
    args = parser.parse_args()
    dir1_files = [i for i in os.listdir(args.dir1) if i.endswith(args.file_suffix)]
    dir2_files = [i for i in os.listdir(args.dir2) if i.endswith(args.file_suffix)]
    pairs = []
    for i,f in enumerate(dir1_files):
        if f in dir2_files:
            match = dir2_files.index(f)
            pairs.append((os.path.join(args.dir1, dir1_files[i]), os.path.join(args.dir2, dir2_files[match])))
        else:
            print(f'no match for {f}')

    if args.check_read_content:
        for p in pairs:
            check_read_content(file1=p[0], file2=p[1])
    
    if args.check_read_order:
        for p in pairs:
            check_read_order(file1=p[0], file2=p[1])