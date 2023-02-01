import combine_seqlanes
import random
import string
import gzip
import os
import shutil
import subprocess
from collections import defaultdict

def gzip_linecount(infile):
    '''Count the lines in a gzipped file.'''
    # i needs to be initialized be avoid pytest error. If no lines, than -1 + 1 = 0
    i = -1
    with gzip.open(infile, 'rb') as f:
        for i, l in enumerate(f):
            pass
    n_lines = i + 1
    return n_lines

def count_lines_gzip(outdict):
    '''
    Check that the number of reads in the output file is the same as the number of reads in each 
    input file.
    '''
    for i in outdict:
        n_lines_original = 0
        for j in outdict[i]:
            n_lines_original += gzip_linecount(j)
        n_lines_out = gzip_linecount(i)
        assert n_lines_out == n_lines_original, f'Expected {n_lines_original}, Result {n_lines_out}, Outfile: {i}'

def check_read_order(outdict):
    '''
    Check that the order of the reads is the same between the R1 and R2 output files (by checking read ID)
    '''
    suffix_r1 = '_R1.fq.gz'
    suffix_r2 = '_R2.fq.gz'
    paired_files = defaultdict(list)
    for i in outdict.keys():
        if len(i.split(suffix_r1)) > 1:
            fname = i.split(suffix_r1)[0]
        elif len(i.split(suffix_r2)) > 1:
            fname = i.split(suffix_r2)[0]
        else:
            pass
        paired_files[fname].append(i)

    for pair in paired_files:
        with gzip.open(paired_files[pair][0], 'rt') as f:
            with gzip.open(paired_files[pair][1], 'rt') as g:
                f1_ids = [line.strip('\n').split(' ')[0] for line in f if line.startswith('@')]
                f2_ids = [line.strip('\n').split(' ')[0] for line in g if line.startswith('@')]
    
    # Check if the IDs are the same and in the same order
    assert f1_ids == f2_ids

def check_read_content(outdict=None, file1=None, file2=None):
    '''
    Check that read ids in specified files are equivalent, but not necessarily in the same order.
    Provide either outdict or file1 and file2
    Check that all the reads contained in the file1 are also in file2
    If outdict is provided, check that reads contained in each outfile are equal to all reads in the 
    input files
    '''
    if outdict is not None:
        pass
    elif (file1 is not None) and (file2 is not None):
        with gzip.open(file1, 'rt') as f:
            with gzip.open(file2, 'rt') as g:
                f1_ids = {line.strip('\n').split(' ')[0] for line in f if line.startswith('@')}
                f2_ids = {line.strip('\n').split(' ')[0] for line in g if line.startswith('@')}
                # check symmetrical!
                assert f1_ids.difference(f2_ids) == {}
    else:
        raise FileNotFoundError('Both file1 and file2 are required if outdict is not passed')

def compare_to_cat(d):
    for sample in d:
        outname = sample
        cat_cmd = ['cat', *d[sample]]
        newname = outname.rstrip('.fq.gz') + '_cat' + '.fq.gz'
        with open(newname, 'w') as f:
            subprocess.run(cat_cmd, stdout=f)
        
        same_file = False
        if os.path.getsize(outname) == os.path.getsize(newname):
            if gzip.open(outname,'r').read() == gzip.open(newname,'r').read():
                same_file = True
        assert same_file

class FakeFastqgz:

    nts = ['A', 'T', 'C', 'G']

    def __init__(self, seqlen=20, n_seqs=100, outname='test_fq', matched_file=None):
        self.seqlen = seqlen
        self.outname = outname
        self.matched_file = matched_file
        # passing a matched file will override default n_seqs
        if self.matched_file is not None:
            self.n_seqs = int(gzip_linecount(self.matched_file)/4)
        else:
            self.n_seqs = n_seqs

    def write(self, output_gz=True):
        '''
        Write fq file. If matched file, then write a file with the same seqids but different seq and qual data
        '''
        self.filename = f'{self.outname}.fq'
        if self.matched_file is not None:
            with gzip.open(self.matched_file, 'rt') as f:
                seq_ids = [i.strip('\n') for i in f if i.startswith('@')]
        else:
            seq_ids = []
            for i in range(self.n_seqs):
                seq_ids.append('@' + ''.join(random.choices(string.ascii_uppercase + string.digits, k=self.seqlen)))

        with open(f'{self.outname}.fq', 'w') as f:
            for i in range(0, self.n_seqs):
                seqid = seq_ids[i]
                seq = ''.join(random.choices(FakeFastqgz.nts, k=self.seqlen))
                qstring = 'E'*self.seqlen
                entry = '\n'.join([seqid, seq, '+', qstring])
                f.write(f'{entry}\n')
        if output_gz:
            self.gzip_file()
            os.remove(self.filename)
            self.filename = f'{self.outname}.fq.gz'
    
    def gzip_file(self):
        outfile = f'{self.outname}.fq.gz'
        with open(f'{self.outname}.fq', 'rb') as f_in:
            with gzip.open(f'{self.outname}.fq.gz', 'wb') as f_out:
                f_out.writelines(f_in)

class TestCombineNextSeq:

    indir = 'testdirs'
    fakefq = 'fakefastq'

    def setup_class(self, n=4):
        '''
        Make three test dirs:
        paired = for paired-end libraries
        single = for single-end libraries
        nested = paired-end libraries in a nested directory structure
        '''
        # Make three output sets: paired, single, nested
        dir_dict = {
            'paired_in': f'{TestCombineNextSeq.indir}/paired_in', 'paired_out': f'{TestCombineNextSeq.indir}/paired_out', 
            'single_in': f'{TestCombineNextSeq.indir}/single_in', 'single_out': f'{TestCombineNextSeq.indir}/single_out', 
            'nested_in': f'{TestCombineNextSeq.indir}/nested_in', 'nested_out': f'{TestCombineNextSeq.indir}/nested_out',
        }

        self.dir_dict = dir_dict
        for dir in dir_dict:
            os.makedirs(dir_dict[dir], exist_ok=True)

        # Make fq files to put in each directory
        self.fq_files = []
      
        # 1) paired
        for i in range(0, n):
            fq1 = FakeFastqgz(outname=os.path.join(dir_dict['paired_in'], f'{TestCombineNextSeq.fakefq}_L{i}_R1'))
            fq1.write()
            fq2 = FakeFastqgz(outname=os.path.join(dir_dict['paired_in'], f'{TestCombineNextSeq.fakefq}_L{i}_R2'), matched_file=fq1.filename)
            fq2.write()
            self.fq_files.append(fq1.filename)
            self.fq_files.append(fq2.filename)

        # 2) single
        for file in os.listdir(dir_dict['paired_in']):
            if 'R1' in file:
                file = shutil.copy(os.path.join(dir_dict['paired_in'], file), dir_dict['single_in'])
                self.fq_files.append(file)

        # 3) nested
        dir1 = os.path.join(dir_dict['nested_in'], f'{os.path.basename(dir_dict["nested_in"])}_1')
        dir2 = os.path.join(dir_dict['nested_in'], f'{os.path.basename(dir_dict["nested_in"])}_2')
        os.makedirs(dir1, exist_ok=True)
        os.makedirs(dir2, exist_ok=True)

        infiles = os.listdir(dir_dict['paired_in'])
        for i in range(len(infiles)):
            this_file = os.path.join(dir_dict['paired_in'], infiles[i])
            if i % 2:
                file = shutil.copy(this_file, dir1)
            else:
                file = shutil.copy(this_file, dir2)
            self.fq_files.append(file)
            i += 1

    def teardown_class(self):
        for i in self.fq_files:
            os.remove(i)

    def test_paired_end(self):
        paired_outdict = combine_seqlanes.concatenate_files(self.dir_dict['paired_in'], self.dir_dict['paired_out'], 
                                                          common_suffix='_L*_R*.fq.gz', lane_marker='L*')
        count_lines_gzip(paired_outdict)
        check_read_order(paired_outdict)
    
    def test_single_end(self):
        single_outdict = combine_seqlanes.concatenate_files(self.dir_dict['single_in'], self.dir_dict['single_out'], 
                                                          common_suffix='_L*_R*.fq.gz', lane_marker='L*')
        count_lines_gzip(single_outdict)
        # cat the files and check if they are the same since cat was used previously
        compare_to_cat(single_outdict)

    def test_nested(self):
        nested_outdict = combine_seqlanes.concatenate_files(self.dir_dict['nested_in'], self.dir_dict['nested_out'], 
                                                          common_suffix='_L*_R*.fq.gz', lane_marker='L*')
        count_lines_gzip(nested_outdict)
        check_read_order(nested_outdict)

    # More tests:
    # works on nested directories
    # works on single and paired-end directories
    # reads are in the same order in files R1 and R2 for paired end directories
    # result is the same as using the cat command