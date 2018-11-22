#!/usr/bin/python
import warnings
warnings.filterwarnings("ignore")
#
import sys
import os
from optparse import OptionParser
import Levenshtein
from multiprocessing import Pool
#
#
opts = OptionParser()
usage = "Trim adapter\nusage: %prog -s source_folder --np 4"
opts = OptionParser(usage=usage, version="%prog 2.1")
opts.add_option("-s", help="Source folder path, should contains {data} folder, where the pair-end sequencing data locate. "
                          +"If you want to use this code to build cell_info.csv file, each fastq file should be names as:"
                          +"type1-001_1.fastq, type1-001_2.fastq, type1-002_1.fastq, type1-002_2.fastq, ... ;"
                          +"type2-001_1.fastq, type2-001_2.fastq, type2-002_1.fastq, type2-002_2.fastq, ... ; etc. "
                          +"{type1, type2, ..., typeN} can be cell-types for your samples, such as {GM, K562, ...}, "
                          +"or you can just use any name you want, but make sure there is no underline(_) or dashline(-) in typeX.")
opts.add_option("--qlen", default=20, help="Query length for adatper trimming, default=20.")
opts.add_option("--aseq", default="CTGTCTCTTATACACATCTGACGCTGCCGACGA", help="Adapter sequence, "
                          +"default=CTGTCTCTTATACACATCTGACGCTGCCGACGA.")
opts.add_option("--np", default=1, help="Number of CPUs used for trimming in parallel, default=1.")
options, arguments = opts.parse_args()
#
#
global query_length, adatper_seq
query_length = options.qlen
adapter_seq = options.aseq
#
#
def mismatch_align(seq1, query_length, read2_rc):
    for s1 in range(len(seq1)-query_length+1, -1, -1):
        temp_read1 = seq1[s1:(s1+query_length)]
        editdist = Levenshtein.distance(temp_read1, read2_rc)
        if editdist<2:
            return s1
    return -1
#
#
def rev_comp_dna(read2_rc):
    temp_read2 = ''
    for i in range(len(read2_rc)-1, -1, -1):
        if (read2_rc[i]=='A') | (read2_rc[i]=='a') :
            temp_read2 += 'T'
        elif (read2_rc[i]=='C') | (read2_rc[i]=='c') :
            temp_read2 += 'G'
        elif (read2_rc[i]=='G') | (read2_rc[i]=='g') :
            temp_read2 += 'C'
        elif (read2_rc[i]=='T') | (read2_rc[i]=='t') :
            temp_read2 += 'A'
        elif read2_rc[i]=='N':
            temp_read2 += 'N'
        else:
            return 'error'
    return temp_read2
#
#
def trim_adapters(fastq):
    cutoff = 50
    fastq1, fastq2 = fastq + '_1.fastq', fastq + '_2.fastq'
    trimed1, trimed2 = fastq + '_1.trim.fastq', fastq + '_2.trim.fastq'
    with open(fastq1) as fa1, open(fastq2) as fa2, open(trimed1, 'w') as out1, open(trimed2, 'w') as out2 :
        nReads, mm0_num_read, mm1_num_read = 0, 0, 0
        while 1:
            seq_header1, seq_header2 = fa1.readline()[:-1], fa2.readline()[:-1]
            seq1, seq2 = fa1.readline()[:-1], fa2.readline()[:-1]
            qual_header1, qual_header2 = fa1.readline()[:-1], fa2.readline()[:-1]
            qual1, qual2 = fa1.readline()[:-1], fa2.readline()[:-1]
            nReads += 1
            if ((not seq_header1) | (not seq_header2) | (not seq1) | (not seq2) |
                (not qual_header1) | (not qual_header2) | (not qual1) | (not qual2)): break
            read2_rc = seq2[:query_length]
            read2_rc = rev_comp_dna(read2_rc)
            s1_pos = -1
            s1_pos_find = seq1.rfind(read2_rc)
            if s1_pos_find > 0 :
                s1_pos = s1_pos_find
                mm0_num_read += 1
            else:
                s1_pos = mismatch_align(seq1, query_length, read2_rc)
                if s1_pos>0: mm1_num_read += 1
            if s1_pos >= 0 :
                seq_len = s1_pos + query_length
                trim_seq1 = seq1[seq_len:]
                adapter_trim_seq = adapter_seq[:len(trim_seq1)]
                if adapter_trim_seq==trim_seq1:
                    seq1 = seq1[:seq_len]
                    seq2 = seq2[:seq_len]
                    qual1 = qual1[:seq_len]
                    qual2 = qual2[:seq_len]
            print >> out1, seq_header1
            print >> out1, seq1[:cutoff]
            print >> out1, qual_header1
            print >> out1, qual1[:cutoff]
            print >> out2, seq_header2
            print >> out2, seq2[:cutoff]
            print >> out2, qual_header2
            print >> out2, qual2[:cutoff]
    return nReads, mm0_num_read, mm1_num_read
#
#
fastqs = [x.split('_')[0] for x in os.listdir(options.s+'/data/') if (x[-6:]=='.fastq')&(x[-11:]!='.trim.fastq')]
fastqs = list(set(fastqs))
fastqs.sort()
pathes = [options.s+'/data/'+x for x in fastqs]
pool = Pool(int(options.np))
read_info = pool.map(trim_adapters, pathes)
pool.close()
pool.join()
#
with open(options.s+'/data/cell_info.csv', 'w') as output:
    print >> output, 'name\tnotes'
    for fastq in fastqs:
        print >> output, fastq + '\t' + '-'.join(fastq.split('-')[:-1])
#
#
#
#
