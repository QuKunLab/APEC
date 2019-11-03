#!/usr/bin/python
#
import warnings
warnings.filterwarnings("ignore")
#
import os,numpy,sys,subprocess
from optparse import OptionParser
import subroutines
#
#
opts = OptionParser()
usage = "Call peaks\nusage: %prog -s project --blist blacklist.bed --fa genome_chr.fa --tss tssFile --logq 3"
opts = OptionParser(usage=usage, version="%prog 1.0")
opts.add_option("-s", help="The project folder.")
opts.add_option("--picard", default="../reference/picard.jar",
                help="The picard.jar file path, default=../reference/picard.jar")
opts.add_option("--blist", default='../reference/hg19_blacklist.JDB.bed',
                help="Blacklist.bed, default=../reference/hg19_blacklist.JDB.bed")
opts.add_option("--fa", default='../reference/hg19_chr.fa',
                help="Genome_chr.fa, default=../reference/hg19_chr.fa")
opts.add_option('--tss', default='../reference/hg19_refseq_genes_TSS.txt',
                help='TSS file, default=../reference/hg19_refseq_genes_TSS.txt')
opts.add_option('--ref', default='hg19', help='Name of genome reference, default=hg19')
opts.add_option('--logq', default='3',
                help='Threshold of -log(p-value) for top peaks, default=3.')
options, arguments = opts.parse_args()
#
workspace_folder = options.s + '/work/'
peak_folder = options.s + '/peak/'
genome_fasta = options.fa
tssFile = options.tss
if not os.path.exists(peak_folder): subprocess.check_call('mkdir ' + peak_folder, shell=True)
#
#
print('!!!!!!  merge all marked bam files  !!!!!!')
bam_folder = [x for x in os.listdir(workspace_folder)]
bam_folder.sort()
print('cells number:', len(bam_folder))
marked_bam = []
#merged_raw = peak_folder + 'mergeAll.raw.bam'
merged_bam = peak_folder + 'mergeAll.bam'
for folder in bam_folder:
    path = workspace_folder + folder + '/'
    if len(folder.split('.'))<=1:
        marked_bam.extend([path + x for x in os.listdir(path) if x[-10:]=='marked.bam'])
if len(marked_bam)<=1000:
    marked_bam = ' '.join(marked_bam)
    subprocess.check_call('samtools merge -f ' + merged_bam + ' ' + marked_bam, shell=True)
else:
    n_batch = len(marked_bam)//1000 + 1
    temps = []
    for i_batch in range(0, n_batch):
        temp_bam = peak_folder+'temp_'+str(i_batch)+'.bam'
        temps.append(temp_bam)
        start, end = i_batch*1000, min((i_batch+1)*1000, len(marked_bam))
        marked = ' '.join(marked_bam[start:end])
        subprocess.check_call('samtools merge -f ' + temp_bam + ' ' + marked, shell=True)
        subprocess.check_call('samtools index ' + temp_bam, shell=True)
    all_temp = ' '.join(temps)
    subprocess.check_call('samtools merge -f ' + merged_bam + ' ' + all_temp, shell=True)
#
subprocess.check_call('samtools index ' + merged_bam, shell=True)
print('!!!!!!  merge done  !!!!!!')
#
hist_log = peak_folder + 'mergeAll.hist.log'
hist_pdf = peak_folder + 'mergeAll.hist.pdf'
subprocess.check_call('java -XX:+UseSerialGC -Xmx1g -jar '+options.picard+' CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT I='
    + merged_bam + ' O=' + hist_log + ' H=' + hist_pdf + ' W=1000', shell=True)
#
refSeqTSS = peak_folder + 'mergeAll.RefSeqTSS'
subroutines.draw_TSS_insert(tssFile, merged_bam, refSeqTSS)
#
print('!!!!!!  call peak by macs2  !!!!!!')
peak_file = peak_folder + 'peaks'
subprocess.check_call('macs2 callpeak --nomodel -t ' + merged_bam + ' -n '
    + peak_file + ' --nolambda --keep-dup all --call-summits', shell=True)
print('!!!!!!  call peak done  !!!!!!')
#
summit = peak_folder + 'peaks_summits.bed'
filtered_peak = peak_folder + 'filtered_peaks.bed'
if options.blist:
    print('!!!!!!  filter peaks  !!!!!!')
    subprocess.check_call('bedtools intersect -v -a ' + summit + ' -b ' + options.blist
        + " | sort -k5 -nr > " + filtered_peak, shell=True)
    print('!!!!!!  filter peaks done  !!!!!!')
else:
    subprocess.check_call('sort -k5 -nr ' + summit + ' > ' + filtered_peak, shell=True)

print('!!!!!!  get top peaks by q-value  !!!!!!')
fold_rank = numpy.loadtxt(filtered_peak, 'str', delimiter='\t')
fold_rank[:, 1] = numpy.array(map(int, fold_rank[:, 1])) - 249  # 250
fold_rank[:, 2] = numpy.array(map(int, fold_rank[:, 2])) + 250
toppeaks = peak_folder + 'temp01.bed'
top_peaks = peak_folder + 'top_peaks.bed'
with open(toppeaks, 'w') as output:
    for peak in fold_rank:
        if float(peak[-1])>=float(options.logq):
#            print >> output, peak[0]+'\t'+peak[1]+'\t'+peak[2]
            output.write(peak[0]+'\t'+peak[1]+'\t'+peak[2]+'\n')
subprocess.check_call('bedtools sort -i ' + toppeaks + ' > ' + top_peaks, shell=True)
print('!!!!!!  get top peaks done  !!!!!!')
#
#
