#!/bin/bash

#runs sspacer and gapfiller on filter_contigs fasta to close gaps
isolate=$1
echo $isolate

Ffastqgz=results/processing/trimmomatic/${isolate}/${isolate}_R1_paired.fastq.gz
Ffastq=results/processing/trimmomatic/${isolate}/${isolate}_R1_paired.fastq
Rfastqgz=results/processing/trimmomatic/${isolate}/${isolate}_R2_paired.fastq.gz
Rfastq=results/processing/trimmomatic/${isolate}/${isolate}_R2_paired.fastq
gunzip $Ffastqgz
gunzip $Rfastqgz

echo $isolate $Ffastq $Rfastq 350 0.75 FR>libraries_sspace.txt
echo $isolate bowtie $Ffastq $Rfastq 350 0.15 FR>libraries_gapfiller.txt
SSPACE_Basic.pl -l libraries_sspace.txt -s results/processing/filtered_contigs/${isolate}_contigs.fasta -b $isolate
echo GapFiller.pl -l libraries_gapfiller.txt -s ${isolate}.final.scaffolds.fasta -b $isolate
GapFiller.pl -l libraries_gapfiller.txt -s ${isolate}.final.scaffolds.fasta -b $isolate
mkdir -pv results/processing/sspace_gapfiller/${isolate}
mv ${isolate}.final.scaffolds.fasta results/processing/sspace_gapfiller/${isolate}/
mv ${isolate}.summaryfile.txt results/processing/sspace_gapfiller/${isolate}/
mv ${isolate}/${isolate}.gapfilled.final.fa results/processing/sspace_gapfiller/${isolate}/
mv ${isolate}/${isolate}.summaryfile.final.txt results/processing/sspace_gapfiller/${isolate}/

gzip $Ffastq
gzip $Rfastq

#clean up
rm -rf ${isolate}
rm -f ${isolate}*
rm -rf bowtieoutput
rm -rf intermediate_results  
rm -rf pairinfo
rm -rf reads
rm -rf ref
rm libraries*