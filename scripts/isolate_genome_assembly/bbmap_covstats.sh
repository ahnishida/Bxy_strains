#!/bin/bash

#script assesses genome coverage using bbmap from bbtools

isolate=$1

mkdir -p results/processing/covstats
Ffastq=results/processing/trimmomatic/${isolate}/${isolate}_R1_paired.fastq.gz
gapfiller_fasta=results/processing/scaffold_assemblies/${isolate}.scaffolds.fasta
covstats_output=results/processing/covstats/${isolate}_covstats.txt
covstats_log=results/processing/covstats/${isolate}_covstats_log.txt

bbmap.sh in=$Ffastq ref=$gapfiller_fasta covstats=$covstats_output 2> $covstats_log
