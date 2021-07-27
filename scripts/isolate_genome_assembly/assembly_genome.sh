#!/bin/bash

# script visualize read quality with fastqc
# quality filters reads using trimmmomatic
# assembly using spades isolate mode
# assess assembly quality with quast
# assess taxonomy of contigs with centrifuge
#(quast requires separate conda env to run py3.6)

isolate=$1

#lane1 reads
L1_folder=SA19175_${isolate}_L001
#lane2 reads
L2_folder=SA19175_${isolate}_L002

L1_forward_fastq=data/isolate_rawfastq/${L1_folder}/*R1_001.fastq.gz
L1_reverse_fastq=data/isolate_rawfastq/${L1_folder}/*R2_001.fastq.gz
L2_forward_fastq=data/isolate_rawfastq/${L2_folder}/*R1_001.fastq.gz
L2_reverse_fastq=data/isolate_rawfastq/${L2_folder}/*R2_001.fastq.gz

#combine reads
mkdir -pv results/processing/tmp_L1_L2_fastq/${isolate}
combined_forward_fastq=results/processing/tmp_L1_L2_fastq/${isolate}/${isolate}_R1.fastq.gz
combined_reverse_fastq=results/processing/tmp_L1_L2_fastq/${isolate}/${isolate}_R2.fastq.gz

cat $L1_forward_fastq $L2_forward_fastq > $combined_forward_fastq
cat $L1_reverse_fastq $L2_reverse_fastq > $combined_reverse_fastq

mkdir -pv results/processing/fastqc/${isolate}
fastqc $combined_forward_fastq  $combined_reverse_fastq -t 8  -o results/processing/fastqc/${isolate}

mkdir -pv results/processing/trimmomatic/${isolate}
trimmomatic PE $combined_forward_fastq $combined_reverse_fastq \
            results/processing/trimmomatic/${isolate}/${isolate}_R1_paired.fastq.gz results/processing/trimmomatic/${isolate}/${isolate}_R1_unpaired.fastq.gz \
            results/processing/trimmomatic/${isolate}/${isolate}_R2_paired.fastq.gz results/processing/trimmomatic/${isolate}/${isolate}_R2_unpaired.fastq.gz \
            LEADING:10 TRAILING:10 SLIDINGWINDOW:5:20 MINLEN:151 \
            -threads 12

mkdir -pv results/processing/assemblies/${isolate}
spades.py -1 results/processing/trimmomatic/${isolate}/${isolate}_R1_paired.fastq.gz -2 results/processing/trimmomatic/${isolate}/${isolate}_R2_paired.fastq.gz \
	-o results/processing/assemblies/${isolate}/spades_error_corrected_reads -t 50 -m 500 --only-error-correction
cor_err_forward_fastq=results/processing/assemblies/${isolate}/spades_error_corrected_reads/corrected/${isolate}_R1_paired.fastq.00.0_0.cor.fastq.gz
cor_err_reverse_fastq=results/processing/assemblies/${isolate}/spades_error_corrected_reads/corrected/${isolate}_R2_paired.fastq.00.0_0.cor.fastq.gz

#spades isolate mode
spades.py -1 $cor_err_forward_fastq -2 $cor_err_reverse_fastq \
          -o results/processing/assemblies/${isolate}/spades_isolate_assembly -t 8 --only-assembler --isolate

source activate quast-env
#quast and centrifuge require separate env to run py3.6 spades requires>py3.7
quast.py results/processing/assemblies/${isolate}/spades_isolate_assembly/scaffolds.fasta -o results/processing/assemblies/${isolate}/spades_isolate_assembly/quast/

centrifuge -f -x bin/centrifuge_taxonomy/p_compressed+h+v \
results/processing/assemblies/${isolate}/spades_isolate_assembly/contigs.fasta \
-S results/processing/assemblies/${isolate}/spades_isolate_assembly/centrifuge_hits.tsv
conda deactivate

rm -r results/processing/tmp_L1_L2_fastq
