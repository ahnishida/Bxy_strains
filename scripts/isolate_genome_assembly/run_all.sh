#!/bin/bash

export PATH=$PATH:bin/sspace_basic
export PATH=$PATH:bin/GapFiller_v1-11_linux-x86_64

#run all processing scripts to generate final assembly from raw fastq

source activate spades-env
input="metadata/isolates_names_head.txt"
while IFS= read -r isolate
do
  echo $isolate
  ./scripts/isolate_genome_assembly/assembly_genome.sh $isolate
  ipython scripts/isolate_genome_assembly/filter_contigs.ipynb $isolate
  ./scripts/isolate_genome_assembly/sspace_gapfiller.sh $isolate
  python scripts/isolate_genome_assembly/compare_sspace_gapfiller.py $isolate
  ./scripts/isolate_genome_assembly/bbmap_covstats.sh $isolate
done < $input

checkm lineage_wf -t 8 -x fasta results/processing/scaffold_assemblies/ results/processing/checkm_prokka/
gtdbtk classify_wf --genome_dir results/processing/scaffold_assemblies/ --out_dir results/processing/gtdbtk -x fasta --cpus 8 --pplacer_cpus 8
python scripts/isolate_genome_assembly/assess_assembly_quality.py
