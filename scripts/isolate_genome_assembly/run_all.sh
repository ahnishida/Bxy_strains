#!/bin/bash

#export PATH=$PATH:bin/sspace_basic
#export PATH=$PATH:bin/GapFiller_v1-11_linux-x86_64

#run all processing scripts to generate final assembly from raw fastq

source activate spades-env
input="metadata/isolates_names.txt"
while IFS= read -r isolate
do
  echo $isolate
  #./scripts/processing_isolates/assembly_genome.sh $isolate
  #ipython scripts/processing_isolates/filter_contigs.ipynb $isolate
  #./scripts/processing_isolates/sspace_gapfiller.sh $isolate
  python scripts/processing_isolates/compare_sspace_gapfiller.py $isolate
  #./scripts/processing_isolates/bbmap_covstats.sh $isolate
done < $input

#checkm lineage_wf -t 8 -x fasta results/processing/scaffold_assemblies/ results/processing/checkm_prokka/
#gtdbtk classify_wf --genome_dir results/processing/scaffold_assemblies/ --out_dir results/processing/gtdbtk -x fasta --cpus 8 --pplacer_cpus 8 
#python scripts/processing_isolates/assess_assembly_quality.py


