#!/bin/bash
python scripts/filter_ncbi_genomes/ncbi_genomes_download.py
python scripts/filter_ncbi_genomes/ncbi_metadata.py
checkm lineage_wf -t 8 -x .fna data/ncbi_genomes/isolate_genomes/ data/ncbi_genomes/checkm_prokka
gtdbtk classify_wf --genome_dir data/ncbi_genomes/isolate_genomes/ --out_dir data/ncbi_genomes/gtdbtk -x fna --cpus 8 --pplacer_cpus 8
python scripts/filter_ncbi_genomes/ncbi_assess_assembly_quality.py
