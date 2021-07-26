#!/bin/bash
species=$1
mkdir -pv results/pangenome/$species/phylogeny
perl bin/fasta2relaxedPhylip.pl \
-f results/pangenome/$species/roary_nosplitparalogs/core_gene_alignment.aln \
-o results/pangenome/$species/phylogeny/core_gene_alignment.phy
raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s results/pangenome/$species/phylogeny/core_gene_alignment.phy -n $species
mv RAxML*$species results/pangenome/$species/phylogeny/