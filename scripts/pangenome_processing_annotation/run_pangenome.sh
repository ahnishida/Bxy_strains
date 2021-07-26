#!/bin/bash
roary -p 50 -i 90 -s -f results/pangenome/Bacteroides_xylanisolvens/roary_nosplitparalogs/ -e -n -v results/pangenome/Bacteroides_xylanisolvens/gff/*.gff
./scripts/pangenome_processing_annotation/raxml_phylogeny.sh Bacteroides_xylanisolvens
mkdir -pv results/pangenome/Bacteroides_xylanisolvens/eggnog_mapper
emapper.py -m diamond --tax_scope 2 --itype CDS -i results/pangenome/Bacteroides_xylanisolvens/roary_nosplitparalogs/pan_genome_reference.fa \
--output_dir results/pangenome/Bacteroides_xylanisolvens/eggnog_mapper \
--override --cpu 50 \
-o pan_genome_reference

roary -p 50 -i 90 -s -f results/pangenome/Bacteroides_fragilis/roary_nosplitparalogs/ -e -n -v results/pangenome/Bacteroides_fragilis/gff/*.gff
./scripts/pangenome_processing_annotation/raxml_phylogeny.sh Bacteroides_fragilis
mkdir -pv results/pangenome/Bacteroides_fragilis/eggnog_mapper
emapper.py -m diamond --tax_scope 2 --itype CDS -i results/pangenome/Bacteroides_fragilis/roary_nosplitparalogs/pan_genome_reference.fa \
--output_dir results/pangenome/Bacteroides_fragilis/eggnog_mapper \
--override --cpu 50 \
-o pan_genome_reference

roary -p 50 -i 90 -s -f results/pangenome/Bacteroides_ovatus/roary_nosplitparalogs/ -e -n -v results/pangenome/Bacteroides_ovatus/gff/*.gff
./scripts/pangenome_processing_annotation/raxml_phylogeny.sh Bacteroides_ovatus
mkdir -pv results/pangenome/Bacteroides_ovatus/eggnog_mapper
emapper.py -m diamond --tax_scope 2 --itype CDS -i results/pangenome/Bacteroides_ovatus/roary_nosplitparalogs/pan_genome_reference.fa \
--output_dir results/pangenome/Bacteroides_ovatus/eggnog_mapper \
--override --cpu 50 \
-o pan_genome_reference

roary -p 50 -i 90 -s -f results/pangenome/Bacteroides_thetaiotaomicron/roary_nosplitparalogs/ -e -n -v results/pangenome/Bacteroides_thetaiotaomicron/gff/*.gff
./scripts/pangenome_processing_annotation/raxml_phylogeny.sh Bacteroides_thetaiotaomicron
mkdir -pv results/pangenome/Bacteroides_thetaiotaomicron/eggnog_mapper
emapper.py -m diamond --tax_scope 2 --itype CDS -i results/pangenome/Bacteroides_thetaiotaomicron/roary_nosplitparalogs/pan_genome_reference.fa \
--output_dir results/pangenome/Bacteroides_thetaiotaomicron/eggnog_mapper \
--override --cpu 50 \
-o pan_genome_reference
