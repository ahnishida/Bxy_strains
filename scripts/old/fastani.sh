#!/bin/bash
species=$1
cd results/pangenome/$species/assemblies/
ls *.fna > genome_list.txt 
fastANI -t 16 --ql genome_list.txt --rl genome_list.txt -o fastani_res.txt
cd ../
mkdir fastani
mv assemblies/*.txt fastani/
