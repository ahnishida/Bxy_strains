#!/bin/bash
cd results/processing/final_assemblies/
ls *.fasta > genome_list.txt
fastANI -t 16 --ql genome_list.txt --rl genome_list.txt -o fastani_res.txt
cd ../
mkdir fastani
mv final_assemblies/*.txt fastani/
