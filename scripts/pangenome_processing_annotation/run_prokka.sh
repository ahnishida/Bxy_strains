#!/bin/bash

source activate prokka
input="results/pangenome/prokka/todo.txt"
while IFS= read -r filepath
do
  echo $filepath
  isolate=${filepath##*/}
  isolate=${isolate%.*}
  echo $isolate
  prokka --force --cpus 50 --prefix $isolate --outdir results/pangenome/prokka/${isolate} $filepath
done < $input



