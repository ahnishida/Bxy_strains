#!/bin/python
import pandas as pd
import os
import sys

#runs prokka pangenome analysis for a single species
species = sys.argv[1]

metadata_file = 'metadata/strain_ncbi_metadata_assembly_passing.txt'
metadata = pd.read_csv(metadata_file,sep='\t')

#subset metadata to isolates belonging to species
species_df =  metadata[metadata['taxonomy_Species'] == species]
print(metadata['taxonomy_Species'].unique())
os.system(f'mkdir -pv results/pangenome/{species}/assemblies/')
print(species,len(species_df),'genomes')

#copy assemblies
def copy_assembly(dataset,isolate):
	if dataset == 'isolate_genomes':
		infile = f'results/processing/scaffold_assemblies/{isolate}.scaffolds.fasta'
	else:
		infile = f'data/ncbi_genomes/isolate_genomes/{isolate}_genomic.fna'	
	assembly = f'results/pangenome/{species}/assemblies/{isolate}.fna'
	os.system(f'cp {infile} {assembly}')
species_df.apply(lambda df: copy_assembly(df.dataset,df.isolate),axis=1)
			
#annotate with prokka because roary requires gff3 format, 
#checkm and gtdbtk prokka outputs missing gene fastas at the end of gff file
for isolate in species_df['isolate']:
	prokka_dir = f'results/pangenome/{species}/prokka/{isolate}/'
	assembly = f'results/pangenome/{species}/assemblies/{isolate}.fna'
	if not os.path.exists(f'{prokka_dir}'):
		os.system(f'prokka --force --prefix {isolate} --outdir {prokka_dir} {assembly}')

