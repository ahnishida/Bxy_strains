#!/bin/python

#adds assembly quality info from checkm,
#taxonomy from gtdbtk, and coverage info from bbmap.sh
#to isolate metadata
 
import pandas as pd
import json
import os

#read-in metadata file
os.chdir('data/ncbi_genomes')
metadata_file = 'ncbi_isolate_genomes_metadata.txt'
metadata = pd.read_csv(metadata_file,sep='\t')

#read-in gtdbtk results
gtdbtk_file = 'gtdbtk/gtdbtk.bac120.summary.tsv'
gtdbtk = pd.read_csv(gtdbtk_file,sep='\t')
gtdbtk['isolate'] = gtdbtk['user_genome'].apply(lambda x: x.split('_genomic')[0])
gtdbtk['taxonomy_Genus'] = gtdbtk['classification'].apply(lambda x: x.split(';')[-2].split('__')[1])
gtdbtk['taxonomy_Species'] = gtdbtk['classification'].apply(lambda x: x.split(';')[-1].split('__')[1])
gtdbtk = gtdbtk[['isolate','classification','taxonomy_Genus','taxonomy_Species']]
print(gtdbtk.head())
metadata = metadata.merge(gtdbtk,on='isolate',how='left')

#read-in checkm results
def read_checkm_tables(checkm_file):
    #parses 2 column in checkm results from dictionary to dataframe
    checkm = pd.read_csv(checkm_file,sep='\t',header=None)
    checkm.columns = ['isolate','stats']
    checkm['isolate']=checkm['isolate'].apply(lambda x: x.split('_genomic')[0])
    checkm['stats'] = checkm['stats'].apply(lambda x: json.loads(x.replace("'",'"'))) #format to dictionary 
    stats = pd.DataFrame(checkm['stats'].tolist())
    checkm = pd.concat([checkm['isolate'], stats], axis=1) 
    return(checkm)
    
checkm_file = 'checkm_prokka/storage/bin_stats_ext.tsv'
checkm = read_checkm_tables(checkm_file)
checkm['isolate'] =checkm['isolate'].apply(lambda x: x.split('_genomic')[0])
checkm = checkm[['isolate','Completeness', 'Contamination', 'GC', 'GC std', 
'Genome size', '# ambiguous bases', '# scaffolds', 
'# contigs', 'Longest scaffold', 'Longest contig', 
'N50 (scaffolds)', 'N50 (contigs)', 'Mean scaffold length', 
'Mean contig length', 'Coding density', 'Translation table', 
'# predicted genes']]
metadata = metadata.merge(checkm,on='isolate',how='left')
print(metadata.head())
metadata.to_csv('ncbi_genome_assembly_quality.txt',sep='\t',index=False)

#metadata = metadata[metadata['taxonomy_Genus'].isin(['Parabacteroides','Bacteroides','Bifidobacterium'])] 
#metadata = metadata[metadata['Completeness']>90]
#metadata.to_csv('metadata/strain_metadata_assembly_quality.txt',sep='\t',index=False)
