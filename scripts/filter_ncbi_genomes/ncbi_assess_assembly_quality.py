#!/bin/python

#adds assembly quality info from checkm,
#taxonomy from gtdbtk, and coverage info from bbmap.sh
#to isolate metadata

import pandas as pd
import json
import os

#inputs/outputs
os.chdir('data/ncbi_genomes')
metadata_file = 'prokaryotes_metadata.txt'
checkm_file = 'checkm_prokka/storage/bin_stats_ext.tsv'
gtdbtk_file = 'gtdbtk/gtdbtk.bac120.summary.tsv'
outfile = 'ncbi_genome_assembly_quality.txt'

#read-in metadata file
metadata = pd.read_csv(metadata_file,sep='\t')

#read-in gtdbtk results
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
metadata.to_csv(outfile,sep='\t',index=False)
