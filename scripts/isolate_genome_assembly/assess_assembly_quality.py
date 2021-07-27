#!/bin/python

#adds assembly quality info from checkm,
#taxonomy from gtdbtk, and coverage info from bbmap.sh
#to isolate metadata

import pandas as pd
import json

metadata_file = 'metadata/strain_metadata.txt'
gtdbtk_file = 'results/processing/gtdbtk/gtdbtk.bac120.summary.tsv'
checkm_file = 'results/processing/checkm_prokka/storage/bin_stats_ext.tsv'
outfile = 'metadata/strain_metadata_assembly_quality.txt'
#read-in metadata file
metadata = pd.read_csv(metadata_file,sep='\t')

#read-in gtdbtk results
gtdbtk = pd.read_csv(gtdbtk_file,sep='\t')
gtdbtk['isolate'] = gtdbtk['user_genome'].apply(lambda x: x.split('.scaf')[0])
gtdbtk['taxonomy_Genus'] = gtdbtk['classification'].apply(lambda x: x.split(';')[-2].split('__')[1])
gtdbtk['taxonomy_Species'] = gtdbtk['classification'].apply(lambda x: x.split(';')[-1].split('__')[1])
gtdbtk = gtdbtk[['isolate','classification','taxonomy_Genus','taxonomy_Species']]
metadata = metadata.merge(gtdbtk,on='isolate',how='left')

#read-in checkm results
def read_checkm_tables(checkm_file):
    #parses 2 column in checkm results from dictionary to dataframe
    checkm = pd.read_csv(checkm_file,sep='\t',header=None)
    checkm.columns = ['isolate','stats']
    checkm['isolate']=checkm['isolate'].apply(lambda x: x.split('_')[0])
    checkm['stats'] = checkm['stats'].apply(lambda x: json.loads(x.replace("'",'"'))) #format to dictionary
    stats = pd.DataFrame(checkm['stats'].tolist())
    checkm = pd.concat([checkm['isolate'], stats], axis=1)
    return(checkm)

checkm = read_checkm_tables(checkm_file)
checkm['isolate'] =checkm['isolate'].apply(lambda x: x.split('.scaf')[0])
checkm = checkm[['isolate','Completeness', 'Contamination', 'GC', 'GC std',
'Genome size', '# ambiguous bases', '# scaffolds',
'# contigs', 'Longest scaffold', 'Longest contig',
'N50 (scaffolds)', 'N50 (contigs)', 'Mean scaffold length',
'Mean contig length', 'Coding density', 'Translation table',
'# predicted genes']]
metadata = metadata.merge(checkm,on='isolate',how='left')

def get_coverage(isolate):
	covstats_file = 'results/processing/covstats/'+isolate+'_covstats_log.txt'
	res = [isolate]
	try:
		with open(covstats_file,'r') as f:
			for l in f:
				if l.startswith('Percent mapped'):
					percent_mapped=l.split('\t')[-1].strip('\n')
					res.append(percent_mapped)
				if l.startswith('Average coverage'):
					coverage=l.split('\t')[-1].strip('\n')
					res.append(coverage)
				if l.startswith('Standard deviation'):
					coverage_std=l.split('\t')[-1].strip('\n')
					res.append(coverage_std)
		covstats=pd.Series(res,index=['isolate','percent_mapped','coverage_mean','coverage_std'],dtype='object')
	except:
		covstats=pd.Series([isolate,'NaN','NaN','NaN'],index=['isolate','percent_mapped','coverage_mean','coverage_std'],dtype='object')
	return(covstats)
covstats = metadata['isolate'].apply(get_coverage)
metadata = metadata.merge(covstats,on='isolate',how='left')

metadata = metadata[metadata['taxonomy_Genus'].isin(['Parabacteroides','Bacteroides','Bifidobacterium'])]
metadata = metadata[metadata['Completeness']>90]
metadata.to_csv(outfile,sep='\t',index=False)
