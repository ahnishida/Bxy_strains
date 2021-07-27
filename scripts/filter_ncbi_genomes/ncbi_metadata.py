#!/bin/python

#determine which ncbi genomes are from isolate, not MAGs,
#outputs these isolate genomes to separate folder
#formats simplified metadata

import pandas as pd
import re
import os

os.chdir('data/ncbi_genomes')
infile = 'prokaryotes_downloaded.txt'
outfile = 'prokaryotes_metadata.txt'

prok_dn = pd.read_csv(infile,sep='\t')
prok_dn = prok_dn[prok_dn['metadata_dn']==True] #subset to downloaded metadata
print(len(prok_dn))

print('Bioproject accessions w/ greatest number of genomes')
print(prok_dn['BioProject Accession'].value_counts()[:10])
#Top 3 bioproject are all isolates
#PRJNA544527: BIO-ML: The Broad Institute-OpenBiome Microbiome Library (>7,000 isolates)
#PRJNA482748: 1,520 reference genomes from cultivated human gut bacteria isolates
#PRJNA637878: 2359 bacterial isolates representing 1255 strains that were isolated

def get_metadata(genome,attribute):
	#returns metadata from ncbi biosample and bioproject data for a given attribute
	df = pd.read_csv(f'metadata/{genome}.txt',sep='\t',index_col=0)
	try:
		return(df.loc[attribute,genome])
	except:
		return('missing')

#fetch metadata infor for host,sample_type, isolation source and location
prok_dn['host'] = prok_dn['genome'].apply(lambda x: get_metadata(x,'host'))
prok_dn['sample_type'] = prok_dn['genome'].apply(lambda x: get_metadata(x,'sample_type'))
prok_dn['isolation_source'] = prok_dn['genome'].apply(lambda x: get_metadata(x,'isolation_source'))
prok_dn['geo_loc_name'] = prok_dn['genome'].apply(lambda x: get_metadata(x,'geo_loc_name'))

print('Problem: there are metagenome assembled genomes in these genomes')
print(prok_dn.loc[:,'sample_type'].value_counts())

#removes genomes labels as metagenomic assembly or lack strain name
prok_iso = prok_dn[prok_dn.loc[:,'sample_type'] != 'metagenomic assembly']
print(len(prok_iso),'ncbi genomes determined to be isolates')
prok_iso = prok_iso[prok_iso.loc[:,'Strain'] != '-']

print(len(prok_iso),'ncbi genomes determined to be isolates')
print('values for host and isolation source are messy, many samples have info missing from these attributes')
print(prok_iso.loc[:,'host'].value_counts())
#clean up host designations
prok_iso.loc[:,'host'] = prok_iso.loc[:,'host'].str.replace(r'ob/ob mouse|mouse', 'Mus musculus',
													flags=re.IGNORECASE,regex=True)
prok_iso.loc[:,'host'] = prok_iso.loc[:,'host'].str.replace(r'not applicable|ATCC strain', 'missing',
													flags=re.IGNORECASE,regex=True)

print('several samples missing host info')
print(prok_iso.loc[:,'host'].value_counts())
def find_missing_host(genome):
	#returns host species for missing samples by matching species names across all metadata attributes
	df = pd.read_csv(f'metadata/{genome}.txt',sep='\t',index_col=0)
	human = df.loc[:,genome].str.count(r'human|homo sapiens',flags=re.IGNORECASE).sum()
	chicken = df.loc[:,genome].str.count(r'chicken|gallus gallus', flags=re.IGNORECASE).sum()
	mouse = df.loc[:,genome].str.count(r'mouse|mus musculus', flags=re.IGNORECASE).sum()
	pig = df.loc[:,genome].str.count(r'pig|sus scrofa', flags=re.IGNORECASE).sum()
	total = pd.Series([human,chicken,mouse,pig],index=['Homo sapiens','Gallus gallus','Mus musculus','Sus scrofa'])
	hits = total[total>0]
	if len(hits)==0: #if no host species found
		return('missing')
	elif len(hits==1):
		return(hits.index[0]) #exactly one host species is found
	else:
		return('multiple') #check to see words matching multiple species are found


prok_iso.loc[:,'host_inferred'] = prok_iso.loc[:,'genome'].apply(find_missing_host)
print('checking method to infer host species is consistent')
print(prok_iso.groupby(['host','host_inferred']).count()) #check to see how host inferred compares to original host designation
print('resolved host for 4 chicken and 12 human samples, 18 still missing')
sp_dict = {'Homo sapiens':'human',
			'Gallus gallus':'chicken',
			'Mus musculus':'mouse',
			'Sus scrofa':'pig',
			'missing':'missing'}
prok_iso.loc[:,'host_inferred_name'] = prok_iso.loc[:,'host_inferred'].apply(lambda x: sp_dict[x]) #get common name


print('several samples missing isolation source info')
print(prok_iso['isolation_source'].value_counts())
def find_missing_isolation_source(genome):
	#infers isolation source for samples by matching key descriptive words
	df = pd.read_csv(f'metadata/{genome}.txt',sep='\t',index_col=0)
	stool = df.loc[:,genome].str.count(r'stool|faecal|feces|fecal',flags=re.IGNORECASE).sum()
	ceacum = df.loc[:,genome].str.count(r'ceacal|cecal|caecum',flags=re.IGNORECASE).sum()
	purulent = df.loc[:,genome].str.count(r'purulent',flags=re.IGNORECASE).sum()
	sewage = df.loc[:,genome].str.count('sewage',flags=re.IGNORECASE).sum()
	dietary_supplement = df.loc[:,genome].str.count(r'dietary supplement|probitic products',flags=re.IGNORECASE).sum()
	colon = df.loc[:,genome].str.count('colon',flags=re.IGNORECASE).sum()
	blood = df.loc[:,genome].str.count('blood',flags=re.IGNORECASE).sum()
	saliva = df.loc[:,genome].str.count('saliva',flags=re.IGNORECASE).sum()
	appendix = df.loc[:,genome].str.count('appendix',flags=re.IGNORECASE).sum()
	total = pd.Series([stool,ceacum,purulent,sewage,dietary_supplement,colon,blood,saliva,appendix],
				index=['stool','ceacum','purulent','sewage','dietary_supplement','colon','blood','saliva','appendix'])
	hits = total[total>0]
	if len(hits)==0:
		return('missing')
	elif len(hits==1):
		return(hits.index[0])
	else:
		return('multiple')

prok_iso.loc[:,'isolation_source_inferred'] = prok_iso.loc[:,'genome'].apply(find_missing_isolation_source)
#check to see how inferred isolation source compares to original isolation source designation
print(prok_iso.groupby(['isolation_source','isolation_source_inferred']).count())
print('reclassified isolation source for multiple species')

#get info on infection status
def infection_status(genome):
	#returns the number of instances infection is mentioned in the sample biodata
	df = pd.read_csv(f'metadata/{genome}.txt',sep='\t',index_col=0)
	infection = df.loc[:,genome].str.count('infection').sum()
	if infection > 0:
		return(True)
	else:
		return(False)
prok_iso.loc[:,'infection_status'] = prok_iso.loc[:,'genome'].apply(infection_status)

#subset to relevant columns to match isolate metadata
prok_iso = prok_iso[['genome','BioSample Accession','geo_loc_name',
'host_inferred_name','host_inferred','isolation_source_inferred','infection_status']]
prok_iso.columns = ['isolate','sample','site','host','genus_sp','isolation_source','infection_status']
prok_iso.to_csv(outfile,sep='\t',index=False)

print('copying',len(prok_iso),'isolates genomes to new folder')
##move over genomes determined to be isolates
os.system('rm -r isolate_genomes')
os.system('mkdir -p isolate_genomes')
for isolate in prok_iso['isolate']:
	os.system(f'cp genomes/{isolate}_genomic.fna.gz isolate_genomes/{isolate}.fna.gz')
