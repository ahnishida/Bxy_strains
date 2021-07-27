import pandas as pd
import requests
import xml.etree.ElementTree as ET
from lxml import etree
import os
import wget

#create outdir for ncbi genomes
os.system('mkdir -p data/ncbi_genomes')
os.chdir('data/ncbi_genomes')

#download master list of all prokaryote genomes
#wget.download('https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt')
prok = pd.read_csv('prokaryotes.txt',sep='\t',dtype=object)

def genus_sp(orgName):
    try:
        genus_sp = orgName.split(' ')[0] + '_' + orgName.split(' ')[1]
        return(genus_sp)
    except:
        return(orgName.split(' ')[0])

def get_strain(orgName):
    try:
        strain = ' '.join(orgName.split(' ')[2:])
        return(strain)
    except:
        return('')

prok['Genus_sp'] = prok['#Organism/Name'].apply(lambda x: genus_sp(x))
prok['strain'] = prok['#Organism/Name'].apply(lambda x: get_strain(x))

isolates_sp = ['Bacteroides_xylanisolvens','Bacteroides_fragilis','Bacteroides_ovatus','Bacteroides_thetaiotaomicron']
prok = prok[prok['Genus_sp'].isin(isolates_sp)] #subset to only bacterial sp. that match isolate genomes sequences
print(prok['Genus_sp'].value_counts(),'genomes in ncbi database') #get idea of how many genomes to download

prok['fna'] = prok['FTP Path'].apply(lambda x : x.split('/')[-1]+'_genomic.fna.gz')
prok['genome']=prok['FTP Path'].apply(lambda x: x.split('/')[-1])

os.system('mkdir -p genomes')
prok_downloaded = os.listdir('genomes')
prok_to_dn = prok[~prok['fna'].isin(prok_downloaded)]
print('genomes to download',prok_to_dn['genome'])

#entrez
BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

def get_attributes_from_biosample(biosample_id,genome):
	#retrieves biosample attributes for a given prokaryote ncbi genome
    url = BASE_URL + "efetch.fcgi"
    parameters = {
        "db": "biosample",
        "id": biosample_id
    }
    response = requests.get(url, params=parameters)
    try:
        parser = etree.XMLParser(recover=True)
        root = ET.fromstring(response.text,parser=parser)
    except AttributeError:
        raise
    attribute_dict = {}
    biosample = root.findall("./BioSample")[0]
    attributes = biosample.findall("./Attributes")[0]
    colnames = []
    values = []
    for attribute in attributes:
    	colnames.append(attribute.attrib['attribute_name'])
    	values.append(attribute.text)
    res = pd.Series(values,index = colnames)
    print(res)
    return(res)

def get_attributes_from_bioproject(bioproject_id,genome):
	#retrieves bioproject title and description a given prokaryote ncbi genome
    url = BASE_URL + "efetch.fcgi"
    parameters = {
        "db": "bioproject",
        "id": bioproject_id
    }
    response = requests.get(url, params=parameters)
    try:
    	parser = etree.XMLParser(recover=True)
    	root = ET.fromstring(response.text,parser=parser)

    except AttributeError:
        raise
    ProjectDescr = root[0][0][1]
    Title = ProjectDescr.find("./Title")
    Description = ProjectDescr.find("./Description")
    res = pd.Series([Title.text,Description.text],index = [Title.tag,Description.tag])
    print(res)
    return(res)

def fetch_metadata(genome,bioproject,biosample):
	#combine and output bioproject and biosample data to a table
	bioproject = get_attributes_from_bioproject(bioproject,genome)
	biosample = get_attributes_from_biosample(biosample,genome)
	biodata = bioproject.append(biosample)
	biodata = biodata.to_frame(name=genome)
	biodata.to_csv('metadata/'+genome+'.txt',sep='\t')

for ftp,bioproject_id,biosample_acc in zip(prok_to_dn['FTP Path'],prok_to_dn['BioProject ID'],prok_to_dn['BioSample Accession']):
	genome=ftp.split('/')[-1]
	fna=ftp+'/'+genome+'_genomic.fna.gz'
	wget.download(fna,out='genomes')
	try:
		fetch_metadata(genome,bioproject_id,biosample_acc)
	except:
		print('metadata download failed')
		continue

metadata_dn = [x.split('.txt')[0] for x in os.listdir('metadata')]
prok['metadata_dn'] = prok['genome'].isin(metadata_dn)
print('metadata_dn')
print(prok['metadata_dn'].value_counts())
prok.to_csv('prokaryotes_downloaded.txt',sep='\t',index=None)
