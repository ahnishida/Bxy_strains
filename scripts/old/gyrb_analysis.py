#!/bin/python
import pandas as pd
from Bio import SeqIO 
from Bio.Seq import Seq
import os
from Bio import pairwise2

##Extract gyrb fasta from gtdbtk intermediate results 
#def extract_gyrb_hit(dataset,isolate):
#	if dataset == 'isolate_genomes':
#		tigrfam_file = f'results/processing/gtdbtk/identify/intermediate_results/marker_genes/{isolate}.scaffolds/{isolate}.scaffolds_tigrfam_tophit.tsv'
#		fasta_file = f'results/processing/gtdbtk//identify/intermediate_results/marker_genes/{isolate}.scaffolds/{isolate}.scaffolds_protein.fna'
#	else:
#		tigrfam_file = f'data/ncbi_genomes/gtdbtk/identify/intermediate_results/marker_genes/{isolate}_genomic/{isolate}_genomic_tigrfam_tophit.tsv'
#		fasta_file = f'data/ncbi_genomes/gtdbtk/identify/intermediate_results/marker_genes/{isolate}_genomic/{isolate}_genomic_protein.fna'
#		
#	print(isolate)
#	with open(tigrfam_file) as f:
#		gyrb_hit = ['no_gyrb_hit_found',1]
#		for l in f:
#			l=l.strip('\n')
#			if 'TIGR01059' in l: #gyrb hmm profile
#				gene = l.split('\t')[0] #get best hit fasta header
#				eval = float(l.split('\t')[1].split(',')[1])
#				if eval < gyrb_hit[1]:
#					gyrb_hit = [gene,eval]
#	if gyrb_hit[0] == 'no_gyrb_hit_found':	
#		return('no_gyrb_hit_found')
#	else:
#		gene = gyrb_hit[0]
#		os.system('mkdir -p results/gyrb/fastas')				
#		for rec in SeqIO.parse(fasta_file, "fasta"): 
#			if gene == str(rec.id):
#				rec.id = isolate 
#				rec.description = ''
#				SeqIO.write(rec, f'results/gyrb/fastas/{isolate}_gyrb.fasta', "fasta")	
#				return(gene)
##test isolates		
#extract_gyrb_hit('isolate_genomes','P17-A2')
#extract_gyrb_hit('ncbi_genomes','GCA_004793765.1_ASM479376v1')
#
#metadata_file = 'metadata/strain_ncbi_metadata_assembly_passing.txt'
#metadata = pd.read_csv(metadata_file,sep='\t')
#print(len(metadata),'isolate and ncbi genomes')
#metadata['gyrb_hit'] = metadata.apply(lambda df: extract_gyrb_hit(df.dataset,df.isolate),axis=1)
#metadata = metadata[metadata['gyrb_hit']!='no_gyrb_hit_found'] #subset to genomes gyrb seqs
#print(len(metadata),'isolate and ncbi genomes with gyrb hits')
#
##copy over gyrb_amplicon datasets
#os.system('cp data/gyrb_amplicon_datasets/* results/gyrb/')
#
####Generate Bacteroidales phylogeny from gyrb amplicon seq
###Generate fasta containing all isolate gyrb sequences for Bacteroides and Parabacteroides
#Bt = metadata[metadata['taxonomy_Genus'].isin(['Bacteroides','Parabacteroides'])]
#with open("results/gyrb/isolate_Bt_gyrb.fasta", "w") as Bt_fasta:
#	for isolate in Bt['isolate']:
#		rec = SeqIO.parse(f'results/gyrb/fastas/{isolate}_gyrb.fasta', "fasta")
#		SeqIO.write(rec, Bt_fasta, "fasta")
#
###Combine Bt isolate gyrb fasta, ASV gyrb fasta
#with open("results/gyrb/isolate_ASV_Bt_gyrb.fasta", "w") as Bt_fasta:
#	rec = SeqIO.parse("results/gyrb/isolate_Bt_gyrb.fasta", "fasta") 
#	SeqIO.write(rec, Bt_fasta, "fasta")
#	ASVs = SeqIO.parse("results/gyrb/physeq_Bacteroidales_asv.fasta", "fasta") 
#	for ASV in ASVs:	
#		ASV.seq = ASV.seq[1:] #trim first nucleo so seq is in frame 1
#		ASVnum = ASV.id.split('_')[1]
#		if int(ASVnum) < 8000:
#			SeqIO.write(ASV, Bt_fasta, "fasta")
#
#os.system('transeq -sequence results/gyrb/isolate_ASV_Bt_gyrb.fasta -outseq results/gyrb/isolate_ASV_Bt_gyrb.faa')
#os.system('mafft results/gyrb/isolate_ASV_Bt_gyrb.faa > results/gyrb/isolate_ASV_Bt_gyrb.faa.aln')
#os.system('tranalign -asequence results/gyrb/isolate_ASV_Bt_gyrb.fasta -bsequence results/gyrb/isolate_ASV_Bt_gyrb.faa.aln -outseq results/gyrb/isolate_ASV_Bt_gyrb.fasta.aln')
#gyrb_seqs = SeqIO.parse("results/gyrb/isolate_ASV_Bt_gyrb.fasta.aln", "fasta")  
#with open("results/gyrb/isolate_ASV_Bt_gyrb.fasta.aln.trimmed",'w') as f:
#	for rec in gyrb_seqs:
#		rec.seq = rec.seq[342:612] #adjust
#		SeqIO.write(rec, f, "fasta")
#os.system('fasttree results/gyrb/isolate_ASV_Bt_gyrb.fasta.aln.trimmed > results/gyrb/isolate_ASV_Bt_gyrb.fasta.aln.trimmed.tre')
#
####Generate Bifidobcaterium phylogeny from gyrb amplicon seq
###Generate fasta containing all isolate gyrb sequences for Bifidobacterium
#Bif = metadata[metadata['taxonomy_Genus']=='Bifidobacterium']
#with open("results/gyrb/isolate_Bif_gyrb.fasta", "w") as Bif_fasta:
#	for isolate in Bif['isolate']:
#		rec = SeqIO.parse(f'results/gyrb/fastas/{isolate}_gyrb.fasta', "fasta")
#		SeqIO.write(rec, Bif_fasta, "fasta")
		
#Combine Bif isolate gyrb fasta, ASV gyrb fasta	
Bif_host = []
biftrans = Seq('DDGRGIPVDEVPGEGVSGVETVMTKLHAGGKFGGGGYAVSGGLHGVGISVVNALSTRVDIEVRRQGFHWTQTYVDQKPTSRLIKGEPMGEEESTGTSVT')	
#bif seqs are not all in the same frame. To be able to translate all seqs in frame1, 
#this script determines the reading frame and trim seq if needed

with open("results/gyrb/isolate_ASV_Bif_gyrb.fasta", "w") as Bif_fasta:
	rec = SeqIO.parse("results/gyrb/isolate_Bif_gyrb.fasta", "fasta") 
	SeqIO.write(rec, Bif_fasta, "fasta")
	ASVs = SeqIO.parse("results/gyrb/Bifidobacteriaceae.fna", "fasta") 
	for ASV in ASVs:	
		Bif_host.append(ASV.id) 
		frame1,frame2,frame3=ASV.seq.translate(),ASV.seq[1:].translate(),ASV.seq[2:].translate()
		best_hit = ['frame']
		frame1_score = pairwise2.align.globalxx(frame1, biftrans,score_only=True)
		if frame1_score > 80: #use 80 best the true frame is always in the 90s and alt frames in 20s 30s
			ASV.seq = ASV.seq	
		frame2_score = pairwise2.align.globalxx(frame2, biftrans,score_only=True)
		if frame2_score > 80:
			ASV.seq=ASV.seq[1:]
		frame3_score = pairwise2.align.globalxx(frame3, biftrans,score_only=True)
		if frame3_score > 80:
			ASV.seq=ASV.seq[2:]
		#print(frame1_score,frame2_score,frame3_score)
		if '*' not in ASV.seq.translate():
			SeqIO.write(ASV, Bif_fasta, "fasta")

#os.system('transeq -trim True -sequence results/gyrb/isolate_ASV_Bif_gyrb.fasta -outseq results/gyrb/isolate_ASV_Bif_gyrb.faa')
#os.system('mafft results/gyrb/isolate_ASV_Bif_gyrb.faa > results/gyrb/isolate_ASV_Bif_gyrb.faa.aln')
#os.system('tranalign -asequence results/gyrb/isolate_ASV_Bif_gyrb.fasta -bsequence results/gyrb/isolate_ASV_Bif_gyrb.faa.aln -outseq results/gyrb/isolate_ASV_Bif_gyrb.fasta.aln')
#gyrb_seqs = SeqIO.parse("results/gyrb/isolate_ASV_Bif_gyrb.fasta.aln", "fasta")  
#with open("results/gyrb/isolate_ASV_Bif_gyrb.fasta.aln.trimmed",'w') as f:
#	for rec in gyrb_seqs:
#		rec.seq = rec.seq[339:636] #adjust
#		SeqIO.write(rec, f, "fasta")
#os.system('fasttree results/gyrb/isolate_ASV_Bif_gyrb.fasta.aln.trimmed > results/gyrb/isolate_ASV_Bif_gyrb.fasta.aln.trimmed.tre')

#generate metadata for moeller_bif seqs
def which_species(read_id):
	if 'Gorilla' in read_id:
		return('wild_gorilla')
	elif 'Chimp' in read_id:
		return('wild_chimp')
	elif 'Bonobo' in read_id:
		return('wild_bonobo')
	else:
		return('human')
read_id = pd.Series(Bif_host)		
host = pd.Series(Bif_host).apply(which_species)
df = pd.concat([read_id,host],axis=1)
df.columns = ['ASV','host']
df.to_csv('results/gyrb/Bifidobacteriaceae_metadata.txt',sep='\t',index=False)