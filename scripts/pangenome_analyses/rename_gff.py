import pandas as pd
import sys 
import os

isolate_old = sys.argv[1]
isolate=isolate_old.replace('-','.')


pres_abs = pd.read_csv("results/pangenome/Bacteroides_xylanisolvens/roary_nosplitparalogs/gene_presence_absence.csv")
pres_abs = pres_abs[['Gene',isolate]]
pres_abs = pres_abs[pres_abs[isolate].notnull()]

GeneID_to_OG = {}
for gene,iso in zip(pres_abs['Gene'],pres_abs[isolate]):
	if '\t' in iso:
		geneids = iso.split('\t')
		for geneid in geneids:
			GeneID_to_OG[geneid]=gene
	else:
		GeneID_to_OG[iso]=gene

df = pd.read_csv("results/pangenome/Bacteroides_xylanisolvens/output/roary_nosplitparalogs_annotation.txt",sep='\t')
print(df.head())
OG_to_Annotation = dict(zip(df.Gene,df.Annotation))
OG_to_Annotation_consensus = dict(zip(df.Gene,df.func_consensus))
OG_to_dbcan = dict(zip(df.Gene,df.DIAMOND))
OG_to_eggnog = dict(zip(df.Gene,df.best_OG_desc)) 

def rename_gff(infile,outfile):
	with open(infile) as f:
		with open(outfile,'w') as g:
			for l in f:
				if len(l.split('\t'))>1:
					l=l.split('\t')
					desc=l[8]
					desc = desc.split(';')	
					try:
						ID=desc[0].split('ID=')[1]	
						gene = GeneID_to_OG[ID]
						print(gene)
						Gene_name = 'gene='+gene
						Name = 'Name='+GeneID_to_OG[ID]
						Prokka = 'Prokka='+OG_to_Annotation[gene] 
						try:
							Annotation_consensus = 'Annotation_consensus='+OG_to_Annotation_consensus[gene] 
						except:
							Annotation_consensus = 'Annotation_consensus=no_hit'	
						try:
							Eggnog = 'Eggnog='+OG_to_eggnog[gene]  
						except:
							Eggnog = 'Eggnog=no_hit'
						print(Eggnog)
						try:
							dbcan =	'dbcan='+OG_to_dbcan[gene] 
						except:
							dbcan =	'dbcan=no_hit'
						desc = ';'.join(['ID='+ID,Gene_name,Name,Annotation_consensus,Prokka,Eggnog,dbcan])
						print(desc)
						l[8] = desc
						l='\t'.join(l)+'\n'
						g.write(l)
					except:
						l='\t'.join(l)
						g.write(l)
				else:
					g.write(l)
				
os.system('mkdir results/pangenome/Bacteroides_xylanisolvens/renamed_gff')
infile = f'results/pangenome/prokka/{isolate_old}/{isolate_old}.gff'
outfile = f"results/pangenome/Bacteroides_xylanisolvens/renamed_gff/{isolate_old}_renamed.gff"
rename_gff(infile,outfile)				