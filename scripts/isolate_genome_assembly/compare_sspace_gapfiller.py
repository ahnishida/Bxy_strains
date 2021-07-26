#compares number of contigs and genome length after running sspace and gapfiller
#some gapfiller fails to run on genomes where sspace does reduce the number of contigs
#copies either gapfiller contig file or sspace contig file to scaffold assemblies folder


from Bio import Seq
from Bio import SeqIO
import sys
import os
import pandas as pd

isolate = sys.argv[1]

sppace_contigs = "results/processing/sspace_gapfiller/"+isolate+"/"+isolate+".final.scaffolds.fasta" 
gapfilled_contigs = "results/processing/sspace_gapfiller/"+isolate+"/"+isolate+".gapfilled.final.fa"
os.system('mkdir -p results/processing/scaffold_assemblies') #outdir
if os.path.isfile(gapfilled_contigs):
	print('gapfilled fasta found')
	sppace_dict = SeqIO.to_dict(SeqIO.parse(sppace_contigs,'fasta'))
	gapfilled_dict = SeqIO.to_dict(SeqIO.parse(gapfilled_contigs,'fasta'))
	print('n contigs sspace:gapfilled',len(sppace_dict),":",len(gapfilled_dict))
	sspace_genome_size = pd.Series(sppace_dict.values()).apply(lambda x: len(x)).sum()
	gapfilled_genome_size = pd.Series(gapfilled_dict.values()).apply(lambda x: len(x)).sum()
	print('genome size ratio sspace:gapfilled',sspace_genome_size,gapfilled_genome_size)
	os.system('cp '+gapfilled_contigs+' results/processing/scaffold_assemblies/'+isolate+'.scaffolds.fasta')
elif os.path.isfile(sppace_contigs):
	print('gapfilled fasta not found, reporting results for sspace results only')
	sppace_dict = SeqIO.to_dict(SeqIO.parse(sppace_contigs,'fasta'))
	print('n contigs sspace',len(sppace_dict))
	sspace_genome_size = pd.Series(sppace_dict.values()).apply(lambda x: len(x)).sum()
	print('genome size sspace',sspace_genome_size)
	os.system('cp '+sppace_contigs+' results/processing/scaffold_assemblies/'+isolate+'.scaffolds.fasta')
else:
	print('gapfilled fasta not found, sspace not found')
		
