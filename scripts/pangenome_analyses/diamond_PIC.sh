for isolate in P21.6E \
P21.11A \
P14.E4 \
GCA.003458755.1.ASM345875v1 \
GCA.000273315.1.Bact.xyla.CL03T12C04.V1 \
GCA.009102805.1.ASM910280v1 \
GCA.009102105.1.ASM910210v1 \
GCA.900114865.1.IMG.taxon.2654588180.annotated.assembly \
GCA.015551805.1.ASM1555180v1 \
GCA.003474645.1.ASM347464v1 \
GCA.009102085.1.ASM910208v1 \
GCA.003468875.1.ASM346887v1 \
GCA.009101945.1.ASM910194v1 \
GCA.900107825.1.IMG.taxon.2623620516.annotated.assembly \
GCA.000210075.1.ASM21007v1 \
GCA.004167295.1.ASM416729v1 \
GCA.003464445.1.ASM346444v1
do
ls results/pangenome/Bacteroides_xylanisolvens/faa/$isolate.faa
./../diamond blastp \
  -q results/pangenome/Bacteroides_xylanisolvens/faa/$isolate.faa \
  -d ../NCBI_db/nr/nr \
  -o results/pangenome/Bacteroides_xylanisolvens/diamond/${isolate}_blast.txt \
  -c1 \
  --fast\
  --max-target-seqs 50 \
  -f 6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore salltitles 
done
