library(seqinr)
dir = 'results/pangenome/Bacteroides_xylanisolvens/roary'
aln <- read.alignment(file.path(dir,'core_gene_alignment.aln'),format='fasta')
alnd <- as.matrix(dist.alignment(aln, matrix = "similarity",gap=0))
write.table(alnd,file = file.path(dir,'core_gene_alignment_dist.txt'))

hc <- hclust(dist(alnd))

cluster <- cutree(hc,h=.01)
cluster <- as.data.frame(cluster) %>% rownames_to_column(var='isolate')
write.table(cluster,file = file.path(dir,'isolate_clusters_99.txt'),sep='\t')

cluster <- cutree(hc,h=.001)
cluster <- as.data.frame(cluster) %>% rownames_to_column(var='isolate')
write.table(cluster,file = file.path(dir,'isolate_clusters_999.txt'),sep='\t')