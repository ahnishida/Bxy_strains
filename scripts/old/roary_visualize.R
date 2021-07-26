rm(list=ls())
library(tidyverse)
library(vegan)
library(ggforce)
library(ggtree)

setwd('/stor/work/Ochman/alex/captive_ape_strain_genomes/')
species <- 'Bacteroides_xylanisolvens'

metadata <- read.table('metadata/strain_ncbi_metadata_assembly_passing.txt',header=TRUE)
metadata <- metadata %>% 
  mutate(taxonomy_Species = str_replace_all(taxonomy_Species,' ','_')) %>%
  filter(taxonomy_Species == species) %>%
  mutate(isolate = str_replace_all(isolate,'[_-]','.'),
  site = recode(site,'USA: Cambridge'='USA','USA:Boston'='USA',
                'China: Shenzhen'='China','USA:Seattle'='USA',
                'not applicable'='siteUnknown','missing'='siteUnknown')) %>%
 unite(host_site, host, site, sep = "_", remove = FALSE)
unique(metadata$host_site)
#
tree_file <- 'results/pangenome/Bacteroides_xylanisolvens/phylogeny/RAxML_bestTree.Bacteroides_xylanisolvens.rooted'
tree <- read.tree(tree_file)
tree_plot <- ggtree(tree) %<+% metadata +
  geom_tippoint(aes(color=host_site), alpha=0.8) 


#add cluster annotations 
cluster99 <- read.table(paste0('results/pangenome/',species,'/roary/isolate_clusters_99.txt'))
cluster99 <- cluster99 %>% mutate(isolate = str_replace_all(isolate,'[_-]','.'))
metadata <- metadata %>% left_join(cluster99,by='isolate')

cluster40 <- metadata %>% filter(cluster==40)
cluster40Node <- getMRCA(tree,cluster40$isolate)
cluster41 <- metadata %>% filter(cluster==41)
cluster41Node <- getMRCA(tree,cluster41$isolate)
cluster42 <- metadata %>% filter(cluster==42)
cluster42Node <- getMRCA(tree,cluster42$isolate)

(tree_plot <- tree_plot  + 
  geom_cladelabel(node=cluster40Node, color='black', offset=.001,label="Mixed-host captive clade") + 
  geom_cladelabel(node=cluster41Node, color='black', offset=.001,label="Captive gorilla clade2") +
  geom_cladelabel(node=cluster42Node, color='black', offset=.001,label="Captive gorilla clade1") +
   theme(legend.position = "none") + 
  xlim(NA, .1))

#NMDS of pangenomes
pan <- read_csv(paste0('results/pangenome/',species,'/roary/gene_presence_absence.csv'))
colnames(pan) <- str_replace_all(colnames(pan),'[_-]','.')
#are there any species genomes in metadata but not in pangenome analysis
setdiff(metadata$isolate,colnames(pan)) 
#species genomes
genomes  <- intersect(metadata$isolate,colnames(pan))

#format roary pres absence table as binary matrix where rows are genomes
#genes are columns and values represent where a gene is present or absent
isblank <- function(x) (ifelse(is.na(x),0,1))
pan <- pan %>% 
  select(-c("Non.unique Gene name":"Avg group size nuc"))  %>% 
  mutate_at(vars(-Gene),isblank) %>% 
  column_to_rownames(var='Gene') %>% 
  as.matrix() %>% t()

pan_dist <- vegdist(pan, method="jaccard") 



get_accessory <- function(pan,min,max){
  #subsets pan genomes matrix to only genes present above min and below max thresholds
  accessory_genes <- colnames(pan)[colSums(pan)<max & colSums(pan)>min]
  pan <- pan[,accessory_genes]
}
dim(pan)
min = 10
max = length(genomes)
acc <- get_accessory(pan,min,max)
dim(acc)

acc_dist <- vegdist(acc, method="jaccard") 
ord_nMDS<- metaMDS(comm= acc_dist,
                   k = 2,
                   maxit = 1000, #have to increase tries for unifrac distances, otherwise they don't converge
                   trymax = 500,
                   wascores = FALSE,
                   autotransform = FALSE,
                   trace = 2,
                   noshare = FALSE,
                   parallel = 4)
ndim_stress = paste0('k = ',ord_nMDS$ndim,', stress = ',round(ord_nMDS$stress,2)) #extract stress and dimensions

ord_nMDS_df <- as.data.frame(scores(ord_nMDS))%>% 
  rownames_to_column(var='isolate') %>% 
  left_join(metadata,by='isolate') 
ord_nMDS_df <- ord_nMDS_df %>% mutate(cluster= paste0('cluster',cluster),
                                      cluster = recode(cluster,'cluster40'="Mixed-host captive clade",
                                                               'cluster41'="Captive gorilla clade2",
                                                               'cluster42'="Captive gorilla clade1",
                                                       ))

nmds <- ggplot(ord_nMDS_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size=3,alpha=.5, aes(colour = host_site)) +
  theme_bw() + 
  geom_mark_ellipse(aes(label=cluster, fill=cluster,
                        filter = cluster %in% c("Mixed-host captive clade","Captive gorilla clade2","Captive gorilla clade1")),
                        label.buffer = unit(5, 'mm'))
library(cowplot)                              
tree_nmds <- plot_grid(tree_plot,nmds)
ggsave(tree_nmds,width=15,height=9,file = 'results/pangenome/Bacteroides_xylanisolvens/fulltree_nmds.pdf')



fastani <- read.table(paste0('results/pangenome/',species,'/fastani/fastani_res.txt'))
colnames(fastani) <- c('genome1','genome2','ANI','v4','v5')
fastani <- fastani %>% 
  mutate(genome1 =  str_replace(genome1,'.fna','')) %>%
  mutate(genome2 =  str_replace(genome2,'.fna',''))
head(fastani)

pan_dist <- vegdist(pan, method="jaccard") 
get_dist <- function(genome1,genome2){
  result = tryCatch({
  as.matrix(pan_dist)[genome1,genome2]
  }, finally = {
  as.matrix(pan_dist)[genome2,genome1] 
  })  
  return(result)
}
get_dist('P21-4G','P21-4E')

GeneNum <- rowSums(pan)
Genomes <- names(GeneNum)
GeneNum <- as.data.frame(cbind(Genomes,GeneNum))

intersect <- function(genome1,genome2) {
  genome1 = 'P21-4G'
  genome2 = 'P19-10A'
  pan_pair <- pan[c(genome1,genome2),]
  intersect = genome1[colSums(pan_pair) == 2]
  print(rowSums(pan_pair))
  print(length(intersect))
}
fastani$jaccard <- mapply(get_dist,fastani$genome1, fastani$genome2)
fastani <- fastani %>% filter(genome1!=genome2) 

ggplot(fastani) + 
  aes(x=ANI,jaccard) +
  geom_point(color='purple',alpha=.3) + 
  theme_bw()

tree = read.tree(paste0('results/pangenome/',species,'/phylogeny/RAxML_bestTree.',species))
ggtree(tree) %<+% metadata +
  geom_tippoint(aes(color=site, size=.1, alpha=.75)) +
  geom_treescale(x=.2,y=0)

