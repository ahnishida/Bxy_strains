library(ape)
library(treeio)
tree_file <- "results/pangenome/Bacteroides_xylanisolvens_outgroup/phylogeny/RAxML_bipartitions.Bacteroides_xylanisolvens_outgroup"
tree <- ape::read.tree(tree_file)
tree$tip.label <- str_replace_all(tree$tip.label,'[_-]','.')
ggtree(tree) + geom_nodelab()

metadata <- read.table('metadata/strain_ncbi_metadata_assembly_passing.txt',header=TRUE)
metadata <- metadata %>% filter(as.numeric(Completeness)>90,as.numeric(Contamination)<5) %>%
                                mutate(taxonomy_Species = str_replace_all(taxonomy_Species,' ','_'),
                                isolate = str_replace_all(isolate,'[_-]','.'),
                                site = recode(site,'USA: Cambridge'='USA','USA:Boston'='USA',
                                'China: Shenzhen'='China','USA:Seattle'='USA',
                                 'not applicable'='siteUnknown','missing'='siteUnknown')) %>%
                                unite(host_site, host, site, sep = "_", remove = FALSE) 
Bov <- metadata %>% filter(taxonomy_Species == 'Bacteroides_ovatus') 
Bov_outgroup <- intersect(tree$tip.label,Bov$isolate)
tree <- root(tree,outgroup=Bov_outgroup) 

tree_plot <- ggtree(tree) %<+% metadata +
  geom_nodepoint(aes(subset = suppressWarnings(as.numeric(label)) > 50),size=.75) +
  geom_tippoint(aes(color=host_site), alpha=0.8)  +
  geom_tiplab() +
  xlim(NA,.5)
write.tree(tree,file=paste0(tree_file,'.rooted'))
ggsave(tree_plot,file='results/pangenome/Bacteroides_xylanisolvens_outgroup/phylogeny.pdf',height=15,width=8)

#check metadata for 3 ancestors, not in the 3 clades
metadata %>% 
  filter(isolate %in% c('GCA.008710235.1.ASM871023v1','GCA.006546965.1.ASM654696v1','GCA.015556195.1.ASM1555619v1'))

#assign superclades A,B,C
#cladeA mostly humans
cladeA_MCRA <- ape::getMRCA(tree,c('GCA.009102525.1.ASM910252v1','GCA.009093695.1.ASM909369v1'))
cladeA_tree <- extract.clade(tree,cladeA_MCRA)
ggtree(cladeA_tree) %<+% metadata +
  geom_nodepoint(aes(subset = suppressWarnings(as.numeric(label)) > 50),size=.75) +
  geom_tippoint(aes(color=host_site), alpha=0.8)  +
  geom_tiplab() +
  xlim(NA,.5)
cladeB_MCRA <- ape::getMRCA(tree,c('P17.D4','GCA.001405055.1.13414.6.23'))
cladeB_tree <- extract.clade(tree,cladeB_MCRA)
cladeB_tree_plot <- ggtree(cladeB_tree) %<+% metadata +
  geom_nodepoint(aes(subset = suppressWarnings(as.numeric(label)) > 50),size=.75) +
  geom_tippoint(aes(color=host_site), alpha=0.8)  +
  geom_tiplab() +
  xlim(NA,.5)
cladeC_MCRA <- ape::getMRCA(tree,c('P21.4E','GCA.009102685.1.ASM910268v1'))
cladeC_tree <- extract.clade(tree,cladeC_MCRA)
ggtree(cladeC_tree) %<+% metadata +
  geom_nodepoint(aes(subset = suppressWarnings(as.numeric(label)) > 50),size=.75) +
  geom_tippoint(aes(color=host_site), alpha=0.8)  +
  geom_tiplab() +
  xlim(NA,.5)


library(estimatr)
metadata_Bxy <- metadata %>% filter(taxonomy_Species == 'Bacteroides_xylanisolvens') 
summary(lm_robust(Genome.size ~ Contamination + Completeness + X..contigs, data= metadata_Bxy))
summary(lm_robust(X..predicted.genes ~ Contamination + Completeness + X..contigs, data= metadata_Bxy))


ggplot(metadata_Bxy,aes(x=Completeness,y=Genome.size)) + geom_point()
ggplot(metadata_Bxy,aes(x=Contamination,y=Genome.size)) + geom_point()
ggplot(metadata_Bxy,aes(x=Completeness,y=X..predicted.genes)) + geom_point()
ggplot(metadata_Bxy,aes(x=Contamination,y=X..predicted.genes)) + geom_point()
ggplot(metadata_Bxy,aes(x=X..contigs,y=Genome.size)) + geom_point()
ggplot(metadata_Bxy,aes(x=X..contigs,y=X..predicted.genes)) + geom_point()


metadata <- metadata %>% mutate(clade = ifelse(isolate %in% cladeA_tree$tip.label,'cladeA',
                                   ifelse(isolate %in% cladeB_tree$tip.label,'cladeB',
                                          ifelse(isolate %in% cladeC_tree$tip.label,'cladeC','unassigned'))))
metadata_clade <- metadata %>% filter(clade !='unassigned') 
lm_robust(X..predicted.genes ~ clade + X..contigs,data = metadata_clade)
kruskal.test(X..predicted.genes ~ clade, data = metadata_clade)
library(PMCMRplus)
posthoc.kruskal.dunn.test(X..predicted.genes ~ clade, data = metadata_clade)
                           
metadata_clade %>% kruskal_test(as.factor(clade) ~ as.numeric(X..predicted.genes))
metadata %>% 
  filter(clade !='unassigned') %>% 
  unite(clade_host_site, clade, host_site, sep = "_", remove = FALSE)  %>%
  ggplot(aes(x=clade,y=Genome.size,fill = host_site)) + 
  geom_dotplot(binaxis = "y",
               stackdir = "center",stackratio=.75, dotsize=.5,
                 aes(color = host_site)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) 
metadata %>% 
  filter(clade !='unassigned') %>% 
  unite(clade_host_site, clade, host_site, sep = "_", remove = FALSE)  %>%
  ggplot(aes(x=clade,y=X..predicted.genes)) + 
  geom_dotplot(binaxis = "y",
               stackdir = "center",stackratio=.75, dotsize=.5,
               aes(fill = host_site)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) 

cladeA_tree_plot
cladeA_genes <- metadata %>% 
  filter(isolate %in% cladeA_tree$tip.label) %>%
  select(isolate,X..predicted.genes) %>%
  column_to_rownames("isolate")
cladeA_tree_heat <- ggtree(cladeA_tree) %<+% metadata + geom_tippoint(aes(color=host_site), alpha=0.8)
gheatmap(cladeA_tree_heat, cladeA_genes) 

cladeB_tree_plot
cladeB_genes <- metadata %>% 
  filter(isolate %in% cladeB_tree$tip.label) %>%
  select(isolate,X..predicted.genes) %>%
  column_to_rownames("isolate")
cladeB_tree_heat <- ggtree(cladeB_tree) %<+% metadata + geom_tippoint(aes(color=host_site), alpha=0.8)
gheatmap(cladeB_tree_heat, cladeB_genes) 

cladeC_genes <- metadata %>% 
  filter(isolate %in% cladeC_tree$tip.label) %>%
  select(isolate,X..predicted.genes) %>%
  column_to_rownames("isolate")
cladeC_tree_heat <- ggtree(cladeC_tree) %<+% metadata + geom_tippoint(aes(color=host_site), alpha=0.8)
gheatmap(cladeC_tree_heat, cladeC_genes) 


