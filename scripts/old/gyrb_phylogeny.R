library(tidyverse)
library(ape)
library(ggtree)
library(phytools)
setwd('/stor/work/Ochman/alex/captive_ape_strain_genomes/')

gyrb_tree <- read.tree('results/gyrb/isolate_ASV_Bt_gyrb.fasta.aln.trimmed.tre')
gyrb_tree$tip.label <- str_replace(gyrb_tree$tip.label,'[-.]','_')

isol <- read.table('metadata/strain_ncbi_metadata_assembly_passing.txt',sep='\t') 
isol <- isol %>% 
  select(isolate,genus_sp,taxonomy_Genus,taxonomy_Species,dataset) %>%
  mutate(isolate = str_replace(isolate,'[-.]','_'))
colnames(isol) <- c('isolate','host','genus','genus_sp','dataset')

asv <- read.table('results/gyrb/gyrb_asv_hr_table.txt',sep='\t',header=TRUE)
unique(asv$Family)
f__F082 <- asv %>% 
  filter(Family=='f__F082')
MRCA <- ape::getMRCA(gyrb_tree,f__F082$ASV)
gyrb_tree <- root(gyrb_tree,node  = MRCA)

asv <- asv %>% select(ASV,HR_type,Genus) %>% 
  mutate(dataset = 'gyrb_amplicon')  %>%
  mutate(Genus = str_replace(Genus,'g__','')) 
colnames(asv) <- c('isolate','host','genus','dataset')
isol_asv <- isol %>% bind_rows(asv) %>% 
  filter(isolate %in% gyrb_tree$tip.label) %>% 
  mutate(host = recode(host,Pan_troglodytes='captive_chimp',
                       HR_human='human',
                       'Homo sapiens'='human',
                       Pongo_pygmaeus='captive_orangutan',
                       Pan_paniscus='captive_bonobo',
                       Gorilla_gorilla_gorilla='captive_gorilla',
                       HR_wild_gorilla='wild_gorilla',
                       HR_wild_bonobo='wild_bonobo',
                       HR_wild_chimp='wild_chimp',
                       'Mus musculus'='mouse',
                       'Gallus gallus'='chicken',
                       Unique_CP='mixed host',
                       MX_human_wild_apes='mixed host',
                       MX_wild_apes='mixed host'
                       ))
unique(isol_asv$host)

Bacteroides_table <- isol_asv %>% 
  filter(dataset == 'isolate_genomes' | dataset == 'ncbi_genomes') %>%
  filter(genus == "Bacteroides")
MRCA <- ape::getMRCA(gyrb_tree,Bacteroides_table$isolate)
Bacteroides_tree <- extract.clade(gyrb_tree,MRCA)

#add Bacteroides species to tree
unique(Bacteroides_table$genus_sp)
Bt_xylan = Bacteroides_table %>% filter(genus_sp == 'Bacteroides xylanisolvens')
Bt_xylan_node = findMRCA(Bacteroides_tree, as.vector(Bt_xylan$isolate))
Bt_ovatus = Bacteroides_table %>% filter(genus_sp == 'Bacteroides ovatus')
Bt_ovatus_node = findMRCA(Bacteroides_tree, as.vector(Bt_ovatus$isolate))
Bt_frag = Bacteroides_table %>% filter(genus_sp == 'Bacteroides fragilis')
Bt_frag_node = findMRCA(Bacteroides_tree, as.vector(Bt_frag$isolate))
Bt_fragA = Bacteroides_table %>% filter(genus_sp == 'Bacteroides fragilis_A')
Bt_fragA_node = findMRCA(Bacteroides_tree, as.vector(Bt_fragA$isolate))

#tree plot
ggtree(Bacteroides_tree)  %<+% isol_asv +
  geom_cladelabel(Bt_xylan_node,label='Bacteroides xylanisolvens',color='magenta') +
  geom_cladelabel(Bt_ovatus_node,label='Bacteroides ovatus',color='cyan') +
  geom_cladelabel(Bt_frag_node,label='Bacteroides fragilis',color='green2') +
  geom_cladelabel(Bt_fragA_node,label='Bacteroides fragilis_A',color='yellow4') +
  geom_tippoint(aes(color=host, shape=dataset, size=.25, alpha=.75)) +
  geom_treescale(x=.2,y=0)
ggsave('results/gyrb/Bacteroides_tree.pdf',height=6, width=6)

Parabacteroides_table <- isol_asv %>% 
  filter(dataset == 'isolate_genomes' | dataset == 'ncbi_genomes') %>%
  filter(genus == "Parabacteroides")
MRCA <- ape::getMRCA(gyrb_tree,Parabacteroides_table$isolate)
Parabacteroides_tree <- extract.clade(gyrb_tree,MRCA)

#add Parabacteroides species to tree
unique(Parabacteroides_table$genus_sp)
Pa_dist = Parabacteroides_table %>% filter(genus_sp == 'Parabacteroides distasonis')
Pa_dist_node = findMRCA(Parabacteroides_tree, as.vector(Pa_dist$isolate))
Pa_distA = Parabacteroides_table %>% filter(genus_sp == 'Parabacteroides distasonis_A')
Pa_distA_node = findMRCA(Parabacteroides_tree, as.vector(Pa_distA$isolate))

ggtree(Parabacteroides_tree) %<+% isol_asv + 
  geom_cladelabel(Pa_dist_node,label='Parabacteroides distasonis',color='magenta') +
  geom_cladelabel(Pa_distA_node,label='Parabacteroides distasonis_A',color='cyan') +
  geom_tippoint(aes(color=host, shape=dataset,size=.25, alpha=.75)) +
  geom_treescale(x=.2,y=0)
ggsave('results/gyrb/Parabacteroides_tree.pdf',height=6, width=6)


Bifidobacterium_tree <- read.tree('results/gyrb/isolate_ASV_Bif_gyrb.fasta.aln.trimmed.tre')
Bifidobacterium_tree$tip.label <- str_replace(Bifidobacterium_tree$tip.label,'[-.]','_')
plot(Bifidobacterium_tree)

bif_meta <- read.table('results/gyrb/Bifidobacteriaceae_metadata.txt',sep='\t',header=TRUE)
bif_meta <- bif_meta %>% select(ASV,host) %>% 
  mutate(dataset = 'gyrb_amplicon') 
colnames(bif_meta) <- c('isolate','host','dataset')

Bifidobacterium_table <- isol %>% bind_rows(bif_meta) %>% 
  filter(isolate %in% Bifidobacterium_tree$tip.label) %>% 
  mutate(host = recode(host,Pan_troglodytes='captive_chimp',
                       HR_human='human',
                       'Homo sapiens'='human',
                       Pongo_pygmaeus='captive_orangutan',
                       Pan_paniscus='captive_bonobo',
                       Gorilla_gorilla_gorilla='captive_gorilla',
                       HR_wild_gorilla='wild_gorilla',
                       HR_wild_bonobo='wild_bonobo',
                       HR_wild_chimp='wild_chimp',
                       'Mus musculus'='mouse',
                       'Gallus gallus'='chicken',
                       Unique_CP='mixed host',
                       MX_human_wild_apes='mixed host',
                       MX_wild_apes='mixed host'))
unique(Bifidobacterium_table$genus_sp)

Bmouk <- Bifidobacterium_table %>% 
  filter(genus_sp=='Bifidobacterium moukalabense')
MRCA <- ape::getMRCA(Bifidobacterium_tree,Bmouk$isolate)
Bifidobacterium_tree <- ape::root(Bifidobacterium_tree,node  = MRCA)

#add Bifidobacterium species to tree
unique(Bifidobacterium_table$genus_sp)
Bf_pseu = Bifidobacterium_table %>% filter(genus_sp == 'Bifidobacterium pseudolongum')
Bf_pseu_node = findMRCA(Bifidobacterium_tree, as.vector(Bf_pseu$isolate))
Bf_anim = Bifidobacterium_table %>% filter(genus_sp == 'Bifidobacterium animalis')
Bf_anim_node = findMRCA(Bifidobacterium_tree, as.vector(Bf_anim$isolate))
Bf_glob = Bifidobacterium_table %>% filter(genus_sp == 'Bifidobacterium globosum')
Bf_glob_node = findMRCA(Bifidobacterium_tree, as.vector(Bf_glob$isolate))
Bf_mouk = Bifidobacterium_table %>% filter(genus_sp == 'Bifidobacterium moukalabense')
Bf_mouk_node = findMRCA(Bifidobacterium_tree, as.vector(Bf_mouk$isolate))

#tree plot
ggtree(Bifidobacterium_tree)  %<+% Bifidobacterium_table +  xlim(NA, .5) +
  geom_cladelabel(Bf_pseu_node,label='Bif. pseudolongum',color='magenta') +
  geom_cladelabel(Bf_anim_node,label='Bif. animalis',color='cyan') +
  geom_cladelabel(Bf_glob_node,label='Bif. globosum',color='green2') +
  geom_cladelabel(Bf_mouk_node,label='Bif.moukalabense',color='yellow4') +
  geom_tippoint(aes(color=host, shape=dataset, size=.1, alpha=.75)) +
  geom_treescale(x=.2,y=0)
ggsave('results/gyrb/Bifidobacterium_tree.pdf',height=6, width=6)
