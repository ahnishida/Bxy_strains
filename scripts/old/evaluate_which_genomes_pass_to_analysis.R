#visualize genome quality based on completeness, contamination, #ofcontigs,
#determine filters for passing genomes to 
library(tidyverse)
library(estimatr)

setwd('/stor/work/Ochman/alex/captive_ape_strain_genomes/')
ncbi <- read.csv('data/ncbi_genomes/ncbi_genome_assembly_quality.txt',sep='\t')
ncbi <- ncbi %>% mutate(dataset = 'ncbi_genomes')
strain <- read.csv('metadata/strain_metadata_assembly_quality.txt',sep='\t')
strain <- strain %>% mutate(dataset = 'isolate_genomes',
                            isolation_source = 'stool',
                            infection_status = 'False')
metadata <- bind_rows(strain,ncbi)
write.table(metadata,'metadata/strain_ncbi_metadata_assembly_all.txt',sep='\t',row.names = F)

#several ncbi genomes have high contamination eliminate outliers
metadata <- metadata %>% filter(Contamination<5 & Completeness>95) %>%
  mutate(taxonomy_Species = str_replace_all(taxonomy_Species,' ','_')) %>%
  filter(taxonomy_Genus=='Bacteroides')
summary(lm_robust(Genome.size ~ Contamination + Completeness + host + taxonomy_Species, data=metadata))
summary(lm_robust(X..predicted.genes ~ Contamination + Completeness + host + taxonomy_Species, data=Bxy_metadata))

Bxy_metadata <- metadata %>% filter(taxonomy_Species== 'Bacteroides_xylanisolvens')
table(Bxy_metadata$host)
summary(lm_robust(Genome.size ~ Contamination + Completeness + host, data=Bxy_metadata))
summary(lm_robust(X..predicted.genes ~ Contamination + Completeness + host, data=Bxy_metadata))

summary <- metadata %>% group_by(taxonomy_Genus,taxonomy_Species,dataset) %>% 
  tally() %>% 
  pivot_wider(names_from = dataset,values_from = n,values_fill = 0)
write.table(metadata,'metadata/strain_ncbi_metadata_assembly_passing.txt',sep='\t',row.names = F,quote=F)



summary <- d %>% group_by(taxonomy_Genus,taxonomy_Species,dataset) %>% 
  tally() %>% 
  pivot_wider(names_from = dataset,values_from = n,values_fill = 0)
library(formattable)
system('mkdir -pv results/assembly_quality/')
genomes_summary <- formattable(summary,list(
  isolate_genomes = color_tile("transparent", "lightpink"),
  ncbi_genomes = color_tile("transparent", "lightblue")))
export_formattable(genomes_summary,file='results/assembly_quality/genomes_summary.pdf')
install.packages('htmltools')
install.packages('webshot')
library(htmltools)
library(webshot)
webshot::install_phantomjs()

para <- d %>%filter(taxonomy_Genus=='Parabacteroides')
lm.para <- lm_robust(X..predicted.genes ~ Completeness + Contamination + taxonomy_Species + N50..contigs.,data=para)
summary(lm.para)

bact <- d %>%filter(taxonomy_Genus=='Bacteroides')
lm.bact <- lm_robust(X..predicted.genes ~ Completeness + Contamination + taxonomy_Species + N50..contigs.,data=bact)
summary(lm.bact)

bif <- d %>%filter(taxonomy_Genus=='Bifidobacterium')
lm.bif <- lm_robust(X..predicted.genes ~ Completeness + Contamination + taxonomy_Species + N50..contigs.,data=bif)
summary(lm.bif)

new_order = c("Bifidobacterium","Bacteroides","Parabacteroides")
d$taxonomy_Genus <- factor(d$taxonomy_Genus, levels = new_order)

library(RColorBrewer)
darkcols <- brewer.pal(8, "Dark2")
set1 <- brewer.pal(3, "Set1")
sp_pal <- c(darkcols,set1)


(complete_v_contam <- d %>%
  ggplot() +
  aes(x=Completeness,y=Contamination,color=taxonomy_Species) + 
  geom_point(size=2,alpha=.5) +
  theme_bw() +
  facet_wrap(~taxonomy_Genus) +
  scale_color_manual(values=sp_pal) + 
  theme(legend.position="bottom"))
ggsave(complete_v_contam,file='results/assembly_quality/complete_v_contam.pdf',width=6,height=4)

(complete_v_genes <- d %>%
    ggplot() +
    aes(x=Completeness,y=X..predicted.genes,color=taxonomy_Species) + 
    geom_point(size=2,alpha=.5) +
    theme_bw() +
    facet_wrap(~taxonomy_Genus) +
    scale_color_manual(values=sp_pal)+ 
    theme(legend.position="bottom"))
ggsave(complete_v_genes,file='results/assembly_quality/complete_v_genes.pdf',width=6,height=4)

(contam_v_genes <- d %>%
    ggplot() +
    aes(x=Contamination,y=X..predicted.genes,color=taxonomy_Species) + 
    geom_point(size=2,alpha=.5) +
    theme_bw() +
    facet_wrap(~taxonomy_Genus) +
    scale_color_manual(values=sp_pal)+ 
    theme(legend.position="bottom"))
ggsave(contam_v_genes,file='results/assembly_quality/contam_v_genes.pdf',width=6,height=4)

(N50_v_genes <- d %>%
    ggplot() +
    aes(x=N50..contigs.,y=X..predicted.genes,color=taxonomy_Species) + 
    geom_point(size=2,alpha=.5) +
    theme_bw() +
    facet_wrap(~taxonomy_Genus) +
    scale_color_manual(values=sp_pal)+ 
    theme(legend.position="bottom"))
ggsave(N50_v_genes,file='results/assembly_quality/N50_v_genes.pdf',width=6,height=4)

(Xcontigs_v_genes <- d %>%
    ggplot() +
    aes(x=X..contigs,y=X..predicted.genes,color=taxonomy_Species) + 
    geom_point(size=2,alpha=.5) +
    theme_bw() +
    facet_wrap(~taxonomy_Genus) +
    scale_color_manual(values=sp_pal)+ 
    theme(legend.position="bottom"))
ggsave(Xcontigs_v_genes,file='results/assembly_quality/Xcontigs_v_genes.pdf',width=6,height=4)

