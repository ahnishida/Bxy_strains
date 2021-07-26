#rm(list=ls())
library(ape)
library(ggtree)
library(treeio)
library(stringr)
library(tidyverse)
library(tidytable)

setwd('/stor/work/Ochman/alex/captive_ape_strain_genomes/')
 
metadata <- read.table('metadata/strain_ncbi_metadata_assembly_all.txt',header=T)
metadata <- metadata %>%
  mutate(taxonomy_Species = str_replace_all(taxonomy_Species,' ','_'),
         isolate = str_replace_all(isolate,'[_-]','.'),
         site = recode(site,'USA: Cambridge'='USA','USA:Boston'='USA',
                       'China: Shenzhen'='China','USA:Seattle'='USA',     
                       'USA:Baltimore'='USA', 'not applicable'='siteUnknown',                          'missing'='siteUnknown')) %>%
  unite(host_site, host, site, sep = "_", remove = FALSE) 

indir = file.path('results/analyses_input_files')
tree_file <- file.path(indir,'Bacteroides_xylanisolvens.tre')
table_file <- file.path(indir,'roary_nosplitparalogs_gene_presence_absence.txt')
outdir <- 'results/pangenome/Bacteroides_xylanisolvens/count'


tree <- read.tree(tree_file)
tree <- makeNodeLabel(tree)
plot(tree)
(tree$tip.label <- stringr::str_replace_all(tree$tip.label,'[_-]','.')) #count doesn't like these
write.tree(tree,file=file.path(outdir,'Bxy.tree'))
           
#format files for count
table <- read_tsv(table_file)
colnames(table) <- stringr::str_replace_all(colnames(table),'[_-]','.')
isblank <-  function(x) (ifelse(is.na(x),0,1))
table <- table %>% 
  select(Gene,all_of(tree$tip.label))  %>% 
  mutate_at(vars(-Gene),isblank)
setdiff(tree$tip.label,colnames(table)) #make sure tree taxa in table
write.table(table,file=file.path(outdir,'pres_abs.txt'),
            sep='\t', quote=FALSE,row.names = FALSE)


runCount <- function(tree,data,output,gain) {
system(paste('java -Xmx2048M -cp bin/Count/Count.jar ca.umontreal.iro.evolution.genecontent.AsymmetricWagner -gain ',
             2,
             tree,
             data,'>',
             output,sep=' '))
system(paste0("grep '# CHANGE' ",output," | sed 's/# //' > ",output,".CHANGE"))
system(paste0("grep '# PRESENT' ",output," | sed 's/# //' > ",output,".PRESENT"))
system(paste0("grep '# FAMILY' ",output," | sed 's/# //' > ",output,".FAMILY"))
present <- read.table(paste0(output,".PRESENT"),sep='\t',header=TRUE)
res <- present %>% 
  mutate(istip = if_else(node %in% metadata$isolate,'tip','node')) %>% 
  rstatix::kruskal_test(genes ~ istip)
res$gain_penalty <- gain
return(res)
}

#compare gain penalty 
outdir <- 'results/pangenome/Bacteroides_xylanisolvens/count'
tree <- file.path(outdir,'Bxy.tree')
data <- file.path(outdir,'pres_abs.txt')
output <- file.path(outdir,'Bxy_count')
runCount(tree,data,output,2)

#gain1 <- runCount(tree,data,output,1)
#gain1.25 <- runCount(tree,data,output,1.25)
gain1.5 <-runCount(tree,data,output,1.5)
gain1.75 <- runCount(tree,data,output,1.75)
#gain2 <- runCount(tree,data,output,2)
#gain2.5 <- runCount(tree,data,output,2.5)
#gain_df <- rbind(gain1,gain1.25,gain1.5,gain1.75,gain2,gain2.5)

#visualize count output

family <- read_tsv(file.path(outdir,'Bxy_count.FAMILY'))
present <- read_tsv(file.path(outdir,'Bxy_count.PRESENT'))
tree <- read.tree(file=file.path(outdir,'Bxy.tree'))
colnames(wag)
wag <- as.data.frame(wag) %>% 
  select(-CHANGE) %>%
  rename(label=node) %>%
  select(label,family_gain,family_loss) %>%
  droplevels.data.frame()

tree_dat = tree %>% 
  treeio::as_tibble() %>% 
  left_join(wag,by='label') %>%
  mutate(family_loss_range = lapply(family_loss, function(x) c(0, x/10000))) %>%
  mutate(family_gain_range = lapply(family_gain, function(x) c(0, x/10000))) %>% 
  left_join(metadata,by=c('label'='isolate')) %>%
  as.treedata()

tree_data <- as_tibble(tree_dat)
tree_plot <- ggtree(tree_dat)  + 
  geom_tiplab(align=TRUE,linetype='dashed', linesize=.3) +
  #geom_range(range='family_gain_range', color="red", size = 1, alpha = 0.5) +
  #geom_range(range='family_loss_range', color="blue", size = 1, alpha = 0.5) +
  geom_text2(aes(x=branch, label=family_gain, subset=(family_gain>400)), vjust=-.1, color='firebrick') +
  geom_text2(aes(x=branch, label=family_loss, subset=(family_gain>400)), vjust=1.3, color='blue')
tree_plot
ggtree(tree) %<+% tree_data +
  geom_tippoint(aes(color=host_site)) +
  geom_text2(aes(x=branch, label=family_gain, subset=(family_gain>400)), size = 3, vjust=-.1, color='firebrick') +
  geom_text2(aes(x=branch, label=family_loss, subset=(family_gain>400)), size = 3, vjust=1.3, color='blue')


cluster99 <- read.table(paste0('results/pangenome/Bacteroides_xylanisolvens/roary/isolate_clusters_99.txt'))
cluster99 <- cluster99 %>% mutate(isolate = str_replace_all(isolate,'[_-]','.'))
metadata <- metadata %>% left_join(cluster99,by='isolate')
metadata <- metadata %>% 
  mutate(isolate = str_replace_all(isolate,'[_-]','.'),
         site = recode(site,'USA: Cambridge'='USA','USA:Boston'='USA',
                       'China: Shenzhen'='China','USA:Seattle'='USA',
                       'not applicable'='siteUnknown','missing'='siteUnknown')) %>%
  unite(host_site, host, site, sep = "_", remove = FALSE)

getNode <- function(NumCluster){
  cl <- metadata %>% filter(cluster == NumCluster) %>% select(isolate,host_site) 
  if (length(cl$isolate)>1) {
    cl$isolate <- stringr::str_replace_all(cl$isolate,'[_-]','.')
    node <- getMRCA(tree,cl$isolate)
    #print(paste('cluster',NumCluster))
    #print(cl$isolate)
    print(paste('node',node,length(cl$isolate),'isolates'))
    print(unique(cl$host_site))
    return(node)
  } else { return(0)}
}
nodes <- lapply(sort(unique(cluster99$cluster)),getNode)
sort(nodes)
tree_plot <- ggtree(tree_dat) +
  geom_tiplab() +
  geom_text2(aes(x=branch, label=family_gain), vjust=-.5, color='firebrick') +
  geom_text2(aes(x=branch, label=family_loss), vjust=1.3, color='blue') 
tree_plot

tree_plot_collapse <- tree_plot %>% 
  collapse(141) %>% 
  collapse(181) %>% 
  collapse(183) %>% 
  collapse(193) %>% 
  collapse(205) %>% 
  collapse(212) %>% 
  collapse(217) %>% 
  collapse(240) %>% 
  collapse(249) %>% 
  collapse(261)  +
  geom_point2(aes(subset=(node==141)), shape=21, size=15, fill='royalblue') +
  geom_text2(aes(subset= (node==141)),label="29 human_USA isolates",hjust = -.25) +
  geom_point2(aes(subset=(node==181)), shape=21, size=2, fill='royalblue')+
  geom_text2(aes(subset= (node==181)),label="2 human_USA isolates",hjust = -.25) +
  geom_point2(aes(subset=(node==183)), shape=21, size=8, fill='royalblue')+
  geom_text2(aes(subset= (node==183)),label="8 human_USA isolates",hjust = -.25) +
  geom_point2(aes(subset=(node==193)), shape=21, size=3, fill='royalblue')+
  geom_text2(aes(subset= (node==193)),label="3 human_USA isolates",hjust = -.25) +
  geom_point2(aes(subset=(node==205)), shape=21, size=5, fill='royalblue') +
  geom_text2(aes(subset= (node==205)),label="5 human_USA isolates",hjust = -.25) +
  geom_point2(aes(subset=(node==212)), shape=21, size=4, fill='royalblue')+
  geom_text2(aes(subset= (node==212)),label="4 human_USA isolates",hjust = -.25) +
  geom_point2(aes(subset=(node==217)), shape=21, size=11, fill='purple')+
  geom_text2(aes(subset= (node==217)),label="Mixed-host captive clade",hjust = -.1) +
  geom_point2(aes(subset=(node==240)), shape=21, size=8, fill='royalblue')+
  geom_text2(aes(subset= (node==240)),label="8 human_USA isolates",hjust = -.1) +
  geom_point2(aes(subset=(node==249)), shape=21, size=10, fill='royalblue') +
  geom_text2(aes(subset= (node==249)),label="13 human_USA isolates",hjust = -.1) +
  geom_point2(aes(subset=(node==261)), shape=21, size=2, fill='royalblue') +
  geom_text2(aes(subset= (node==261)),label="2 human_USA isolates",hjust = -.1) 

tree_plot_collapse
