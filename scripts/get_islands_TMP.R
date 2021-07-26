library(nice)
library(tidyverse)
library("Rcpp")
library(Biostrings)
library(ape)
library(harrietr)
library(vegan)
library(parallel)
library(data.table)

set.my.priority(priority = 15)
source('scripts/genomic_island_function.R')

species = 'Bacteroides_ovatus'

run_all_pairwise_gene_gain = function(species) {

  dir = paste0('results/pangenome/',species)
  metadata = read_tsv(file.path(dir,'metadata.txt'),col_types = cols())
  outgroup = metadata %>% filter(taxonomy_Species!=species) %>% pull(isolate)
  pres_abs <- read_csv(file.path(dir,'roary_nosplitparalogs/gene_presence_absence.csv'),col_types = cols())
  pw_outdir = file.path(dir,'pairwise_gene_gain')
  print(pw_outdir)
  window_size = 5
  tree_file = file.path(dir,paste0("phylogeny/RAxML_bipartitions.",species))
  tree = read.tree(tree_file)
  tree = drop.tip(tree,outgroup)
  tree_dist <- cophenetic.phylo(tree)
  tree_distm <- melt_dist(tree_dist) %>% 
    dplyr::rename(tree_dist=dist) 
  
  #dereplication closely related strains
  #set.seed(123)
  #hc = hclust(as.dist(tree_dist))
  #cl = cutree(hc,k=5)
  #df = data.frame(cl) %>% 
  #  rownames_to_column(var='isolate') %>%
  #  group_by(cl) %>% sample_n(1)
  tree_distm_filt <- tree_distm %>%
    filter(iso1 %in% df$isolate, iso2 %in% df$isolate)
  
  #exclude outgroup from comparisons
  metadata_subset <- metadata %>% dplyr::select(isolate,host,sample,taxonomy_Species)
  tree_distm  <- tree_distm  %>%    
    left_join(metadata_subset,by=c('iso1'='isolate')) %>%  #add metadata for individual1 and 
    left_join(metadata_subset,by=c('iso2'='isolate'),suffix = c(".iso1", ".iso2")) %>%
    filter(taxonomy_Species.iso1==species,taxonomy_Species.iso2==species)

  #subset to those left to do
  #files = list.files(file.path(pw_outdir),'summary.txt')
  #files = paste0(pw_outdir,'/',files)
  #print(length(files))
  #complete = matrix(NA,ncol=15,nrow=length(files))
  #for (i in 1:length(files)){
  #  one_row = read.csv(files[i],sep='\t')
  #  complete[i,]=as.matrix(one_row[1,])
  #}
  #complete = as.data.frame(complete)
  #colnames(complete)=colnames(one_row)
  
  #to_do = tree_distm %>% 
  #    left_join(complete,by=c('iso1','iso2')) %>% 
  #    filter(is.na(window_size))
  #print('left to do:')
  #print(nrow(to_do))

  pairwise_gene_gain_runner = function(iso1,iso2) {
      #print(species)
      iso1_old = metadata$isolate_old[metadata$isolate==iso1]
      iso2_old = metadata$isolate_old[metadata$isolate==iso2]
      iso1_gff_file = paste0('results/pangenome/prokka/',iso1_old,'/',iso1_old,'.gff')
      iso2_gff_file = paste0('results/pangenome/prokka/',iso2_old,'/',iso2_old,'.gff')
      #print(iso1)
      #print(iso2)
      pairwise_gene_gain(
        pres_abs = pres_abs,
        iso1=iso1,
        iso2=iso2,
        iso1_gff_file = iso1_gff_file,
        iso2_gff_file = iso2_gff_file,
        window_size = window_size,
        pw_outdir = pw_outdir)
  }
  
  #for (i in nrow(to_do)){
  #  row = complete[i,]
  #  pairwise_gene_gain_runner(as.character(row$iso1),as.character(row$iso2))
   # }
  #pairwise_gene_gain_runner(to_do$iso1,to_do$iso2,mc.cores = 20)
  
  #combine results into one
  files = list.files(file.path(pw_outdir),'summary.txt')
  files = paste0(pw_outdir,'/',files)
  print(length(files))
  
  #COMMENT OUT
  #files=files[1:100]
  
  res = matrix(NA,ncol=15,nrow=length(files))
  for (i in 1:length(files)){
    one_row = read.csv(files[i],sep='\t')
    res[i,]=as.matrix(one_row[1,])
  }
  res = as.data.frame(res)
  colnames(res)=colnames(one_row)

  #join with tree dist
  tree_distm_pw = tree_distm %>% 
    left_join(res,by=c('iso1','iso2'))%>%
    filter(!is.na(window_size))
  write_tsv(tree_distm_pw,file=file.path(dir,'pairwise_gene_gain.txt'))
  print(paste(species,'complete'))
  return(tree_distm_pw)
}

bov = run_all_pairwise_gene_gain('Bacteroides_ovatus')
bfr = run_all_pairwise_gene_gain('Bacteroides_fragilis')
bth = run_all_pairwise_gene_gain('Bacteroides_thetaiotaomicron')
bth = run_all_pairwise_gene_gain('Bacteroides_xylanisolvens')





