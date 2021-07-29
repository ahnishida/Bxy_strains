library(nice)
library(tidyverse)
library("Rcpp")
library(Biostrings)
library(ape)
library(harrietr)
library(vegan)
library(parallel)
library(data.table)
source('scripts/genomic_island_function.R')
set.my.priority(priority=15)

#produces output for figures showing pairwise gene content differences,
# of events and mean island size

final_outdir = 'results/analyses_input_files/pairwise_gene_islands'
dir.create(final_outdir)
window_size = 5


get_phylo_independent_comps = function(species,cutoff) {
  #identify phylogenetically independent pairs of strains
  #given the min cutoff tree dist
  
  #read in input files
  print(species)
  dir = paste0('results/pangenome/',species)
  metadata = read_tsv(file.path(dir,'metadata.txt'),col_types = cols())
  pres_abs <- read_csv(file.path(dir,'roary_nosplitparalogs/gene_presence_absence.csv'),col_types = cols())
  tree = read.tree(file.path(dir,paste0(species,'.tre')))
  
  #outdir
  pw_outdir = file.path(dir,'pairwise_gene_gain')
  print(pw_outdir)
  outgroup = metadata %>% filter(taxonomy_Species!=species) %>% pull(isolate)
  tree = drop.tip(tree,outgroup)
  tree_dist <- cophenetic.phylo(tree)
  
  pw_df <- melt_dist(tree_dist) %>% 
    dplyr::rename(tree_dist=dist) %>% 
    filter(tree_dist > as.numeric(cutoff))
  
  #add metadata to pw comparisons
  metadata_subset <- metadata %>% dplyr::select(isolate,host,sample,taxonomy_Species)
  pw_df_meta  <- pw_df  %>%    
    left_join(metadata_subset,by=c('iso1'='isolate')) %>%  #add metadata for individual1 and 
    left_join(metadata_subset,by=c('iso2'='isolate'),suffix = c(".iso1", ".iso2")) %>%
    filter(host.iso1 != 'missing',host.iso2 != 'missing') %>%
    mutate(diff_host = if_else(host.iso1 != host.iso2, 1, 0)) %>%
    mutate(ape_host = if_else(host.iso1 %in% c('chimpanzee','bonobo','orangutan','gorilla')|
                                host.iso2 %in% c('chimpanzee','bonobo','orangutan','gorilla'),1,0)) %>%
    arrange(desc(diff_host),tree_dist) %>% 
    mutate(comp = paste(iso1,iso2,sep='_')) 
  
  identify_taxa <- function(iso1,iso2) {
    #determine all isolates have the same MRCA as two input isolates
    node = getMRCA(tree,c(iso1,iso2))
    subtree = extract.clade(tree, node)
    return(subtree$tip.label)
  }
  
  taxa = c()
  rows = c()
  for (row in 1:nrow(pw_df_meta)) { #loop through tree distm
    iso1=pw_df_meta[row,1]
    iso2=pw_df_meta[row,2]
    subtree_taxa = identify_taxa(iso1,iso2)
    #add row to list if there's no phylogenetic overlap with isolate comparison already present in the list
    if (length(intersect(subtree_taxa,taxa))==0){ 
      taxa = c(taxa,subtree_taxa)  #keep track of isolates included in comparisons            
      rows=c(rows,row)
      print(row)
    }
  }
  pw_df_PI = pw_df_meta[rows,]
  return(pw_df_PI) 
}

Bxy = get_phylo_independent_comps('Bacteroides_xylanisolvens',.008)
Bov = get_phylo_independent_comps('Bacteroides_ovatus',.0001)
Bfr = get_phylo_independent_comps('Bacteroides_fragilis',.0001)
Bth = get_phylo_independent_comps('Bacteroides_thetaiotamicron',.0001)


#sanity check to be sure the pairs are independent
test = pw_df_PI %>% 
  mutate(pair = 1:nrow(pw_df_PI)) %>% 
  tidyr::pivot_longer(cols=c(iso1,iso2), names_to='1vs2') %>%
  dplyr::rename('iso'='value') %>%
  select(iso,everything())
ggtree(tree) %<+% test + 
  geom_tippoint(aes(subset=(label%in%test$iso))) +
  geom_tiplab(aes(subset=(label%in%test$iso)))

identify_pw_50strain = function(species) {
  #reads in relevant data files
  #identify 50strains based on hclust of tree distances 
  #returns flattened dist matrix of all 1225 pw comparison and metadata
  
  #read in input files
  print(species)
  dir = paste0('results/pangenome/',species)
  metadata = read_tsv(file.path(dir,'metadata.txt'),col_types = cols())
  pres_abs <- read_csv(file.path(dir,'roary_nosplitparalogs/gene_presence_absence.csv'),col_types = cols())
  tree = read.tree(file.path(dir,paste0("phylogeny/RAxML_bipartitions.",species)))
  #outdir
  pw_outdir = file.path(dir,'pairwise_gene_gain')
  dir.create(pw_outdir)
  #create distance matrix from tree
  outgroup = metadata %>% filter(taxonomy_Species!=species) %>% pull(isolate)
  tree = drop.tip(tree,outgroup)
  tree_dist <- cophenetic.phylo(tree)
  
  #dereplication closely related strains
  set.seed(123)
  hc = hclust(as.dist(tree_dist))
  cl = cutree(hc,k=50)
  df = data.frame(cl) %>% 
    rownames_to_column(var='isolate') %>%
    group_by(cl) %>% sample_n(1)
  pw_df <- melt_dist(tree_dist) %>% 
    dplyr::rename(tree_dist=dist) %>%
    filter(iso1 %in% df$isolate, iso2 %in% df$isolate)
  pw_df = pw_df %>% 
    mutate(comp = paste(iso1,iso2,sep='_')) 
  #add metadata to pw comparisons
  metadata_subset <- metadata %>% dplyr::select(isolate,host,sample,taxonomy_Species)
  pw_df  <- pw_df  %>%    
    left_join(metadata_subset,by=c('iso1'='isolate')) %>%  #add metadata for individual1 and 
    left_join(metadata_subset,by=c('iso2'='isolate'),suffix = c(".iso1", ".iso2")) 
  return(pw_df)
}
Bxy_pw = identify_pw_50strain('Bacteroides_xylanisolvens')
Bov_pw = identify_pw_50strain('Bacteroides_ovatus')
Bfr_pw = identify_pw_50strain('Bacteroides_fragilis')
Bth_pw = identify_pw_50strain('Bacteroides_thetaiotaomicron')

pairwise_gene_gain_runner = function(iso1,iso2,species) {
  #set input filepaths to run pairwise_gene_gain
  dir = paste0('results/pangenome/',species)
  metadata = read_tsv(file.path(dir,'metadata.txt'),col_types = cols())
  pres_abs <- read_csv(file.path(dir,'roary_nosplitparalogs/gene_presence_absence.csv'),col_types = cols())
  pw_outdir = file.path('results/pangenome',species,'pairwise_gene_gain')
  iso1_old = metadata$isolate_old[metadata$isolate==iso1]
  iso2_old = metadata$isolate_old[metadata$isolate==iso2]
  iso1_gff_file = paste0('results/pangenome/prokka/',iso1_old,'/',iso1_old,'.gff')
  iso2_gff_file = paste0('results/pangenome/prokka/',iso2_old,'/',iso2_old,'.gff')
  #print(iso1)
  #print(iso2)
  pairwise_gene_gain(
    #loaded from functions script
    pres_abs = pres_abs,
    iso1=iso1,
    iso2=iso2,
    iso1_gff_file = iso1_gff_file,
    iso2_gff_file = iso2_gff_file,
    window_size = window_size,
    pw_outdir = pw_outdir)
}
pairwise_gene_gain_runner(Bxy_pw$iso1[1],Bxy_pw$iso2[1],'Bacteroides_xylanisolvens')

run_pw_on_left_to_do <- function(pw_df) {
  #subset to those pairwise comparisons left to do
  #run pairwise_gene_gain_runner
  species = pw_df$taxonomy_Species.iso1[1]
  pw_outdir = file.path('results/pangenome',species,'pairwise_gene_gain')
  dir.create(pw_outdir)
  files = list.files(file.path(pw_outdir),'summary.txt')
  files = paste0(pw_outdir,'/',files)
  completed = sapply(strsplit(files, "/"), "[", 5)
  completed = sapply(strsplit(completed,"_window"),"[", 1)
  completed = ''
  to_do = pw_df %>%
      filter(!comp %in% completed)
  print('left to do:')
  print(nrow(to_do))
  if (nrow(to_do)>0){
    mcmapply(pairwise_gene_gain_runner,to_do$iso1,to_do$iso2,species,mc.cores=20)
    }
  }
run_pw_on_left_to_do(Bxy_pw)
run_pw_on_left_to_do(Bov_pw)
run_pw_on_left_to_do(Bfr_pw)
run_pw_on_left_to_do(Bth_pw)

aggregate_summary_files <- function(pw_df) {
  #concatenate individual pairwise comparison summary files  
  species = pw_df$taxonomy_Species.iso1[1]
  pw_outdir = file.path('results/pangenome',species,'pairwise_gene_gain')
  summary_files = unlist(paste0(file.path(pw_outdir),'/',pw_df$comp,'_window5_summary.txt'))
  summary_all = rbindlist(lapply(summary_files,function(x){read.csv(x,sep='\t')}))
  summary_all = pw_df  %>%  #join with tree dist
    mutate(comp=paste0(iso1,'_',iso2)) %>% 
    left_join(summary_all,by=c('iso1','iso2','comp')) %>%
    mutate(mean_treedist = mean(tree_dist),
           sd_treedist = sd(tree_dist),
           tree_dist_norm = ((tree_dist-mean_treedist)/sd_treedist)) 
  write_tsv(summary_all,file=file.path(final_outdir,paste0(species,'_summary.txt')))
  return(summary_all)
}
Bxy_sum = aggregate_summary_files(Bxy_pw)
Bov_sum = aggregate_summary_files(Bov_pw)
Bfr_sum = aggregate_summary_files(Bfr_pw)
Bth_sum = aggregate_summary_files(Bth_pw)

aggregate_gff_island_files <- function(pw_df) {
  #concatenate individual pairwise comparison gff_island files, which show all genes 
  #in every event, needed for seeing distribution of events by size/genes transferred
  species = pw_df$taxonomy_Species.iso1[1]
  pw_outdir = file.path('results/pangenome',species,'pairwise_gene_gain')
  gff_island_files = unlist(paste0(file.path(pw_outdir),'/',pw_df$comp,'_window5_island_gff.txt'))
  gff_island_all = rbindlist(lapply(gff_island_files,function(x){read.csv(x,sep='\t')}))
  gff_island_all = gff_island_all %>%   #join with tree dist
    left_join(pw_df,by=c('iso1','iso2','comp')) 
  gff_island_all$comp[1]
  pw_df$comp[1]
  intersect(pw_df$iso2,gff_island_all$iso2)
  print(nrow(gff_island_all))
  print(head(gff_island_all))
  write_tsv(gff_island_all,file=file.path(final_outdir,paste0(species,'_events.txt')))
}


#determine cutoff to use
Bac = Bxy_sum %>% add_row(Bov_sum) %>% add_row(Bfr_sum) %>% add_row(Bth_sum)
Bac_2sd = Bac %>% 
  select(taxonomy_Species.iso1,mean_treedist,sd_treedist) %>% distinct() %>%
  mutate(mean_treedist-(2*sd_treedist))

Bxy_pi = get_phylo_independent_comps('Bacteroides_xylanisolvens',0.01214)
run_pw_on_left_to_do(Bxy_pi)
aggregate_gff_island_files(Bxy_pi)
Bov_pi = get_phylo_independent_comps('Bacteroides_ovatus',0.0138)
run_pw_on_left_to_do(Bov_pi)
aggregate_gff_island_files(Bov_pi)
Bth_pi = get_phylo_independent_comps('Bacteroides_thetaiotaomicron',0.0100)
run_pw_on_left_to_do(Bth_pi)
aggregate_gff_island_files(Bth_pi)
Bfr_pi = get_phylo_independent_comps('Bacteroides_fragilis',0.00548)
run_pw_on_left_to_do(Bfr_pi)
aggregate_gff_island_files(Bfr_pi)

 

