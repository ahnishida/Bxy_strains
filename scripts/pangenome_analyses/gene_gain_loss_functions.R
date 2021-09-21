set.seed(1324)

library(Biostrings)
library(harrietr)
library(vegan)
library(tidyverse)
library(ape)

#script contains all functions to cluster genes into genomic islands and 
#conduct pairwise genome comparisons, to test out functions download test data
#and uncomment function test runs.

##test inputs
#setwd(#path_to_repo)
#test_dir = 'scripts/test_data'
#list.files(test_dir)
#pres_abs <- read_csv(file.path(test_dir,'gene_presence_absence.csv'),col_types = cols())
#captive_convergent_genes = readRDS(file.path(test_dir,'captive_convergent_genes.RDS'))

parse_gff_window <- function(gff_file,window_size){
  #reads in gff file determines start and end positions of genes that are 
  #upstream and downstream on the same contigs based on the window size provided
  gff <- read_tsv(gff_file, comment='##',col_names =F,col_types = cols())
  colnames(gff) <- c('contig','source','feature',
                     'start','end','score','strand','frame','attribute')
  gff <- gff %>% 
    separate(col=c('attribute'),sep=';',into = c('Gene.ID'),extra='drop') %>%
    separate(col=c('Gene.ID'),sep='=',into = c(NA,'Gene.ID'),extra='drop') %>% #pull out GeneID
    mutate(start = as.numeric(start),
           end = as.numeric(end)) %>%
    filter(!is.na(start)) %>% #get rid of the fasta lines at the end of file
    as.data.frame() #need to assert this after separating/mutating/grouping, otherwise get mixed class object
  contig_df = gff %>% group_by(contig) %>%  #get start and end positions of all contigs
    summarise(contig_start = min(start,na.rm=T),
              contig_end = max(end,na.rm=T)) %>% 
    as.data.frame()
  gff_window = gff %>% 
    group_by(contig) %>% #performs action on each contig
    left_join(contig_df,by='contig') %>% #add contig start and end info
    mutate(upstream = lag(start,n=window_size), #
           downstream = lead(end,n=window_size),
           upstream = ifelse(is.na(upstream),contig_start,upstream),
           downstream = ifelse(is.na(downstream),contig_end,downstream)) %>% 
    as.data.frame()
  return(gff_window)
}

#gff_window = parse_gff_window(file.path(test_dir,'P17.A3.gff'),window_size=5)

get_genomic_island <- function(pres_abs,list_of_genes,isolate,gff_file,window_size,outfile) {
  #given a list of Genes as found in roary presence abs table output, parse isolate gff to determine
  #how many genomic islands the genes fall into given a specified window size, 
  #for example with a window size of 5, 2 genes that are within 5 genes of eachother are
  #considered part of the same island
  df = dplyr::select(pres_abs,Gene,all_of(isolate)) %>% 
    tidyr::separate_rows(isolate, sep = '[\t]',convert = FALSE) %>% #subset to isolate of interest split paralogs into own rows
    filter(Gene %in% list_of_genes) %>% #subset to genes specified
    as.data.frame() 
  #print('gene families more than one copy')
  #print(df$Gene[table(df$Gene)>1]) #some genes might be paralogs
  gff_window = parse_gff_window(gff_file,window_size) 
  gff_filt = gff_window %>% 
    right_join(df,by=c('Gene.ID'=all_of(isolate))) %>% #add Gene info to gff_window and subset to genes of interest
    filter(!is.na(Gene.ID)) 
  if (nrow(gff_filt) == 0) {
    print('Warning no Genes from list_of_genes found in isolate')
    return(data.frame())
    
  }
  if (nrow(gff_filt) == 1) {
    print('Warning single Gene from list_of_genes found in isolate')
    gff_islands = gff_filt %>% mutate(size = 1,
                                      rank = 1,
                                      gene_span = paste0(Gene,'_',Gene),
                                      name = paste0(gene_span,'_',rank))
  } else {
    #get info of previous gene
    gff_filt$prev_gene = lag(gff_filt$Gene.ID,1) 
    gff_filt$prev_contig = lag(gff_filt$contig,1)
    gff_filt$prev_upstream = lag(gff_filt$upstream,1) 
    gff_filt$prev_downstream = lag(gff_filt$downstream,1) 
    gff_filt= gff_filt %>% mutate(in_window = ifelse((contig==prev_contig)& #assess if in same window
                                                       (start>=prev_upstream)&
                                                       (end<=prev_downstream),1,0))  
    gff_filt$in_window[1] <- gff_filt$in_window[2] #first gene will always be NA so apply value from second gene
    gff_filt = as.data.frame(gff_filt)
    n=1
    rank = c()
    for (x in gff_filt$in_window){ #iterate over in_window vales
      if(x==0){n=n+1} #increase rank by 1 everytime there is a 0 
      rank <- c(rank,n)
    }
    gff_filt$rank <- rank
    gff_islands = gff_filt %>% 
      group_by(rank) %>%
      add_count(rank,name = 'size') %>% #determine island size
      mutate(gene_start = Gene[start==min(start)], #determine gene most upstream
             gene_end = Gene[end==max(end)], #determine gene most upstream
             gene_span = 
               paste(pmin(gene_start, gene_end), #label gene span as sorted order some seqs are in rc orientation
                     pmax(gene_start, gene_end),sep='_'),
             name = paste0(gene_span,'_',rank)) %>% #get unique name for gene island
      dplyr::select(-gene_start,-gene_end) %>% as.data.frame() #remove excess variables
  }
  #if outdir included in function output island gff df
  if(!missing(outfile)){
    dir.create(dirname(outfile))
    write_tsv(gff_islands,file=outfile)
  }
  closing_message = paste('Found ',as.character(length(unique(gff_islands$Gene))),
                          'of',as.character(length(unique(list_of_genes))),'genes families provided.',
                          as.character(length(gff_islands$Gene.ID)),'gene copies belonging to this families found on',
                          as.character(length(unique(gff_islands$name))), 'genomic islands using window size',
                          as.character(window_size),sep = ' ')
  print(closing_message)
  return(gff_islands)
}

#P17.A3_GItest = get_genomic_island(pres_abs=pres_abs,
#                               list_of_genes=captive_convergent_genes,
#                               isolate='P17.A3',
#                               gff_file=file.path(test_dir,'P17.A3.gff'),
#                               window_size = 5,
#                               outfile = file.path(test_dir,'output/P17.A3_window5_GIoutput.txt'))

output_island_fasta <- function(isolate,genomic_islands_outfile,gene_span_input,
                                gff_file,updown_size,fna_file,island_fasta_outfile) {
  #outputs cluster fna after running extract_clusters functions 
  #define how many genes upstream and downstream to extract
  #set to 0 to just get the cluster
  print(c(isolate,gene_span_input)) 
  
  #get contig and list of geneIDs on genomic island
  genomic_islands_output = read_tsv(genomic_islands_outfile,col_types = cols())
  island_gff = genomic_islands_output %>% filter(gene_span == gene_span_input) 
  island_GeneIDs = island_gff %>% pull(Gene.ID)
  island_contig = island_gff$contig[1]
  
  #define start and end position on contig for cluster given the windown size
  island_updown = parse_gff_window(gff_file,updown_size) %>%
    filter(Gene.ID %in% island_GeneIDs) 
  start = min(island_updown$upstream)
  end = max(island_updown$downstream)
  
  #extract x genes upstream and downstream of island
  fna = readDNAStringSet(fna_file)
  island_fna = fna[names(fna)==island_contig]
  names(island_fna) = paste(isolate,gene_span_input,sep='_')
  island_fna = subseq(island_fna, start=start, end=end)
  system(paste0('mkdir -pv ',dirname(island_fasta_outfile)))
  writeXStringSet(island_fna,file=island_fasta_outfile)
  
  closing_message = paste('output genomic island',gene_span_input,'containing',as.character(length(island_GeneIDs)),
                          'genes plus',
                          as.character(updown_size),'genes on either side, total length is',as.character(end - start),
                          sep=' ')
  return(closing_message)
}

#output_island_fasta(
#  isolate='P17.A3',
#  genomic_islands_outfile = file.path(test_dir,'output/P17.A3_window5_GIoutput.txt'),
#  island_fasta_outfile = file.path(test_dir,'output/fna/P17.A3_group_20589_group_20604_window5.fna'),
#  gene_span = 'group_20589_group_20604',
#  updown_size = 10,
#  fna_file = file.path(test_dir,'P17.A3.fna'),
#  gff_file = file.path(test_dir,'P17.A3.gff')
#)


get_copynum_table = function(pres_abs,iso1,iso2) {
  #reads in pres abs table to calc copy numbers and merges with subset table
  geneid_table = pres_abs %>% 
    dplyr::select(Gene,all_of(iso1),all_of(iso2))
  colnames(geneid_table) <- c('Gene','Gene.ID_iso1','Gene.ID_iso2')
  get_copy_num <- function(x) {as.numeric(str_count(x, pattern = "_"))}
  copynum_table <- geneid_table  %>%
    mutate_at(vars(-Gene),get_copy_num) 
  colnames(copynum_table) <- c('Gene','copynum_iso1','copynum_iso2')
  copynum_table = copynum_table %>% 
    mutate(iso1=iso1,
           iso2=iso2,
           copynum_iso1 = tidyr::replace_na(copynum_iso1,0),
           copynum_iso2 = tidyr::replace_na(copynum_iso2,0),
           diff = ifelse(copynum_iso1>copynum_iso2,iso1,
                         ifelse(copynum_iso2>copynum_iso1,iso2,NA)),
           copy_diff = abs(copynum_iso1-copynum_iso2)) 
  copynum_table = geneid_table %>% 
    left_join(copynum_table,by='Gene') %>% 
    as.data.frame()
  return(copynum_table)
}
#copynum_table = get_copynum_table(pres_abs,'P21.6G','P21.6D')


#Count output
get_mrca = function(countFAMILYoutput_file,countTREE_file,copynum_table) {
  #Select mrca node of iso1 and iso2 from count output table showing gene family 
  #pres/abs across the phylogeny, labels if gene copy number in iso1 or iso2 represent
  #a loss/gain/nochange from mrca
  iso1 = copynum_table$iso1[1]
  iso2 = copynum_table$iso2[1]
  countTREE = read.tree(countTREE_file)
  tree_data = countTREE %>% 
    treeio::as_tibble() 
  nodeNumber = ape::getMRCA(countTREE,c(iso1,iso2))
  nodeName = tree_data %>% filter(node==nodeNumber) %>% pull(label)
  countFAMILYoutput = read_tsv(countFAMILYoutput_file,col_types = cols())
  mrca_gene_counts = countFAMILYoutput %>% select(Gene,all_of(nodeName)) %>% as.data.frame()
  colnames(mrca_gene_counts) = c('Gene','mrca')
  copynum_mrca_table = copynum_table %>% 
    left_join(mrca_gene_counts,by='Gene') 
  copynum_mrca_table = copynum_mrca_table %>% 
    mutate(MRCA_iso1 = if_else((copynum_iso1 == mrca), 'nochange',
                               if_else(copynum_iso1 > mrca,'gain','loss'))) %>% 
    mutate(MRCA_iso2 = if_else((copynum_iso2 == mrca), 'nochange',
                               if_else(copynum_iso2 > mrca,'gain','loss'))) 
  return(copynum_mrca_table)
}

#countFAMILYoutput_file = file.path(test_dir,"countOutput_Bacteroides_xylanisolvens_gainpenalty2.FAMILY")
#countTREE_file = file.path(test_dir,"count_Bacteroides_xylanisolvens.tre")
#copynum_mrca_table = get_mrca(countFAMILYoutput_file,
#                              countTREE_file,
#                            get_copynum_table(pres_abs,'P15.E10','P21.6D'))


get_copy_diff = function(pres_abs,
                         subset_mrca_table_output,
                         iso1,iso2,iso1_gff_file,iso2_gff_file,window_size) {
  #determine copy number differences between gene families and how many genomic island these
  #gained genes reside on given a specified window size

  copydiff_table = subset_mrca_table_output %>% 
    filter(!is.na(diff)) %>% #genes that show some copy number difference
    mutate(Gene.ID = ifelse(diff==iso1,Gene.ID_iso1,Gene.ID_iso2)) %>%
    dplyr::select(-Gene.ID_iso1,-Gene.ID_iso2) %>%
    separate_rows(Gene.ID, sep = '[\t]') %>%
    as.data.frame()
  
  #divide up gene ids into four categories and see if genes are adjacent 
  iso1_gain =  copydiff_table$Gene[copydiff_table$geneGainLoss=='iso1_gain']  
  iso1_loss =  copydiff_table$Gene[copydiff_table$geneGainLoss=='iso1_loss']
  iso2_gain =  copydiff_table$Gene[copydiff_table$geneGainLoss=='iso2_gain']
  iso2_loss =  copydiff_table$Gene[copydiff_table$geneGainLoss=='iso2_loss']
  print('iso1_gain_events')
  iso1_gain_islands = get_genomic_island(pres_abs=pres_abs,
                                    list_of_genes=iso1_gain,
                                    isolate=iso1,
                                    gff_file=iso1_gff_file,
                                    window_size=window_size) 
  print('iso1_loss_events')
  iso1_loss_islands = get_genomic_island(pres_abs=pres_abs,
                                         list_of_genes=iso1_loss,
                                         isolate=iso2,
                                         gff_file=iso2_gff_file,
                                         window_size=window_size) 
  print('iso2_gain_events')
  iso2_gain_islands = get_genomic_island(pres_abs=pres_abs,
                                         list_of_genes=iso2_gain,
                                         isolate=iso2,
                                         gff_file=iso2_gff_file,
                                         window_size=window_size) 
  print('iso2_loss_events')
  iso2_loss_islands = get_genomic_island(pres_abs=pres_abs,
                                         list_of_genes=iso2_loss,
                                         isolate=iso1,
                                         gff_file=iso1_gff_file,
                                         window_size=window_size)  
  islands = rbind(iso1_gain_islands,iso1_loss_islands,iso2_gain_islands,iso2_loss_islands)  %>%
    as.data.frame() %>% 
    dplyr::select(Gene.ID,contig,name,size,rank)
  copydiff_island = copydiff_table %>% 
    left_join(islands, by='Gene.ID') %>%
    as.data.frame() 
  return(copydiff_island)
  
}

#iso1='P21.6D'
#iso2='P17.A3'
#iso1_gff_file=file.path(test_dir,'P21.6D.gff')
#iso2_gff_file=file.path(test_dir,'P17.A3.gff')

#copynum_mrca_table = get_mrca(countFAMILYoutput_file,
#                              countTREE_file,
#                              get_copynum_table(pres_abs,iso1,iso2))

##filter out HGGs where one of the isolate copy number is not equal to mcra
#sum(copynum_mrca_table$mrca)
#subset_mrca_table = copynum_mrca_table %>%
#  filter(MRCA_iso1 == 'nochange' | MRCA_iso2 == 'nochange') %>% #genes that match mrca 
#  mutate(geneGainLoss = ifelse(MRCA_iso1 == 'gain', 'iso1_gain',
#                               ifelse(MRCA_iso1 == 'loss', 'iso1_loss',
#                                      ifelse(MRCA_iso2 == 'gain','iso2_gain',
#                                             ifelse(MRCA_iso2 == 'loss', 'iso2_loss','nochange'))))) 
#sum(subset_mrca_table$mrca)
#copydiff_island = get_copy_diff(pres_abs=pres_abs,
#                                subset_mrca_table_output=subset_mrca_table,
#                                iso1=iso1,
#                                iso2=iso2,
#                                iso1_gff_file=iso1_gff_file,
#                                iso2_gff_file=iso2_gff_file,
#                                window_size=1)

is_dup <- function(copydiff_island) {
  #adds is_dup column to df 
  #determine islands that may represent duplication events
  #criterion all HGG on island have to have multiple copies
  dup = copydiff_island %>% 
    group_by(name,size,Gene) %>% 
    tally() %>%
    as.data.frame() %>%
    filter(n>1) %>%
    group_by(name,Gene) %>%
    mutate(dup_size = sum(n)) %>%
    filter(size == dup_size) %>%
    as.data.frame()
  copydiff_island = copydiff_island %>% mutate(is_dup = ifelse(name %in% dup$name, 1,0))
  return(copydiff_island)
}
#copydiff_island = is_dup(copydiff_island)

  
pick_genes_on_clust_size <- function(copydiff_island,diff){
  #in cases where gene family is present in both isolate but copy number differs
  # we select geneID that was gain/lost by optimizing the few number of transfers 
  #select the gene copy that is located on the largest island
  df = copydiff_island %>% filter(copy_diff == UQ(diff))
  df <- df  %>%
    group_by(Gene) %>%
    slice_max(n=diff,order_by=size) #gets biggest cluster
  df  <- df  %>%
    group_by(Gene) %>%
    slice_min(n=diff,order_by=Gene.ID) #picks just one if hits to same cluster,
  return(df)                           #usually the case with single duplication events
}
#res = pick_genes_on_clust_size(copydiff_island,1)

pairwise_gene_gain = function(pres_abs_file,iso1,iso2,iso1_gff_file,iso2_gff_file,
                              window_size,pw_outdir,countTREE_file,countOutput_file,gainPenalty) {
  #given two isolates, compare gene content differences in roary pres abs table
  #then loads gffs clusters genes into genomic islands given a specified window size
  #output stats about differences: pangenome dist, # of events, mean island size
  comp = paste(as.character(iso1),as.character(iso2),sep='_')
  comparison = c(as.character(iso1),as.character(iso2),comp)
  names(comparison) = c('iso1','iso2','comp')
  print(comparison)
  
  #get stats on pangenome dist between isolates
  pres_abs <- read_csv(pres_abs_file,col_types = cols())
  copynum_table = get_copynum_table(pres_abs,iso1,iso2) 
  genenum_iso1 = sum(copynum_table$copynum_iso1)
  genenum_iso2 = sum(copynum_table$copynum_iso2)
  tab = copynum_table %>% 
    column_to_rownames(var='Gene') %>% 
    dplyr::select(copynum_iso1,copynum_iso2) %>% 
    t()
  (dist_bray <- vegdist(tab, method="bray"))
  genenum_all = c(genenum_iso1,genenum_iso2,round(dist_bray,3))
  names(genenum_all) = c('genenum_iso1','genenum_iso2','dist_bray')
  print(genenum_all)
  
  #subset to genes that have same copy number as mrca
  copynum_mrca_table = get_mrca(countOutput_file,
                                countTREE_file,
                                get_copynum_table(pres_abs,iso1,iso2))
  
  #filter out HGGs where one of the isolate copy number is not equal to mcra
  sum(copynum_mrca_table$mrca)
  subset_mrca_table = copynum_mrca_table %>%
    filter(MRCA_iso1 == 'nochange' | MRCA_iso2 == 'nochange') %>% #genes that match mrca 
    mutate(geneGainLoss = ifelse(MRCA_iso1 == 'gain', 'iso1_gain',
                                 ifelse(MRCA_iso1 == 'loss', 'iso1_loss',
                                        ifelse(MRCA_iso2 == 'gain','iso2_gain',
                                               ifelse(MRCA_iso2 == 'loss', 'iso2_loss','nochange'))))) 
  sum(subset_mrca_table$mrca)
  subset_genenum_iso1 = sum(subset_mrca_table$copynum_iso1)
  subset_genenum_iso2 = sum(subset_mrca_table$copynum_iso2)
  subset_genenum_mrca = sum(subset_mrca_table$mrca)
  tab = subset_mrca_table %>% 
    column_to_rownames(var='Gene') %>% 
    dplyr::select(copynum_iso1,copynum_iso2) %>% 
    t()
  (subset_dist_bray <- vegdist(tab, method="bray"))
  genenum_subset = c(subset_genenum_iso1,subset_genenum_iso2,subset_genenum_mrca,round(subset_dist_bray,3))
  names(genenum_subset) = c('subset_genenum_iso1','subset_genenum_iso2','subset_genenum_mrca','subset_dist_bray')
  print(genenum_subset)
 
 #summarize gene gains and losses since mrca
 summary = subset_mrca_table %>% group_by(geneGainLoss) %>% summarise(count = sum(copy_diff))
 subset_n_gene_diff_bray = sum(subset_mrca_table$copy_diff)
 iso1_gain = summary$count[summary$geneGainLoss == 'iso1_gain']
 iso1_loss = summary$count[summary$geneGainLoss == 'iso1_loss']
 iso2_gain = summary$count[summary$geneGainLoss == 'iso2_gain']
 iso2_loss = summary$count[summary$geneGainLoss == 'iso2_loss']
 #sanity check,
 sum(c(iso1_gain,iso1_loss,iso2_gain,iso2_loss))==subset_n_gene_diff_bray

 geneGainLoss = c(subset_n_gene_diff_bray,iso1_gain,iso1_loss,iso2_gain,iso2_loss)
 names(geneGainLoss)  = c('subset_n_gene_diff_bray','iso1_gain','iso1_loss','iso2_gain','iso2_loss')
 print(geneGainLoss)
 
 #get events
  copydiff_island = get_copy_diff(pres_abs=pres_abs,
                                  subset_mrca_table_output=subset_mrca_table,
                                  iso1=iso1,
                                  iso2=iso2,
                                  iso1_gff_file=iso1_gff_file,
                                  iso2_gff_file=iso2_gff_file,
                                  window_size=window_size)
  #label islands that may represent duplications
  copydiff_island = is_dup(copydiff_island)
  
  #estimate which gene copies where transferred based on island size
  print(paste0(subset_n_gene_diff_bray,' genes differences but ',nrow(copydiff_island),
               ' genes identified in gene families due to multiple gene copies'))
  gene_family_gain = copydiff_island %>% #gene family is present in one isolate not the other
    filter(copynum_iso1 == '0'|copynum_iso2 == '0') 
  gene_copy_gain = copydiff_island %>% #gene family is present in both isolate but copy number differs
    filter(!Gene %in% gene_family_gain$Gene) 
  all_genes = gene_family_gain
  for (diff in unique(gene_copy_gain$copy_diff)) {
    res = pick_genes_on_clust_size(gene_copy_gain,diff)
    all_genes = all_genes %>% add_row(res)
  }
  if (nrow(all_genes) !=subset_n_gene_diff_bray){ #now we've select same number of gene copies to explain extact diff in copy num
    stop("ERROR: gene gain df doesn't match number of unique genes")
  }
  #resize clusters after choosing gene copies
  all_genes = all_genes %>% 
    add_count(name,name = 'cluster_size') %>% 
    arrange(desc(cluster_size)) %>%
    select(-size) %>%
    mutate('iso1'=iso1,'iso2'=iso2,'comp'=paste0(iso1,'_',iso2))
  iso1_G50 = all_genes %>% 
    filter(geneGainLoss == 'iso1_gain') %>%
    slice_max(order_by = cluster_size, prop = .5) %>%
    pull(cluster_size) %>% min() %>% as.numeric()
  iso1_L50 = all_genes %>% 
    filter(geneGainLoss == 'iso1_loss') %>%
    slice_max(order_by = cluster_size, prop = .5) %>%
    pull(cluster_size) %>% min() %>% as.numeric()
  iso2_G50 = all_genes %>% 
    filter(geneGainLoss == 'iso2_gain') %>%
    slice_max(order_by = cluster_size, prop = .5) %>%
    pull(cluster_size) %>% min() %>% as.numeric()
  iso2_L50 = all_genes %>% 
    filter(geneGainLoss == 'iso2_loss') %>%
    slice_max(order_by = cluster_size, prop = .5) %>%
    pull(cluster_size) %>% min() %>% as.numeric()
  P50 = c(iso1_G50,iso1_L50,iso2_G50,iso2_L50)
  names(P50) = c('iso1_G50','iso1_L50','iso2_G50','iso2_L50')
  
  #output islands table 
  island_output = all_genes %>% 
    select(iso1,iso2,comp,name,geneGainLoss,cluster_size,is_dup) %>% 
    distinct() 
  
  #summarize gene gain and loss EVENTS since mrca
  total_events = nrow(island_output)
  summary_events = island_output %>%
    group_by(geneGainLoss) %>% tally()
  
  iso1_gain_events = summary_events$n[summary_events$geneGainLoss == 'iso1_gain']
  iso1_loss_events = summary_events$n[summary_events$geneGainLoss == 'iso1_loss']
  iso2_gain_events = summary_events$n[summary_events$geneGainLoss == 'iso2_gain']
  iso2_loss_events = summary_events$n[summary_events$geneGainLoss == 'iso2_loss']
  #sanity check
  total_events == sum(c(iso1_gain_events,iso1_loss_events,iso2_gain_events,iso2_loss_events))
  events = c(total_events,iso1_gain_events,iso1_loss_events,iso2_gain_events,iso2_loss_events)
  names(events) = c('total_events','iso1_gain_events','iso1_loss_events','iso2_gain_events','iso2_loss_events')
  print(events)
  
  #output island and their sizes
  system(paste0('mkdir -pv ',pw_outdir))
  filename = paste0(comp,'_window',as.character(window_size),'_gp',gainPenalty)
  res=c(comparison,genenum_all,genenum_subset,geneGainLoss,events,P50)
  print(res)
  
  write_tsv(as.data.frame(t(res)),file.path(pw_outdir,paste0(filename,'_summary.txt'))) 
  write_tsv(all_genes,file.path(pw_outdir,paste0(filename,'_gff.txt'))) 
  write_tsv(copynum_mrca_table,file.path(pw_outdir,paste0(filename,'_copynum_mrca.txt')))
  write_tsv(island_output, file.path(pw_outdir,paste0(filename,'_island.txt')))
}


# pres_abs_file = file.path(test_dir,'gene_presence_absence.csv')
# countFAMILYoutput_file = file.path(test_dir,"countOutput_Bacteroides_xylanisolvens_gainpenalty2.FAMILY")
# countTREE_file = file.path(test_dir,"count_Bacteroides_xylanisolvens.tre")
# 
# w1_g2 = pairwise_gene_gain(
#    pres_abs_file = pres_abs_file,
#    iso1='P21.6D',
#    iso2='P17.A3',
#    iso1_gff_file=file.path(test_dir,'P21.6D.gff'),
#    iso2_gff_file=file.path(test_dir,'P17.A3.gff'),
#    window_size = 1,
#    pw_outdir = file.path(test_dir,'output'),
#    countTREE_file=countTREE_file,
#    countOutput_file=countFAMILYoutput_file,
#    gainPenalty=2)
# 
# 
# w5_g2 = pairwise_gene_gain(
#   pres_abs_file = pres_abs_file,
#   iso1='P21.6D',
#   iso2='P17.A3',
#   iso1_gff_file=file.path(test_dir,'P21.6D.gff'),
#   iso2_gff_file=file.path(test_dir,'P17.A3.gff'),
#   window_size = 5,
#   pw_outdir = file.path(test_dir,'output'),
#   countTREE_file=countTREE_file,
#   countOutput_file=countFAMILYoutput_file,
#   gainPenalty=2)
