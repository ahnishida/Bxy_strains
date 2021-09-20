#### Calculate ANI between 50 representative strains for each Bacteroides species

output_genome_list = function(df_pw) {
  #outputs a list of 50 represent genomes in a format that fastANI will take
  species = df_pw$taxonomy_Species.iso1[1]
  dir = paste0('results/pangenome/',species)
  metadata = read_tsv(file.path(dir,'metadata.txt'),col_types = cols())
  rep50 = unique(c(df_pw$iso1,df_pw$iso2))
  rep50_old = lapply(rep50, function(x) {metadata$isolate_old[metadata$isolate==x]})
  rep50_old_files = ifelse(str_detect(rep50_old,'GCA'),
                           paste0('data/ncbi_genomes/isolate_genomes/',rep50_old,'.fna.gz'),
                           paste0('results/processing/scaffold_assemblies/',rep50_old,'.scaffolds.fasta')
  )
  write.table(rep50_old_files,file.path(dir,'gene_gain_loss/rep50_genomes.txt'),
              row.names = F, quote=F, col.names = F)
}

output_genome_list(Bxy_pw)
output_genome_list(Bov_pw)
output_genome_list(Bfr_pw)
output_genome_list(Bth_pw)

#run fastANI on commandline
#species=Bacteroides_fragilis
#fastANI --refList results/pangenome/$species/gene_gain_loss/rep50_genomes.txt --ql results/pangenome/$species/gene_gain_loss/rep50_genomes.txt -o results/pangenome/$species/gene_gain_loss/rep50_genomes_fastANI.txt -t 15

#combine fastANI with pw_df
combine_fastANI_w_braycurtis <- function(species) {
  dir = paste0('results/pangenome/',species)
  pw_df_file = file.path(dir,'gene_gain_loss/pw_50strain.txt')
  fastANI_file = file.path(dir,'gene_gain_loss/rep50_genomes_fastANI.txt')
  metadata = read_tsv(file.path(dir,'metadata.txt'),col_types = cols())
  fastANI = read_tsv(fastANI_file,col_names = F, col_types = cols()) %>%
    separate(X1,'/',into = c(NA,NA,NA,'iso1')) %>%
    separate(X2,'/',into = c(NA,NA,NA,'iso2')) %>%
    mutate(iso1 = str_remove(iso1, '.fna.gz')) %>%
    mutate(iso1 = str_remove(iso1, '.scaffolds.fasta')) %>%
    mutate(iso2 = str_remove(iso2, '.fna.gz')) %>%
    mutate(iso2 = str_remove(iso2, '.scaffolds.fasta')) %>%
    rename('ANI' = 'X3') %>%
    select(iso1,iso2,ANI) %>%
    as.data.frame() 
  fastANI$iso1 = lapply(fastANI$iso1, function(x) {metadata$isolate[metadata$isolate_old==x]})
  fastANI$iso1 = as.character(fastANI$iso1)
  fastANI$iso2 = lapply(fastANI$iso2, function(x) {metadata$isolate[metadata$isolate_old==x]})
  fastANI$iso2 = as.character(fastANI$iso2)
  
  pw_df = read_tsv(pw_df_file,col_types = cols()) %>%
    left_join(as.data.frame(fastANI),by=c('iso1','iso2'))
  pw_df = write_tsv(pw_df, file = pw_df_file)  
}

combine_fastANI_w_braycurtis('Bacteroides_xylanisolvens')
combine_fastANI_w_braycurtis('Bacteroides_ovatus')
combine_fastANI_w_braycurtis('Bacteroides_fragilis')
combine_fastANI_w_braycurtis('Bacteroides_thetaiotaomicron')

