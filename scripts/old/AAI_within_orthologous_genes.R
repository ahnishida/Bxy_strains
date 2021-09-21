#### Determine AAI within orthologous genes 
Shows leg work behind how AAI is calculated, requires all_prot.faa generated from processing_pangenome
final table with AAI among sulfatase genes is uploaded 
```{r echo=TRUE, message=FALSE, warning=FALSE}
AAI_outdir = file.path(dir,'gene_AAI')
dir.create(AAI_outdir)

get_ave_AAI <- function(gene){
  #given a orthologous gene extracts faa of all seqs, performs alignment,
  #writes matrix of AAI among seqs
  GeneIDs <- pres_abs %>% filter(Gene==gene) %>% select(metadata$isolate) %>% as.character()
  GeneIDs <- unlist(strsplit(GeneIDs, "[\t]"))
  GeneIDs = GeneIDs[!is.na(GeneIDs)]
  GeneID_seqs = all_prot[names(all_prot) %in% GeneIDs]
  print(paste(gene,length(GeneID_seqs)))
  if (length(GeneID_seqs)>1){
    file = file.path(AAI_outdir,paste0(gene,'.faa'))
    writeXStringSet(GeneID_seqs,file)
    system(paste0('mafft --auto ',file,' > ',file,'.align'))
    prot_align = read.alignment(file=paste0(file,'.align'),format='fasta')
    dist = dist.alignment(prot_align,"identity")
    dist = as.matrix(dist)
    dist[lower.tri(dist,diag = TRUE)] <- NA
    dist = 1-dist
    dist_file = file.path(AAI_outdir,paste0(gene,'_dist.txt'))
    write_tsv(as.data.frame(dist),file=dist_file)
    ave_AAI = mean(dist,na.rm=TRUE)*100
    return(c(gene,ave_AAI))}
  else {
    return(NA)}}
#get_ave_AAI(gene = "group_683") 

get_AAI_from_dist <- function(gene){
  #given gene reads dist matrix and 
  #return average amino acid identity 
  dist_file = file.path(AAI_outdir,paste0(gene,'_dist.txt'))
  dist = as.matrix(read_tsv(file=dist_file),col_types=cols())
  ave_AAI = mean(dist,na.rm=TRUE)*100
  return(data.frame(gene,ave_AAI))
}
#get_AAI_from_dist(gene = "group_683") 


#completed = list.files(AAI_outdir,pattern = '_dist.txt')
#completed = unique(sapply(str_split(completed,'_dist'), `[`, 1))
#sulfatase_to_do = setdiff(sulfatase_genes$Gene,completed)
#mclapply(sulfatase_to_do,get_ave_AAI,mc.cores=20)
#sulfatase_OG_AAI = lapply(sulfatase_genes$Gene,get_AAI_from_dist) %>% bind_rows()
#write_tsv(sulfatase_OG_AAI,file=file.path(sulfa_outdir,'sulfatase_OG_AAI.txt'))
sulfatase_OG_AAI = read_tsv(file=file.path(sulfa_outdir,'sulfatase_OG_AAI.txt'))
summary(sulfatase_OG_AAI$ave_AAI)
```