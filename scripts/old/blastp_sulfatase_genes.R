#### blastp of sulfatase genes to sulfatlas
requires download of sulfatlas as show in processing_pangenome, blastp results uploads
```{bash}
#cd results/pangenome/Bacteroides_xylanisolvens
#blastp -query sulfatase_circle/sulfatase.faa -db sulfatlas/sulfatlas_v1.3.faa -qcov_hsp_perc 80 -outfmt "6 qseqid sseqid #salltitles pident evalue" -out sulfatase_circle/sulfatase_blastp.txt
```

```{r echo=TRUE, message=FALSE, warning=FALSE}
#read in blast results
sulfa_blastp = read_tsv(file.path(sulfa_outdir,'sulfatase_blastp.txt'),
                        comment='#',col_names = FALSE, col_types = cols())

colnames(sulfa_blastp) <- c('Gene.ID', 'sseqid','salltitles', 'pident', 'evalue_sulfa')
sulfa_blastp <- sulfa_blastp %>%
  #mutate(Gene.ID=str_sub(Gene.ID, end=-3)) %>% 
  filter(pident > 85) %>% 
  separate(col = 'salltitles',into=c(NA,'id','full_descripton'),extra='drop',sep='[|]') %>%
  separate(col='id',into=c('lcl','sulfatase_family','sulfatase_subfamily'),sep='_') %>%
  separate(col='full_descripton',into=c('description','species'),sep='OS=') %>%
  separate(col='description',into=c('front','description'),sep='BACE |BACOV |BACO1 |_BACSE ')  %>%
  separate(col='species',into=c('species',NA),sep='OX=') %>% 
  separate(col='species',into=c('Genus','species'),sep=' ') %>% 
  mutate(Genus_sp = paste(Genus,species,sep='_')) %>%
  select('Gene.ID', 'sseqid','pident', 'evalue_sulfa',
         'sulfatase_family','sulfatase_subfamily','description','Genus_sp') 

#get best blastp for each Bacteroides species 
sulfa_blastp <- sulfa_blastp %>% 
  select(Gene.ID,Genus_sp,pident) %>% 
  distinct() %>%  
  group_by(Gene.ID,Genus_sp) %>%
  slice_max(order_by=pident,n = 1) %>%
  as.data.frame() %>%
  pivot_wider(names_from=Genus_sp,values_from=pident,values_fill=NA)
head(sulfa_blastp)

#combine with AAI results
table_sulfadata_captive = sulfadata %>% 
  select(Gene,Gene.ID,Annotation,sulfatlas_subfamily,mixedhost_cladeA,gorilla1_cladeA,gorilla2_cladeB,human_cladeA,human_cladeB,human_cladeC)  %>%
  left_join(sulfatase_OG_AAI,by=c('Gene'='gene')) %>% 
  left_join(sulfa_blastp,by=c('Gene.ID')) %>%
  select(Gene:ave_AAI,Bacteroides_xylanisolvens,Bacteroides_ovatus,Bacteroides_thetaiotaomicron,
         Bacteroides_finegoldii,Bacteroides_caecimuris,Bacteroides_faecis,
         gut_metagenome,Parabacteroides_merdae,everything()) 

#get max pident for other Bacteroides species
Bt_sp = table_sulfadata_captive %>%
  select(Bacteroides_sp.:Bacteroides_uniformis) %>% 
  replace(is.na(.), 0)  %>% 
  mutate_all(as.numeric) %>% 
  rowwise() %>% 
  summarize(max = max(c_across(),na.rm=T)) %>%
  pull(max)

table_sulfadata_captive = table_sulfadata_captive %>% 
  select(Gene:ave_AAI,Bacteroides_xylanisolvens,Bacteroides_ovatus,Bacteroides_thetaiotaomicron,
         Bacteroides_finegoldii,Bacteroides_caecimuris,Bacteroides_faecis,
         gut_metagenome,Parabacteroides_merdae)
table_sulfadata_captive$Bt_sp = Bt_sp
head(table_sulfadata_captive)
write_tsv(table_sulfadata_captive,file.path(outdir,'TableS5_sulfatase_blastp.txt'))
```