Figures and Tables
================

``` r
outdir='figures_tables'
system('rm -r figures_tables')
dir.create(outdir)
```

### Copy over Table S1 & S2

``` r
system('cp metadata/TableS1_isolate_genomes.txt figures_tables/TableS1_isolate_genomes.txt')
system('cp metadata/TableS2_ncbi_genomes.txt figures_tables/TableS2_ncbi_genomes.txt')
```

### Figure 1, tree plus cols A and B

``` r
species = 'Bacteroides_xylanisolvens'
dir = paste0('results/pangenome/',species)
tree_file = file.path(dir,paste0(species,'.tre'))
tree = read.tree(tree_file)

#drop outgroups from tree
metadata_file = file.path(dir,'metadata.txt')
metadata = read_tsv(metadata_file,col_types = cols())
outgroup = metadata %>% filter(taxonomy_Species!=species) %>% pull(isolate)
tree = drop.tip(tree,outgroup) 
metadata %>% 
  filter(taxonomy_Species==species) %>%
  group_by(clade,captive_clade) %>% 
  tally()
```

    ## # A tibble: 7 x 3
    ## # Groups:   clade [4]
    ##   clade      captive_clade     n
    ##   <chr>      <chr>         <int>
    ## 1 cladeA     gorilla1          4
    ## 2 cladeA     mixedhost        22
    ## 3 cladeA     unassigned       17
    ## 4 cladeB     gorilla2         13
    ## 5 cladeB     unassigned        9
    ## 6 cladeC     unassigned       60
    ## 7 unassigned unassigned        6

``` r
#add PIC to tree
PI_comps_file = file.path(dir,'gene_gain_loss/PI_comps.txt')
PI_comps = read_tsv(PI_comps_file)
PI_comps_long = PI_comps %>% 
    mutate(compNum = 1:nrow(PI_comps)) %>% 
    tidyr::pivot_longer(cols=c(iso1,iso2), names_to='1vs2') %>%
    dplyr::rename('iso'='value') %>%
    select(iso,everything())
fig1_v1 <- ggtree(tree) %<+% PI_comps_long + 
    geom_tippoint(aes(subset=(label%in%PI_comps_long$iso),
                      color =  factor(compNum)),show.legend = FALSE) +
    geom_nodepoint(aes(subset = suppressWarnings(as.numeric(label)) > 90),size=.75) + #add bootstrap 
    geom_treescale() #add scale

#add clade labels
cladeA_MCRA <- ape::getMRCA(tree,metadata$isolate[metadata$clade=='cladeA'])
cladeB_MCRA <- ape::getMRCA(tree,metadata$isolate[metadata$clade=='cladeB'])
cladeC_MCRA <- ape::getMRCA(tree,metadata$isolate[metadata$clade=='cladeC'])
mixedhostNode <- ape::getMRCA(tree,metadata$isolate[metadata$captive_clade=='mixedhost'])
gorilla1Node <- ape::getMRCA(tree,metadata$isolate[metadata$captive_clade=='gorilla1'])
gorilla2Node <- ape::getMRCA(tree,metadata$isolate[metadata$captive_clade=='gorilla2'])
fig1_v2 = fig1_v1 + #add clade labels
   geom_cladelabel(node=cladeA_MCRA,label="cladeA",offset = .005)+
   geom_cladelabel(node=cladeB_MCRA,label="cladeB",offset = .005)+
   geom_cladelabel(node=cladeC_MCRA,label="cladeC",offset = .005)+
   geom_cladelabel(node=mixedhostNode,label="mixedhost",offset = .01) +
   geom_cladelabel(node=gorilla1Node,label="gorilla1",offset = .01)+ 
   geom_cladelabel(node=gorilla2Node,label="gorilla2",offset = .01) +
   geom_tippoint(aes(subset=(label=='GCA.000210075.1.ASM21007v1')),color='red') 

#column A
HOST_table <- metadata %>% 
  select(isolate,host) %>%
  column_to_rownames(var='isolate') 
offset_val1 = .022

get_color_palette <- function(tips,metadata) {
  metadata <- metadata %>% filter(isolate %in% tips)
  vec <- sort(unique(metadata$host))
  return(recode(vec,
                          'human'='cadetblue4',
                          'rumen' = 'brown4',
                          'bonobo'='red2',
                          'chimpanzee'='orange2',
                          'orangutan'='purple4',
                          'gorilla'='green3',
                          'chicken'='tan',
                          'mouse' = 'yellow2',
                          'pig' = 'pink',
                          'missing'='white'))}

fig1_v3 <- gheatmap(fig1_v2 + ylim(-10,NA),offset=.05,
    HOST_table,width=.1,colnames_angle=90,hjust=1)  + 
    scale_fill_manual(values=get_color_palette(fig1_v2$data$label,metadata))
#column B
predicted.genes = metadata %>% select(isolate,predicted.genes) %>%
  mutate(predicted.genes=as.numeric(predicted.genes)) %>%
  column_to_rownames(var='isolate')
fig1_v3 <- fig1_v3 + new_scale_fill()
fig1_v4 <- gheatmap(fig1_v3,predicted.genes,
                               colnames_angle=90,hjust=1,width=.1,offset = .06) + 
    scale_fill_viridis_c(direction = -1, option="D")
fig1_v4 <- fig1_v4 + new_scale_fill()
fig1_v4
```

![](figures_tables_files/figure-gfm/figure1-1.png)<!-- -->

### Table S3 and sulfatase functional groups convergently enriched

``` r
#Table S3: Functional groups enriched in all three captive ape clade since their mcra with closest human-associated strain 

copyNumTable = read_tsv(file.path(dir,'gene_gain_loss/PIC_gene_gain_loss_summary/Bacteroides_xylanisolvens_window5_gp2_copynum.txt'))
annotation = read_tsv(file= file.path(dir,'Bxy_roary_nosplitparalogs_annotation.txt'),
                      col_types = cols())

copyNumTable_annotation = copyNumTable  %>% #add annotation info to copy number table 
  left_join(select(annotation,Gene,func_group,func_annot,Annotation),by='Gene')
copyNum_funcGroup = copyNumTable_annotation %>% #
  group_by(func_group,func_annot,iso1,iso2) %>% 
  summarise(copyChangeMrca_iso1 = sum(copynum_iso1)-sum(mrca),
            copyChangeMrca_iso2 = sum(copynum_iso2)-sum(mrca)) %>% 
  as.data.frame()
iso1 = copyNum_funcGroup %>% dplyr::select(func_group,func_annot,iso1,copyChangeMrca_iso1) %>% 
  dplyr::rename(isolate = iso1,copyChangeMrca = copyChangeMrca_iso1)
iso2 = copyNum_funcGroup %>% dplyr::select(func_group,func_annot,iso2,copyChangeMrca_iso2) %>%
  dplyr::rename(isolate = iso2,copyChangeMrca = copyChangeMrca_iso2)
TableS3 = rbind(iso1,iso2) %>% 
  as.data.frame() %>% 
  pivot_wider(names_from = isolate, values_from = copyChangeMrca) %>% 
  select(func_group,func_annot,P21.11A,P14.E4,P21.6E) %>% 
  filter(P21.11A>0,P14.E4>0,P21.6E>0) %>%
  as.data.frame() %>% 
  filter(func_group != 'NA_NA_NA') %>% 
  mutate(total = as.numeric(P21.11A+P14.E4+P21.6E))%>% 
  dplyr::rename(gorilla1_P21.11A = P21.11A,
                gorilla2_P21.6E = P21.6E,
                mixedhost_P14.E4 = P14.E4) %>%
  arrange(desc(total))

head(TableS3)
```

    ##                      func_group
    ## 1      COG3119@2|Bacteria_NA_NA
    ## 2 COG3250@2|Bacteria_K01190_GH2
    ## 3     COG0438@2|Bacteria_NA_GT4
    ## 4      COG0582@2|Bacteria_NA_NA
    ## 5        2ZA7H@2|Bacteria_NA_NA
    ## 6  COG1435@2|Bacteria_K21572_NA
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    func_annot
    ## 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               N-acetylgalactosamine-6-O-sulfatase,Endo-4-O-sulfatase,Delta 4,5-hexuronate-2-O-sulfatase,Alkaline phosphatase PafA,hypothetical protein,N-acetylglucosamine-6-O-sulfatase,Arylsulfatase,Ulvan-active sulfatase,Bifunctional sulfatase/alpha-L-rhamnosidase,Choline-sulfatase
    ## 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           Beta-galactosidase,Evolved beta-galactosidase subunit alpha,Beta-galactosidase BoGH2A,Putative beta-glucuronidase
    ## 3 Alpha-D-kanosaminyltransferase,D-inositol-3-phosphate glycosyltransferase,Putative teichuronic acid biosynthesis glycosyltransferase TuaC,N,N'-diacetylbacillosaminyl-diphospho-undecaprenol alpha-1,3-N-acetylgalactosaminyltransferase,GDP-mannose-dependent alpha-(1-6)-phosphatidylinositol monomannoside mannosyltransferase,Glycogen synthase,O-antigen biosynthesis glycosyltransferase WbnH,hypothetical protein,UDP-N-acetylglucosamine--peptide N-acetylglucosaminyltransferase GtfA subunit,Putative glycosyltransferase EpsF,Alpha-maltose-1-phosphate synthase,N-acetylgalactosamine-N,N'-diacetylbacillosaminyl-diphospho-undecaprenol 4-alpha-N-acetylgalactosaminyltransferase,GalNAc-alpha-(1->4)-GalNAc-alpha-(1->3)-diNAcBac-PP-undecaprenol alpha-1,4-N-acetyl-D-galactosaminyltransferase,N-acetyl-alpha-D-glucosaminyl L-malate synthase,Putative glycosyltransferase EpsD,Spore coat protein SA,D-inositol 3-phosphate glycosyltransferase,Alpha-monoglucosyldiacylglycerol synthase,GDP-mannose-dependent alpha-mannosyltransferase
    ## 4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    hypothetical protein,Tyrosine recombinase XerC,Tyrosine recombinase XerD
    ## 5                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        hypothetical protein
    ## 6                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     hypothetical protein,SusD-like protein P25,SusD-like protein P2,SusD-like protein,SusD-like protein P38
    ##   gorilla1_P21.11A mixedhost_P14.E4 gorilla2_P21.6E total
    ## 1                9               24               8    41
    ## 2                5               13               4    22
    ## 3                4               11               2    17
    ## 4                2                9               4    15
    ## 5                4                7               2    13
    ## 6                1                7               5    13

``` r
write_tsv(TableS3,file.path(outdir,'TableS3_genefunc_gained.txt'))

#Break up func group COG3119 to see which sulfatase Annotations are convergently gained
sulfatase_Annotation = copyNumTable_annotation  %>% #
  filter(str_detect(func_group,'COG3119')) %>%
  group_by(Annotation,iso1,iso2) %>% 
  summarise(copyChangeMrca_iso1 = sum(copynum_iso1)-sum(mrca),
            copyChangeMrca_iso2 = sum(copynum_iso2)-sum(mrca)) %>% 
  as.data.frame()
iso1 = sulfatase_Annotation %>% dplyr::select(Annotation,iso1,copyChangeMrca_iso1) %>% 
  dplyr::rename(isolate = iso1,copyChangeMrca = copyChangeMrca_iso1)
iso2 = sulfatase_Annotation %>% dplyr::select(Annotation,iso2,copyChangeMrca_iso2) %>%
  dplyr::rename(isolate = iso2,copyChangeMrca = copyChangeMrca_iso2)
sulfatase_convergent_gain = rbind(iso1,iso2) %>% 
  as.data.frame() %>% 
  pivot_wider(names_from = isolate, values_from = copyChangeMrca) %>% 
  select(Annotation,P21.11A,P14.E4,P21.6E) %>% 
  filter(P21.11A>0,P14.E4>0,P21.6E>0) %>%
  as.data.frame() %>% 
  mutate(total = P21.11A+P14.E4+P21.6E) %>% 
  dplyr::rename(gorilla1_P21.11A = P21.11A,
                gorilla2_P21.6E = P21.6E,
                mixedhost_P14.E4 = P14.E4)
print(sulfatase_convergent_gain)
```

    ##                            Annotation gorilla1_P21.11A mixedhost_P14.E4
    ## 1                       Arylsulfatase                3                5
    ## 2 N-acetylgalactosamine-6-O-sulfatase                2                3
    ## 3   N-acetylglucosamine-6-O-sulfatase                2                4
    ## 4              Ulvan-active sulfatase                6                7
    ##   gorilla2_P21.6E total
    ## 1               2    10
    ## 2               2     7
    ## 3               2     8
    ## 4               5    18

### Figure 1 columns C & D

Visualize function groups and sulfatase annotations convergently
enriched

``` r
#inputs: annotation,gene table
pres_abs <- read_csv(file.path(dir,'roary_nosplitparalogs/gene_presence_absence.csv'),col_types = cols())
isblank <- function(x) {as.numeric(str_count(x, pattern = "_"))}
gene_table <- pres_abs %>% 
      select(Gene,all_of(metadata$isolate))  %>%
      mutate_at(vars(-Gene),isblank) %>%
      column_to_rownames(var = 'Gene') %>% 
      as.matrix()
gene_table[is.na(gene_table)] <- 0
gene_table <- as.data.frame(gene_table) %>% rownames_to_column(var='Gene')

get_func_copynumber_table <- function(list_of_func){
  #given list of function groups returns table with gene copy numbers for each isolate
  func_copynumber_table  <- annotation %>% 
    left_join(gene_table,by='Gene') %>%
    filter(func_group %in% list_of_func) %>%
    pivot_longer(cols=metadata$isolate,
               names_to='isolate',values_to='present') %>%
    select(func_group,isolate,present) %>%
    group_by(func_group,isolate) %>%
    summarise(count = sum(present)) %>%
    pivot_wider(names_from = func_group, values_from = count) %>% 
    column_to_rownames(var='isolate')
    func_copynumber_table[func_copynumber_table==0] <- NA
  return(func_copynumber_table)
}  

mucin <- get_func_copynumber_table(c(
                                     'COG3119@2|Bacteria_K01565_NA', 
                                     'COG1649@2|Bacteria_K05970_NA', 
                                     'COG5434@2|Bacteria_NA_GH110' , 
                                     'COG3250@2|Bacteria_K01190_GH2' 
                                   )) 

mucin <- mucin %>% as.data.frame() %>%
                                    dplyr::rename(
                                           'N-sulfoglucosamine sulfohydrolase(K01565)'='COG3119@2|Bacteria_K01565_NA',
                                           'sialate O-acetylesterase(K05970)'='COG1649@2|Bacteria_K05970_NA',
                                           'alpha-1,3-galactosidase(GH110)'='COG5434@2|Bacteria_NA_GH110',
                                           'beta-galactosidase(K01190)'='COG3250@2|Bacteria_K01190_GH2')

susCD <- get_func_copynumber_table(c(
                                     'COG1435@2|Bacteria_K21572_NA', 
                                     'COG4206@2|Bacteria_K21573_NA'
                                   )) 
susCD<- susCD %>% as.data.frame() %>%
                                    dplyr::rename(
                                      'SusD family protein(K21572)' = 
                                        'COG1435@2|Bacteria_K21572_NA',
                                       'SusC family protein(K21573)' = 
                                        'COG4206@2|Bacteria_K21573_NA')
                                      
carrageenen <- get_func_copynumber_table(c('COG1874@2|Bacteria_NA_GH167', 'NA_NA_GH167',
                                           '33PQM@2|Bacteria_NA_GH150', 'COG5434@2|Bacteria_NA_GH82')) 

carrageenen <- carrageenen %>% as.data.frame() %>%
                                        dplyr::rename('Lambda-carrageenase(GH150)'='33PQM@2|Bacteria_NA_GH150',
                                               'Iota-carrageenase(GH82)'=
                                                 'COG5434@2|Bacteria_NA_GH82',
                                               'GH167a'='COG1874@2|Bacteria_NA_GH167',
                                               'GH167b'='NA_NA_GH167') %>%
                                        mutate(beta_carrageenase=GH167a+GH167b) %>%
                                        dplyr::rename('beta_carrageenase(GH167)'='beta_carrageenase') %>%
                                        select(-GH167a,-GH167b)

#get gene counts for each func group for all isolates 
get_annot_copynumber_table <- function(list_of_annot){
  #given list of function groups returns table with gene copy numbers for each isolate
  annot_copynumber_table  <- annotation %>% 
    left_join(gene_table,by='Gene') %>%
    filter(Annotation %in% list_of_annot) %>%
    pivot_longer(cols=metadata$isolate,
               names_to='isolate',values_to='present') %>%
    select(Annotation,isolate,present) %>%
    group_by(Annotation,isolate) %>%
    summarise(count = sum(present)) %>%
    pivot_wider(names_from = Annotation, values_from = count) %>% 
    column_to_rownames(var='isolate')
    annot_copynumber_table[annot_copynumber_table==0] <- NA
  return(annot_copynumber_table)
}  
sulfatase = get_annot_copynumber_table(c('N-acetylgalactosamine-6-O-sulfatase',
                                         'N-acetylglucosamine-6-O-sulfatase',
                                         'Ulvan-active sulfatase','Arylsulfatase')) 
sulfatase = sulfatase %>% as.data.frame() %>%
                  dplyr::rename('N-acetylgalactosamine-6-O-sulfatase(K01132)'='N-acetylgalactosamine-6-O-sulfatase',
                         'N-acetylglucosamine-6-O-sulfatase(K01137)'='N-acetylglucosamine-6-O-sulfatase'
                         )
sulfa_mucin = cbind(susCD,sulfatase,mucin)   
sulfa_mucin = sulfa_mucin[c(
              'SusC family protein(K21573)',
              'SusD family protein(K21572)',
              "beta-galactosidase(K01190)",
              "N-acetylgalactosamine-6-O-sulfatase(K01132)",
              "N-acetylglucosamine-6-O-sulfatase(K01137)",
              'Arylsulfatase',
              'Ulvan-active sulfatase',
              "sialate O-acetylesterase(K05970)",
              "alpha-1,3-galactosidase(GH110)"
              )]
fig1_v5 <- gheatmap(fig1_v4,sulfa_mucin,colnames_angle=90,hjust=1,offset=.07,width=.9) + 
    scale_fill_viridis_c(direction = -1, option="B")
fig1_v5  <- fig1_v5 + new_scale_fill()
(fig1_v6 <- gheatmap(fig1_v5,carrageenen,colnames_angle=90,hjust=1,offset=.13,width=.3)+ 
    scale_fill_viridis_c(option="D",direction = -1))
```

![](figures_tables_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
ggsave(fig1_v6,file = file.path(outdir,'Figure1_long.pdf'),width = 10, height=25)
ggsave(fig1_v6,file = file.path(outdir,'Figure1_short.pdf'),width = 10, height=6)
```

### Figure S1, S2, S3

``` r
Fig1_sup <- function(species, offset_val1, offset_val2) {
dir = paste0('results/pangenome/',species)

#read in tree
tree_file = file.path(dir,paste0(species,'.tre'))
tree = read.tree(tree_file)

#drop outgroups from tree
metadata_file = file.path(dir,'metadata.txt')
metadata = read_tsv(metadata_file,col_types = cols())
outgroup = metadata %>% filter(taxonomy_Species!=species) %>% pull(isolate)
tree = drop.tip(tree,outgroup) 
print(c(species,length(tree$tip.label)))

#add PIC to tree
PI_comps_file = file.path(dir,'gene_gain_loss/PI_comps.txt')
PI_comps = read_tsv(PI_comps_file)
PI_comps_long = PI_comps %>% 
    mutate(compNum = 1:nrow(PI_comps)) %>% 
    tidyr::pivot_longer(cols=c(iso1,iso2), names_to='1vs2') %>%
    dplyr::rename('iso'='value') %>%
    select(iso,everything())
fig1_v1 <- ggtree(tree) %<+% PI_comps_long + 
    geom_tippoint(aes(subset=(label%in%PI_comps_long$iso),
                      color = factor(compNum)),show.legend = FALSE) +
    geom_nodepoint(aes(subset = suppressWarnings(as.numeric(label)) > 90),size=.75) + #add bootstrap 
    geom_treescale() +#add scale
    theme(legend.position="none")
fig1_v1<- fig1_v1 + new_scale_fill()
#column A
HOST_table <- metadata %>% 
  select(isolate,host) %>%
  column_to_rownames(var='isolate') 

fig1_v2 <- gheatmap(fig1_v1 + ylim(-10,NA),colnames = FALSE,offset = offset_val1,
    HOST_table,width=.1,colnames_angle=90,hjust=1)  + 
    scale_fill_manual(values=get_color_palette(fig1_v1$data$label,metadata))

#column B
predicted.genes = metadata %>% select(isolate,predicted.genes) %>%
  mutate(predicted.genes=as.numeric(predicted.genes)) %>%
  column_to_rownames(var='isolate')
fig1_v2 <- fig1_v2 + new_scale_fill()
(fig1_v3 <- gheatmap(fig1_v2,predicted.genes,
                               colnames_angle=90,hjust=1,width=.1,offset = offset_val2) + 
    scale_fill_viridis_c(direction = -1, option="D"))
return(fig1_v3)
}

(Fig1_Bov = Fig1_sup('Bacteroides_ovatus',0,.012))
```

    ## [1] "Bacteroides_ovatus" "95"

![](figures_tables_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
ggsave(Fig1_Bov,file=file.path(outdir,'FigureS1_BovTree.pdf'),height=5,width=5)
(Fig1_Bfr = Fig1_sup('Bacteroides_fragilis',0,.001))
```

    ## [1] "Bacteroides_fragilis" "182"

![](figures_tables_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
ggsave(Fig1_Bfr,file=file.path(outdir,'FigureS2_BfrTree.pdf'),height=5,width=5)
(Fig1_Bth = Fig1_sup('Bacteroides_thetaiotaomicron',0,.004))
```

    ## [1] "Bacteroides_thetaiotaomicron" "74"

![](figures_tables_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

``` r
ggsave(Fig1_Bth,file=file.path(outdir,'FigureS3_BthTree.pdf'),height=5,width=5)
```

### Figure 2A

``` r
Bt_metadata = 
  read_tsv('results/pangenome/Bacteroides_xylanisolvens/metadata.txt',col_types = cols()) %>%
  sjmisc::add_rows(read_tsv('results/pangenome/Bacteroides_ovatus/metadata.txt',col_types = cols())) %>%
  sjmisc::add_rows(read_tsv('results/pangenome/Bacteroides_fragilis/metadata.txt',col_types = cols())) %>%
  sjmisc::add_rows(read_tsv('results/pangenome/Bacteroides_thetaiotaomicron/metadata.txt',col_types = cols())) %>%
  filter(taxonomy_Species != 'Bacteroides_fragilis_A')

Bt_rep = 
  read_tsv('results/pangenome/Bacteroides_xylanisolvens/gene_gain_loss/pw_50strain.txt',col_types = cols()) %>%
  sjmisc::add_rows(read_tsv('results/pangenome/Bacteroides_ovatus/gene_gain_loss/pw_50strain.txt',col_types = cols())) %>%
  sjmisc::add_rows(read_tsv('results/pangenome/Bacteroides_fragilis/gene_gain_loss/pw_50strain.txt',col_types = cols())) %>%
  sjmisc::add_rows(read_tsv('results/pangenome/Bacteroides_thetaiotaomicron/gene_gain_loss/pw_50strain.txt',col_types = cols())) 

summary(kwAllPairsDunnTest(predicted.genes~ as.factor(taxonomy_Species),
         data=Bt_metadata,
         method="fdr")) 
```

    ##                                                               z value Pr(>|z|)
    ## Bacteroides_ovatus - Bacteroides_fragilis == 0                 13.780  < 2e-16
    ## Bacteroides_thetaiotaomicron - Bacteroides_fragilis == 0        9.827  < 2e-16
    ## Bacteroides_xylanisolvens - Bacteroides_fragilis == 0          12.551  < 2e-16
    ## Bacteroides_thetaiotaomicron - Bacteroides_ovatus == 0          2.414 0.047376
    ## Bacteroides_xylanisolvens - Bacteroides_ovatus == 0             2.320 0.047376
    ## Bacteroides_xylanisolvens - Bacteroides_thetaiotaomicron == 0   0.449 0.653258
    ##                                                                  
    ## Bacteroides_ovatus - Bacteroides_fragilis == 0                ***
    ## Bacteroides_thetaiotaomicron - Bacteroides_fragilis == 0      ***
    ## Bacteroides_xylanisolvens - Bacteroides_fragilis == 0         ***
    ## Bacteroides_thetaiotaomicron - Bacteroides_ovatus == 0          *
    ## Bacteroides_xylanisolvens - Bacteroides_ovatus == 0             *
    ## Bacteroides_xylanisolvens - Bacteroides_thetaiotaomicron == 0

``` r
PIC = function(species) {
  dir = paste0('results/pangenome/',species)
  tree_file = file.path(dir,paste0(species,'.tre'))
  tree = read.tree(tree_file)
  metadata = read_tsv(file.path(dir,'metadata.txt'),col_types = cols()) 
  outgroup = metadata %>% filter(taxonomy_Species!=species) %>% pull(isolate)
  tree = drop.tip(tree,outgroup)
  metadata = metadata %>% filter(!isolate %in% outgroup)
  pred.genes <- metadata$predicted.genes
  names(pred.genes) <- metadata$isolate
  
  pic = pic(pred.genes,tree,var.contrasts=TRUE) %>% 
    as.data.frame() %>% 
    mutate(species = species)
  return(pic)
}
Bxy_pic = PIC('Bacteroides_xylanisolvens')
Bfr_pic = PIC('Bacteroides_fragilis')
Bov_pic = PIC('Bacteroides_ovatus')
Bth_pic = PIC('Bacteroides_thetaiotaomicron')
df = rbind(Bxy_pic,Bfr_pic,Bov_pic,Bth_pic) %>% as.data.frame()
df %>% ggplot(aes(x=species,y=variance)) + geom_boxplot()
```

![](figures_tables_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
fit <- aov(df$variance ~ df$species)
shapiro.test(residuals(fit))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(fit)
    ## W = 0.7566, p-value < 2.2e-16

``` r
kruskal.test(variance ~ species,data=df)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  variance by species
    ## Kruskal-Wallis chi-squared = 37.727, df = 3, p-value = 3.229e-08

``` r
summary(kwAllPairsDunnTest(variance ~ as.factor(species),
         data=df,
         method="fdr"))
```

    ##                                                               z value
    ## Bacteroides_ovatus - Bacteroides_fragilis == 0                  2.704
    ## Bacteroides_thetaiotaomicron - Bacteroides_fragilis == 0        4.304
    ## Bacteroides_xylanisolvens - Bacteroides_fragilis == 0           1.808
    ## Bacteroides_thetaiotaomicron - Bacteroides_ovatus == 0          1.622
    ## Bacteroides_xylanisolvens - Bacteroides_ovatus == 0             4.074
    ## Bacteroides_xylanisolvens - Bacteroides_thetaiotaomicron == 0   5.501
    ##                                                                 Pr(>|z|)    
    ## Bacteroides_ovatus - Bacteroides_fragilis == 0                0.02057069   *
    ## Bacteroides_thetaiotaomicron - Bacteroides_fragilis == 0      8.3734e-05 ***
    ## Bacteroides_xylanisolvens - Bacteroides_fragilis == 0         0.14129746    
    ## Bacteroides_thetaiotaomicron - Bacteroides_ovatus == 0        0.14129746    
    ## Bacteroides_xylanisolvens - Bacteroides_ovatus == 0           0.00018496 ***
    ## Bacteroides_xylanisolvens - Bacteroides_thetaiotaomicron == 0 2.2601e-07 ***

``` r
Bt_rep_metadata = Bt_metadata %>% 
  filter(isolate %in% c(Bt_rep$iso1,Bt_rep$iso2))

colors = recode(sort(unique(Bt_rep_metadata$host)),
                          'human'='cadetblue4',
                          'rumen' = 'brown4',
                          'bonobo'='red2',
                          'chimpanzee'='orange2',
                          'orangutan'='purple4',
                          'gorilla'='green3',
                          'chicken'='tan',
                          'mouse' = 'yellow2',
                          'pig' = 'pink',
                          'missing'='white')

Figure2A = Bt_metadata %>% 
  filter(isolate %in% c(Bt_rep$iso1,Bt_rep$iso2)) %>% 
  ggplot(aes(x =taxonomy_Species,y=predicted.genes,fill=host)) +
  geom_dotplot(binaxis = "y", stackdir = "centerwhole") +
  theme_bw() +
  scale_fill_manual(values=colors)

Bt_rep = Bt_rep %>% group_by(taxonomy_Species.iso1) %>% 
  mutate(norm_tree_dist = (tree_dist-mean(tree_dist)) /sd(tree_dist)) %>%
  filter(norm_tree_dist<1,norm_tree_dist>-1)

#compare models with quadratic and linear relationship with tree dist
(m1 = summary(lm_robust(bray_curtis ~ I(norm_tree_dist^2) + norm_tree_dist + taxonomy_Species.iso1 + 
                    norm_tree_dist * taxonomy_Species.iso1, data = Bt_rep)))
```

    ## 
    ## Call:
    ## lm_robust(formula = bray_curtis ~ I(norm_tree_dist^2) + norm_tree_dist + 
    ##     taxonomy_Species.iso1 + norm_tree_dist * taxonomy_Species.iso1, 
    ##     data = Bt_rep)
    ## 
    ## Standard error type:  HC2 
    ## 
    ## Coefficients:
    ##                                                                   Estimate
    ## (Intercept)                                                       0.212668
    ## I(norm_tree_dist^2)                                              -0.001132
    ## norm_tree_dist                                                    0.013265
    ## taxonomy_Species.iso1Bacteroides_ovatus                           0.082961
    ## taxonomy_Species.iso1Bacteroides_thetaiotaomicron                 0.045791
    ## taxonomy_Species.iso1Bacteroides_xylanisolvens                    0.086627
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_ovatus            0.019303
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_thetaiotaomicron  0.007200
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_xylanisolvens     0.034454
    ##                                                                  Std. Error
    ## (Intercept)                                                        0.001342
    ## I(norm_tree_dist^2)                                                0.002063
    ## norm_tree_dist                                                     0.001847
    ## taxonomy_Species.iso1Bacteroides_ovatus                            0.001553
    ## taxonomy_Species.iso1Bacteroides_thetaiotaomicron                  0.001625
    ## taxonomy_Species.iso1Bacteroides_xylanisolvens                     0.001668
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_ovatus             0.002806
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_thetaiotaomicron   0.002806
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_xylanisolvens      0.002885
    ##                                                                   t value
    ## (Intercept)                                                      158.4854
    ## I(norm_tree_dist^2)                                               -0.5488
    ## norm_tree_dist                                                     7.1808
    ## taxonomy_Species.iso1Bacteroides_ovatus                           53.4177
    ## taxonomy_Species.iso1Bacteroides_thetaiotaomicron                 28.1830
    ## taxonomy_Species.iso1Bacteroides_xylanisolvens                    51.9269
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_ovatus             6.8802
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_thetaiotaomicron   2.5657
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_xylanisolvens     11.9412
    ##                                                                    Pr(>|t|)
    ## (Intercept)                                                       0.000e+00
    ## I(norm_tree_dist^2)                                               5.832e-01
    ## norm_tree_dist                                                    8.438e-13
    ## taxonomy_Species.iso1Bacteroides_ovatus                           0.000e+00
    ## taxonomy_Species.iso1Bacteroides_thetaiotaomicron                1.401e-157
    ## taxonomy_Species.iso1Bacteroides_xylanisolvens                    0.000e+00
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_ovatus            7.060e-12
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_thetaiotaomicron  1.034e-02
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_xylanisolvens     3.055e-32
    ##                                                                   CI Lower
    ## (Intercept)                                                       0.210037
    ## I(norm_tree_dist^2)                                              -0.005178
    ## norm_tree_dist                                                    0.009643
    ## taxonomy_Species.iso1Bacteroides_ovatus                           0.079916
    ## taxonomy_Species.iso1Bacteroides_thetaiotaomicron                 0.042606
    ## taxonomy_Species.iso1Bacteroides_xylanisolvens                    0.083356
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_ovatus            0.013802
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_thetaiotaomicron  0.001698
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_xylanisolvens     0.028797
    ##                                                                  CI Upper   DF
    ## (Intercept)                                                      0.215299 3479
    ## I(norm_tree_dist^2)                                              0.002913 3479
    ## norm_tree_dist                                                   0.016887 3479
    ## taxonomy_Species.iso1Bacteroides_ovatus                          0.086006 3479
    ## taxonomy_Species.iso1Bacteroides_thetaiotaomicron                0.048977 3479
    ## taxonomy_Species.iso1Bacteroides_xylanisolvens                   0.089897 3479
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_ovatus           0.024804 3479
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_thetaiotaomicron 0.012702 3479
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_xylanisolvens    0.040111 3479
    ## 
    ## Multiple R-squared:  0.6695 ,    Adjusted R-squared:  0.6688 
    ## F-statistic: 989.6 on 8 and 3479 DF,  p-value: < 2.2e-16

``` r
#use linear because quad tree_dist term large, 
(m2 = summary(lm_robust(bray_curtis ~  norm_tree_dist + taxonomy_Species.iso1 + norm_tree_dist * taxonomy_Species.iso1, 
                       data = Bt_rep)))
```

    ## 
    ## Call:
    ## lm_robust(formula = bray_curtis ~ norm_tree_dist + taxonomy_Species.iso1 + 
    ##     norm_tree_dist * taxonomy_Species.iso1, data = Bt_rep)
    ## 
    ## Standard error type:  HC2 
    ## 
    ## Coefficients:
    ##                                                                  Estimate
    ## (Intercept)                                                      0.212300
    ## norm_tree_dist                                                   0.013433
    ## taxonomy_Species.iso1Bacteroides_ovatus                          0.083064
    ## taxonomy_Species.iso1Bacteroides_thetaiotaomicron                0.045835
    ## taxonomy_Species.iso1Bacteroides_xylanisolvens                   0.086652
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_ovatus           0.018823
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_thetaiotaomicron 0.007066
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_xylanisolvens    0.034074
    ##                                                                  Std. Error
    ## (Intercept)                                                        0.001137
    ## norm_tree_dist                                                     0.001833
    ## taxonomy_Species.iso1Bacteroides_ovatus                            0.001534
    ## taxonomy_Species.iso1Bacteroides_thetaiotaomicron                  0.001622
    ## taxonomy_Species.iso1Bacteroides_xylanisolvens                     0.001665
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_ovatus             0.002630
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_thetaiotaomicron   0.002798
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_xylanisolvens      0.002754
    ##                                                                  t value
    ## (Intercept)                                                      186.760
    ## norm_tree_dist                                                     7.329
    ## taxonomy_Species.iso1Bacteroides_ovatus                           54.146
    ## taxonomy_Species.iso1Bacteroides_thetaiotaomicron                 28.257
    ## taxonomy_Species.iso1Bacteroides_xylanisolvens                    52.042
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_ovatus             7.156
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_thetaiotaomicron   2.526
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_xylanisolvens     12.373
    ##                                                                    Pr(>|t|)
    ## (Intercept)                                                       0.000e+00
    ## norm_tree_dist                                                    2.872e-13
    ## taxonomy_Species.iso1Bacteroides_ovatus                           0.000e+00
    ## taxonomy_Species.iso1Bacteroides_thetaiotaomicron                2.517e-158
    ## taxonomy_Species.iso1Bacteroides_xylanisolvens                    0.000e+00
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_ovatus            1.005e-12
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_thetaiotaomicron  1.158e-02
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_xylanisolvens     1.927e-34
    ##                                                                  CI Lower
    ## (Intercept)                                                      0.210071
    ## norm_tree_dist                                                   0.009839
    ## taxonomy_Species.iso1Bacteroides_ovatus                          0.080056
    ## taxonomy_Species.iso1Bacteroides_thetaiotaomicron                0.042654
    ## taxonomy_Species.iso1Bacteroides_xylanisolvens                   0.083387
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_ovatus           0.013666
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_thetaiotaomicron 0.001581
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_xylanisolvens    0.028674
    ##                                                                  CI Upper   DF
    ## (Intercept)                                                       0.21453 3480
    ## norm_tree_dist                                                    0.01703 3480
    ## taxonomy_Species.iso1Bacteroides_ovatus                           0.08607 3480
    ## taxonomy_Species.iso1Bacteroides_thetaiotaomicron                 0.04901 3480
    ## taxonomy_Species.iso1Bacteroides_xylanisolvens                    0.08992 3480
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_ovatus            0.02398 3480
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_thetaiotaomicron  0.01255 3480
    ## norm_tree_dist:taxonomy_Species.iso1Bacteroides_xylanisolvens     0.03947 3480
    ## 
    ## Multiple R-squared:  0.6695 ,    Adjusted R-squared:  0.6688 
    ## F-statistic:  1129 on 7 and 3480 DF,  p-value: < 2.2e-16

``` r
coef = as.data.frame(m2$coefficients) %>% rownames_to_column(var = 'coefficient')
write_tsv(coef,file=file.path(outdir,'TableS4_model_coef.txt'))

(Figure2B = Bt_rep  %>% 
  ggplot() + 
  aes(x=norm_tree_dist,y=bray_curtis) +
  geom_point(aes(color=taxonomy_Species.iso1),alpha=1,size=2) +
  theme_bw() +
  scale_colour_manual(values=c('red3','orange','green4','blue4')) +
  stat_smooth(data = filter(Bt_rep ,taxonomy_Species.iso1 == 'Bacteroides_fragilis'),
             method = "lm",color='red3', formula = y ~ x, size = 1, se=TRUE)+
  stat_smooth(data = filter(Bt_rep ,taxonomy_Species.iso1 == 'Bacteroides_ovatus'),
             method = "lm",color='orange', formula = y ~ x, size = 1, se=TRUE)+
  stat_smooth(data = filter(Bt_rep ,taxonomy_Species.iso1 == 'Bacteroides_thetaiotaomicron'),
             method = "lm",color='green4', formula = y ~ x, size = 1, se=TRUE) +
  stat_smooth(data = filter(Bt_rep ,taxonomy_Species.iso1 == 'Bacteroides_xylanisolvens'),
              method = "lm",color='blue4', formula = y ~  x, size = 1, se=TRUE)+
  ylab('Pangenome distance (Bray-Curtis)')+
  xlab('Normalized phylogenetic distance'))
```

![](figures_tables_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
Figure2C <- get_legend(Figure2B)

left = Figure2A+theme(axis.title.x= element_blank())+theme(legend.position = 'bottom')
                 
right =  plot_grid(Figure2B+theme(legend.position="none"),Figure2C,ncol=1,rel_heights = c(1,.5))
Figure2 = plot_grid(left,right,ncol=2,rel_widths = c(1,.75))
ggsave(Figure2,file=file.path(outdir,'Figure2_geneContentBt.pdf'),width = 8)
```

### Figure 3

``` r
species = 'Bacteroides_xylanisolvens'
suffix = 'window5_gp2'
dir = file.path('results/pangenome/',species)
summary_filename = paste0('gene_gain_loss/PIC_gene_gain_loss_summary/',species,'_',suffix,'_summary.txt')
summary_df = read_tsv(file.path(dir,summary_filename))
island_filename = paste0('gene_gain_loss/PIC_gene_gain_loss_summary/',species,'_',suffix,'_island.txt')
island_df = read_tsv(file.path(dir,island_filename))

#Panel A: tree 
tree_file = file.path(dir,paste0(species,'.tre'))
tree = read.tree(tree_file)

PI_comps_file = file.path(dir,'gene_gain_loss/PI_comps.txt')
PI_comps = read_tsv(PI_comps_file)
PI_comps_long = PI_comps %>% 
    mutate(compNum = 1:nrow(PI_comps)) %>% 
    tidyr::pivot_longer(cols=c(iso1,iso2), names_to='1vs2') %>%
    dplyr::rename('iso'='value') %>%
    select(iso,everything())
tree_PI = drop.tip(tree,setdiff(tree$tip.label,PI_comps_long$iso))  
tree_PI
```

    ## 
    ## Phylogenetic tree with 18 tips and 17 internal nodes.
    ## 
    ## Tip labels:
    ##   GCA.009102105.1.ASM910210v1, GCA.015547545.1.ASM1554754v1, GCA.009102805.1.ASM910280v1, GCA.009101945.1.ASM910194v1, GCA.015551805.1.ASM1555180v1, GCA.000210075.1.ASM21007v1, ...
    ## Node labels:
    ##   60, 61, 100, 100, 100, 100, ...
    ## 
    ## Rooted; includes branch lengths.

``` r
fig3_v1 <- ggtree(tree_PI ) %<+% PI_comps_long + 
    geom_tippoint(aes(subset=(label%in%PI_comps_long$iso),
                      color =  factor(compNum)),show.legend = FALSE) +
    geom_nodepoint(aes(subset = suppressWarnings(as.numeric(label)) > 90),size=.75) + #add bootstrap 
    geom_treescale() #add scale

HOST_table <- metadata %>% 
  select(isolate,host) %>%
  filter(isolate %in% tree_PI$tip.label) %>%
  column_to_rownames(var='isolate') 

fig3_v2 <- gheatmap(fig3_v1 + ylim(-10,NA),colnames = FALSE,offset = 0,
    HOST_table,width=.1,colnames_angle=90,hjust=1)  + 
    scale_fill_manual(values=get_color_palette(fig3_v1$data$label,metadata))


#Panel B: Number of gene gain/loss panel
iso1 = summary_df %>% select(iso1,iso1_gain_events,iso1_gain,iso1_loss_events,iso1_loss)
colnames(iso1) = c('isolate','gain_events','gain_numGenes','loss_events','loss_numGenes')
iso2 = summary_df %>% select(iso2,iso2_gain_events,iso2_gain,iso2_loss_events,iso2_loss)
colnames(iso2) = c('isolate','gain_events','gain_numGenes','loss_events','loss_numGenes')
df = rbind(iso1,iso2) %>% 
  as.data.frame() %>% 
  mutate(diff_gain_numGenes = gain_numGenes-gain_events,
         diff_loss_numGenes = loss_numGenes-loss_events) %>%
  select(isolate,gain_events,diff_gain_numGenes,loss_events,diff_loss_numGenes) %>%
  pivot_longer(cols=c(gain_events,diff_gain_numGenes,loss_events,diff_loss_numGenes),
               names_to='cat',values_to='count') %>%
  mutate(count=as.numeric(count),cat=as.factor(cat),
         count = ifelse(cat %in% c('loss_events','diff_loss_numGenes'),-count,count)) %>%
  as.data.frame()
fig3_v3 <- facet_plot(fig3_v2, panel = 'Gene Loss/Gain', data = df, 
                geom = geom_barh, 
                mapping = aes(x = count, fill = as.factor(cat)), 
                stat='identity') +
  scale_fill_manual(values=c("orange2", 'pink','lightblue','firebrick',"green3", "cadetblue4",'blue4',"brown4"))

get_color_palette(fig3_v1$data$label,metadata)
```

    ## [1] "orange2"    "green3"     "cadetblue4" "brown4"

``` r
#Panel C: G50/L50
iso1 = select(summary_df,iso1,iso1_G50,iso1_L50)
colnames(iso1) <- c('isolate','G50','L50')
iso2 = select(summary_df,iso2,iso2_G50,iso2_L50)
colnames(iso2) <- c('isolate','G50','L50')
G50L50 = rbind(iso1,iso2) %>% 
  pivot_longer(cols = c('G50','L50'),values_to = 'count') %>% 
  as.data.frame()
fig3_v4 <- facet_plot(fig3_v3, panel = 'G50/L50', data = G50L50 , 
                geom = geom_point, 
                mapping = aes(x = count, color = as.factor(name)), 
                stat='identity')
fig3_v5 = fig3_v4  + 
  theme_bw()  +
  theme(legend.position="none")

#Panels D & E event size distribution
island_df2 = island_df  %>% 
    separate(col=geneGainLoss,into=c('iso','gainloss'),remove=F,sep='_') %>%
    mutate(iso_category = 
          dplyr::if_else(iso == 'iso1',paste0(iso1,gainloss),paste0(iso2,gainloss)), 
    size_category=cut(cluster_size, breaks=c(0,5,10,15,20,25,50,250)), 
    size_category = recode(size_category,
                                '(0,5]'='01-5',
                                '(5,10]'='06-10',
                                '(10,15]'='11-15',
                                '(15,20]'='16-20',
                                '(20,25]'='21-25',
                                '(25,50]'='26-50',
                                '(50,250]'='50-250'))
#get total number of event/genes by size cat 
island_df3 = island_df2 %>% 
  group_by(iso_category,gainloss,size_category) %>%
  summarise(number_events = n(),
  number_genes = sum(cluster_size)) 
#convert to proportion for each isolate 
island_df4 = island_df3 %>% group_by(gainloss,iso_category) %>%
         mutate(total_events = sum(number_events),
                prop_events = number_events / total_events,
                total_genes = sum(number_genes),
                prop_genes = number_genes / total_genes) 
#pivot magic, some isolates don't have events of a size category, need to put in a 0.
island_df5 = island_df4 %>% 
  select(iso_category,size_category,prop_events,prop_genes) %>%
  pivot_longer(cols = c(prop_events,prop_genes), names_to = 'category', values_to = 'proportion') %>%
  pivot_wider(names_from = size_category,values_from = proportion,values_fill = 0 ) %>%
  pivot_longer(cols = unique(island_df4$size_category), names_to = 'size_category', values_to = 'proportion')
island_df5$gainloss <- factor(island_df5$gainloss, levels = c('loss','gain'))

#Panel D: proportion of gain/loss EVENTS by event size
gain_loss_events_prop = island_df5  %>% 
      filter(category=='prop_events') %>% 
      ggplot(aes(x=size_category,y=proportion,fill=gainloss)) +
      geom_boxplot(position=position_dodge(width = 1)) + 
      theme_bw() +
      xlab('Event size')+
      scale_y_continuous(
        name = "Proportion of events",
        ) + 
      theme(
        legend.position = c(.95, .95), 
        legend.justification = c(.95, .95)
        ) +
      scale_fill_manual(values = c('blue4','firebrick'))
#Panel E: proportion of GENES gained/lost by event size
gain_loss_gene_prop = island_df5 %>% 
      filter(category=='prop_genes') %>% 
      ggplot(aes(x=size_category,y=proportion,fill=gainloss)) +
      geom_boxplot(position=position_dodge(width = .9)) + 
      theme_bw() +
      xlab('Event size')+
      scale_y_continuous(
        name = "Proportion of genes",
        ) + 
      theme(
        legend.position = c(.95, .95), 
        legend.justification = c(.95, .95)
        ) +
      scale_fill_manual(values = c('lightblue','pink'))
bottom = plot_grid(gain_loss_events_prop, gain_loss_gene_prop)

(Figure3 = plot_grid(fig3_v5,bottom,ncol=1,height=10,width=6)) 
```

![](figures_tables_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
ggsave(fig3_v5,file=file.path('Figure3top_geneGainLoss.pdf'),height=4)
ggsave(bottom,file=file.path('Figure3bottom_geneGainLoss.pdf'),height=4)
```

### Duplications

``` r
dup = island_df2  %>% 
    group_by(iso_category,gainloss,size_category,is_dup) %>%
    summarise(number_events = n(),
              number_genes = sum(cluster_size)) 

#what fraction of gain events are duplications?
dup %>% filter(gainloss == 'gain') %>%
  group_by(is_dup) %>% 
  summarize(number_events = sum(number_events)) %>%
  mutate(proportion = number_events/sum(number_events))
```

    ## # A tibble: 2 x 3
    ##   is_dup number_events proportion
    ##    <dbl>         <int>      <dbl>
    ## 1      0          1585      0.860
    ## 2      1           259      0.140

### G50/L50 - sliding window 1 vs sl

``` r
species = 'Bacteroides_xylanisolvens'
suffix = 'window5_gp2'
dir = file.path('results/pangenome/',species)
summary_filename = paste0('gene_gain_loss/PIC_gene_gain_loss_summary/',species,'_',suffix,'_summary.txt')
window5_summary_df = read_tsv(file.path(dir,summary_filename)) %>% mutate(window = 'window5')
suffix = 'window1_gp2'
summary_filename = paste0('gene_gain_loss/PIC_gene_gain_loss_summary/',species,'_',suffix,'_summary.txt')
window1_summary_df = read_tsv(file.path(dir,summary_filename)) %>% mutate(window = 'window1')
summary_df = window5_summary_df %>% sjmisc::add_rows(window1_summary_df)

iso1 = select(summary_df,window,iso1,iso1_G50,iso1_L50)
colnames(iso1) <- c('window','isolate','G50','L50')
iso2 = select(summary_df,window,iso2,iso2_G50,iso2_L50)
colnames(iso2) <- c('window','isolate','G50','L50')
G50L50 = rbind(iso1,iso2) %>% 
  pivot_longer(cols = c('G50','L50'),values_to = 'count') %>% 
  mutate(name = as.factor(name),
         count = as.numeric(count)) %>%
  as.data.frame()

(G50L50_window1_vs_window5 = G50L50 %>% ggplot(aes(x=name,y=count,fill=window)) + 
  geom_boxplot() +
  theme_bw())
```

![](figures_tables_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
ggsave(G50L50_window1_vs_window5,file=file.path(outdir,'FigureSX_G50L50_window1_vs_window5.pdf'))

suffix = 'window1_gp2'
island_filename = paste0('gene_gain_loss/PIC_gene_gain_loss_summary/',species,'_',suffix,'_island.txt')
island_df = read_tsv(file.path(dir,island_filename))
island_df2 = island_df  %>% 
    separate(col=geneGainLoss,into=c('iso','gainloss'),remove=F,sep='_') %>%
    mutate(iso_category = 
          dplyr::if_else(iso == 'iso1',paste0(iso1,gainloss),paste0(iso2,gainloss)), 
    size_category=cut(cluster_size, breaks=c(0,5,10,15,20,25,50,250)), 
    size_category = recode(size_category,
                                '(0,5]'='01-5',
                                '(5,10]'='06-10',
                                '(10,15]'='11-15',
                                '(15,20]'='16-20',
                                '(20,25]'='21-25',
                                '(25,50]'='26-50',
                                '(50,250]'='50-250'))
#get total number of event/genes by size cat 
island_df3 = island_df2 %>% 
  group_by(iso_category,gainloss,size_category) %>%
  summarise(number_events = n(),
  number_genes = sum(cluster_size)) 
#convert to proportion for each isolate 
island_df4 = island_df3 %>% group_by(gainloss,iso_category) %>%
         mutate(total_events = sum(number_events),
                prop_events = number_events / total_events,
                total_genes = sum(number_genes),
                prop_genes = number_genes / total_genes) 
#pivot magic, some isolates don't have events of a size category, need to put in a 0.
island_df5 = island_df4 %>% 
  select(iso_category,size_category,prop_events,prop_genes) %>%
  pivot_longer(cols = c(prop_events,prop_genes), names_to = 'category', values_to = 'proportion') %>%
  pivot_wider(names_from = size_category,values_from = proportion,values_fill = 0 ) %>%
  pivot_longer(cols = unique(island_df4$size_category), names_to = 'size_category', values_to = 'proportion')
island_df5$gainloss <- factor(island_df5$gainloss, levels = c('loss','gain'))

#Panel E: proportion of GENES gained/lost by event size
gain_loss_gene_prop_window1 = island_df5 %>% 
      filter(category=='prop_genes') %>% 
      ggplot(aes(x=size_category,y=proportion,fill=gainloss)) +
      geom_boxplot(position=position_dodge(width = .9)) + 
      theme_bw() +
      xlab('Event size')+
      scale_y_continuous(
        name = "Proportion of genes",
        ) + 
      theme(
        legend.position = c(.95, .95), 
        legend.justification = c(.95, .95)
        ) +
      scale_fill_manual(values = c('lightblue','pink'))

ggsave(gain_loss_gene_prop_window1,file=file.path(outdir,'FigureSX_gain_loss_gene_prop_window1.pdf'))
```

``` r
species = 'Bacteroides_xylanisolvens'
suffix = 'window5_gp2'
dir = file.path('results/pangenome/',species)
summary_filename = paste0('gene_gain_loss/PIC_gene_gain_loss_summary/',species,'_',suffix,'_summary.txt')
Bxy_summary_df = read_tsv(file.path(dir,summary_filename)) %>% mutate(taxonomy_Species = species)

species = 'Bacteroides_ovatus'
dir = file.path('results/pangenome/',species)
summary_filename = paste0('gene_gain_loss/PIC_gene_gain_loss_summary/',species,'_',suffix,'_summary.txt')
Bov_summary_df = read_tsv(file.path(dir,summary_filename)) %>% mutate(taxonomy_Species = species)

species = 'Bacteroides_fragilis'
dir = file.path('results/pangenome/',species)
summary_filename = paste0('gene_gain_loss/PIC_gene_gain_loss_summary/',species,'_',suffix,'_summary.txt')
Bfr_summary_df = read_tsv(file.path(dir,summary_filename)) %>% mutate(taxonomy_Species = species)

species = 'Bacteroides_thetaiotaomicron'
dir = file.path('results/pangenome/',species)
summary_filename = paste0('gene_gain_loss/PIC_gene_gain_loss_summary/',
                          species,'_',suffix,'_summary.txt')
Bth_summary_df = read_tsv(file.path(dir,summary_filename)) %>% mutate(taxonomy_Species = species)
library(sjmisc)
summary_df = Bxy_summary_df %>% 
  sjmisc::add_rows(Bov_summary_df) %>% 
  sjmisc::add_rows(Bfr_summary_df) %>% 
  sjmisc::add_rows(Bth_summary_df)

iso1 = select(summary_df,taxonomy_Species,iso1,iso1_G50,iso1_L50)
colnames(iso1) <- c('taxonomy_Species','isolate','G50','L50')
iso2 = select(summary_df,taxonomy_Species,iso2,iso2_G50,iso2_L50)
colnames(iso2) <- c('taxonomy_Species','isolate','G50','L50')
G50L50 = rbind(iso1,iso2) %>% 
  pivot_longer(cols = c('G50','L50'),values_to = 'count') %>% 
  mutate(name = as.factor(name),
         count = as.numeric(count)) %>%
  as.data.frame()

G50L50_stats = G50L50 %>% 
  group_by(taxonomy_Species) %>% 
  rstatix::kruskal_test(count ~ name) 
G50L50_stats$p_fdr_adj = p.adjust(G50L50_stats$p, method = 'fdr')
print(G50L50_stats)
```

    ## # A tibble: 4 x 8
    ##   taxonomy_Species      .y.       n statistic    df         p method   p_fdr_adj
    ##   <chr>                 <chr> <int>     <dbl> <int>     <dbl> <chr>        <dbl>
    ## 1 Bacteroides_fragilis  count    68     26.8      1   2.31e-7 Kruskal…   9.24e-7
    ## 2 Bacteroides_ovatus    count    28      7.30     1   6.91e-3 Kruskal…   6.91e-3
    ## 3 Bacteroides_thetaiot… count    32     11.2      1   8.36e-4 Kruskal…   1.11e-3
    ## 4 Bacteroides_xylaniso… count    36     16.4      1   5.01e-5 Kruskal…   1.00e-4

``` r
(G50L50_Btspecies = G50L50 %>% ggplot(aes(x=taxonomy_Species,y=count,fill=name)) + 
  geom_boxplot() +
  theme_bw()+
  scale_fill_manual(values = c('firebrick','blue4')))
```

![](figures_tables_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
ggsave(G50L50_Btspecies,file=file.path(outdir,'FigureSX_G50L50_Btspecies.pdf'))
```

## Potential source of gained genes

``` r
#read in gff gain/loss summary 
dir = 'results/pangenome/Bacteroides_xylanisolvens'
gff_gain_loss_file <- file.path(dir,
  'gene_gain_loss/PIC_gene_gain_loss_summary/Bacteroides_xylanisolvens_window5_gp2_gff.txt')

#subset to gained genes
gff = read_tsv(gff_gain_loss_file) %>% 
        filter(diff %in% c('P14.E4','P21.6E','P21.11A'),str_detect(geneGainLoss,'gain')) %>% 
        mutate(isolate = if_else(iso1 %in% c('P14.E4','P21.6E','P21.11A'),iso1,iso2)) %>%
        select(isolate,Gene,Gene.ID)
print(head(gff))
```

    ## # A tibble: 6 x 3
    ##   isolate Gene        Gene.ID       
    ##   <chr>   <chr>       <chr>         
    ## 1 P14.E4  group_12802 BKHAPHKO_01314
    ## 2 P14.E4  group_3031  BKHAPHKO_01379
    ## 3 P14.E4  group_7705  BKHAPHKO_01362
    ## 4 P14.E4  rsrIM       BKHAPHKO_01374
    ## 5 P14.E4  group_17075 BKHAPHKO_01375
    ## 6 P14.E4  group_17076 BKHAPHKO_01376

``` r
#read in diamond blast results for 3 representative captive ape lineages in PICs
diamond_dir = file.path(dir,'diamond')        
diamond = read_tsv(file.path(diamond_dir,'diamond_blast_taxonomy.txt'))

#merge gff with diamond blast
gff_diamond = gff %>% 
  left_join(diamond, by = c(Gene.ID = 'GeneID')) 
gff_diamond$group = replace_na(gff_diamond$group,value='no_hit')

#reclassify infrequent Bacteroides species to Bacteroides other
other_Bt = c('Bacteroides_caccae','Bacteroides_clarus','Bacteroides_congonensis',
             'Bacteroides_faecis','Bacteroides_finegoldii','Bacteroides_fragilis',
             'Bacteroides_intestinalis','Bacteroides_salyersiae','Bacteroides_uniformis',
             'Bacteroides_intestinigallinarum','Bacteroides_zhangwenhongi') 
gff_diamond$group[gff_diamond$group%in%other_Bt] <- 'Bacteroides_other'

#reclassify hits to single Bacteroidales genus to Bacteroidales
Bacteroidales = c('Parabacteroides_merdae','Prevotella_sp.','Prevotella_copri',"Butyricimonas_vaginalis")
gff_diamond$group[gff_diamond$group%in%Bacteroidales] <- 'Bacteroidales'

#reclassify hits to Bacteria
gff_diamond$group[gff_diamond$group == "Campylobacter_jejuni"] <- 'Campylobacterota'
Phyla = c('Proteobacteria','Bacteriophage','Bacteroidetes','Campylobacterota')
gff_diamond$group[gff_diamond$group%in%Phyla] <- 'Bacteria'

#count occurences of gained genes in various phyla
table(gff_diamond$phylalist)
```

    ## 
    ##                                                               Bacteriophage 
    ##                                                                           7 
    ##                                                               Bacteroidetes 
    ##                                                                        1282 
    ##                   Bacteroidetes, Actinobacteria, Firmicutes, Proteobacteria 
    ##                                                                           2 
    ##                               Bacteroidetes, Actinobacteria, Proteobacteria 
    ##                                                                           4 
    ##                                                Bacteroidetes, Bacteriophage 
    ##                                                                         272 
    ##                 Bacteroidetes, Campylobacterota, Proteobacteria, Firmicutes 
    ##                                                                           2 
    ##                                                   Bacteroidetes, Firmicutes 
    ##                                                                           3 
    ##                                    Bacteroidetes, Firmicutes, Bacteriophage 
    ##                                                                           4 
    ## Bacteroidetes, Firmicutes, Campylobacterota, Actinobacteria, Proteobacteria 
    ##                                                                           1 
    ##                                   Bacteroidetes, Firmicutes, Proteobacteria 
    ##                                                                           1 
    ##                    Bacteroidetes, Firmicutes, Proteobacteria, Bacteriophage 
    ##                                                                           1 
    ##                                               Bacteroidetes, Proteobacteria 
    ##                                                                          22 
    ##                                Bacteroidetes, Proteobacteria, Bacteriophage 
    ##                                                                          20 
    ##                                   Bacteroidetes, Proteobacteria, Firmicutes 
    ##                                                                           5 
    ##                                                            Campylobacterota 
    ##                                                                           3 
    ##                                 Firmicutes, Campylobacterota, Bacteroidetes 
    ##                                                                           1 
    ##                                                  Firmicutes, Proteobacteria 
    ##                                                                           1 
    ##                                                              Proteobacteria 
    ##                                                                          28 
    ##                                            Verrucomicrobiota, Bacteroidetes 
    ##                                                                           1

``` r
gff_diamond$phylalist= gff_diamond$phylalist %>% replace_na(value='no_hit')
count_Phyla <- function(a_phylum) {
  return(length(gff_diamond$phylalist[str_detect(gff_diamond$phylalist,a_phylum)==T]))
}

phyla = c('Bacteriophage','Bacteroidetes','Proteobacteria','Campylobacterota',
          'Firmicutes','Actinobacteria','no_hit')
phyla_counts = lapply(phyla,count_Phyla)
names(phyla_counts) <- phyla
print(phyla_counts)
```

    ## $Bacteriophage
    ## [1] 304
    ## 
    ## $Bacteroidetes
    ## [1] 1621
    ## 
    ## $Proteobacteria
    ## [1] 87
    ## 
    ## $Campylobacterota
    ## [1] 7
    ## 
    ## $Firmicutes
    ## [1] 21
    ## 
    ## $Actinobacteria
    ## [1] 7
    ## 
    ## $no_hit
    ## [1] 442

``` r
gff_diamond_grouped = gff_diamond %>% 
  group_by(isolate,group,cat) %>% 
  tally() %>% 
  as.data.frame() 
unique(gff_diamond_grouped$group)
```

    ## [1] "Bacteria"                     "Bacteroidales"               
    ## [3] "Bacteroides"                  "Bacteroides_other"           
    ## [5] "Bacteroides_ovatus"           "Bacteroides_sp."             
    ## [7] "Bacteroides_thetaiotaomicron" "Bacteroides_xylanisolvens"   
    ## [9] "no_hit"

``` r
gff_diamond_grouped$group= factor(gff_diamond_grouped$group, levels =
                                   c("no_hit", "Bacteria", "Bacteroidales",
                                     "Bacteroides","Bacteroides_sp.",
                                      "Bacteroides_other","Bacteroides_ovatus",
                                     "Bacteroides_thetaiotaomicron","Bacteroides_xylanisolvens"))
gff_diamond_grouped$isolate = factor(gff_diamond_grouped$isolate, levels =
                                   c("P14.E4","P21.6E","P21.11A"))
(figure_HGTsource = gff_diamond_grouped %>% 
  ggplot(aes(x=isolate,y=n,fill=group)) + 
  geom_col(position = "fill") + 
  theme_bw()+
  scale_fill_manual(values = brewer.pal(15, "Paired")))
```

![](figures_tables_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
ggsave(figure_HGTsource,file=file.path(outdir,'Figure4_HGTsource.pdf'))
```

``` r
species = 'Bacteroides_xylanisolvens'
suffix = 'window5_gp2'
dir = file.path('results/pangenome/',species)
gff_filename = paste0('gene_gain_loss/PIC_gene_gain_loss_summary/',species,'_',suffix,'_gff.txt')
Bxy_gff_df = read_tsv(file.path(dir,gff_filename)) %>% mutate(taxonomy_Species = species)

mixedhost_gains = Bxy_gff_df %>% filter(iso2 == "P14.E4", geneGainLoss == 'iso2_gain') 
gorilla1_gains = Bxy_gff_df %>% filter(iso1 == "P21.11A", geneGainLoss == 'iso1_gain') 
gorilla2_gains = Bxy_gff_df %>% filter(iso1 == "P21.6E", geneGainLoss == 'iso1_gain') 

mixedhost_HGGs = unique(mixedhost_gains$Gene)
gorilla1_HGGs = unique(gorilla1_gains$Gene)
gorilla2_HGGs = unique(gorilla2_gains$Gene)

HGGs_gained = as.data.frame(table(c(mixedhost_HGGs,gorilla1_HGGs,gorilla2_HGGs)))
colnames(HGGs_gained) <- c('HGG','count') 
HGGs_gained_2lineages = HGGs_gained %>% filter(count>=2)

#gene table
pres_abs <- read_csv(file.path(dir,'roary_nosplitparalogs/gene_presence_absence.csv'),col_types = cols())
isblank <- function(x) {as.numeric(str_count(x, pattern = "_"))}
gene_table <- pres_abs %>% 
      select(Gene,all_of(metadata$isolate))  %>%
      mutate_at(vars(-Gene),isblank) %>%
      column_to_rownames(var = 'Gene') %>% 
      as.matrix()
gene_table[is.na(gene_table)] <- 0
gene_table <- as.data.frame(gene_table) %>% rownames_to_column(var='Gene')

#generate table genes distributed across captive ape clades
captive_clade_table <- gene_table %>% 
  pivot_longer(cols=metadata$isolate,
               names_to='isolate',values_to='present') %>% 
  filter(present>0) %>%
  left_join(select(metadata,isolate,host,human_ape,captive_clade),by='isolate') %>%
  group_by(Gene,captive_clade) %>%
  tally() %>%
  pivot_wider(names_from=captive_clade,values_from=n,values_fill=0) %>%
  filter(unassigned<3)

HGG_table <- gene_table %>% 
  filter(Gene %in% HGGs_gained_2lineages$HGG) %>%
  filter(Gene %in% captive_clade_table$Gene) %>%
  column_to_rownames(var='Gene') %>%
  t()

#generate host table
HOST_table <- metadata %>% 
  dplyr::select(isolate,host_site) %>%
  column_to_rownames(var='isolate') 

get_color_palette <- function(tips) {
  metadata <- metadata %>% filter(isolate %in% tips)
  vec <- sort(unique(metadata$host_site))
  return(recode(vec,
                          'human_USA_Europe'='cadetblue4',
                          'human_China' = 'deepskyblue',
                          'rumen_USA' = 'brown4',
                          'bonobo_Columbus_Zoo'='red2',
                          'chimpanzee_Houston_Zoo'='orange2',
                          'orangutan_Houston_Zoo'='purple4',
                          'gorilla_Columbus_Zoo'='green3',
                          'chicken_Europe'='tan',
                          'missing_siteUnknown'='brown4'))}
#visualize genes
p <- gheatmap(ggtree(tree) + ylim(-10,NA),
              HOST_table,width=.1,colnames_angle=90,hjust=1)  +
              scale_fill_manual(values=get_color_palette(tree$tip.label)) 
p <- p + new_scale_fill()
(p2 <- gheatmap(p,HGG_table,offset=.05,colnames_angle=90,hjust=1) +
                scale_fill_viridis_c(option="C"))
```

![](figures_tables_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
captive_reps_178 = pres_abs %>% 
  filter(Gene %in% colnames(HGG_table)) %>% 
  select(Gene,Annotation,P21.11A,P21.6E,P14.E4) %>%
  separate(P14.E4,sep = '\t',into = c('representative'),remove=F) %>%
  left_join(select(gff_diamond,-isolate,-Gene),by=c('representative' = 'Gene.ID')) 

dir.create(file.path(dir,'output_PUL_fasta'))
saveRDS(captive_reps_178$Gene,file=file.path(dir,'output_PUL_fasta','gene_list.RDS'))


length(captive_reps_178$species_list[str_detect(captive_reps_178$species_list,'xylan')])
```

    ## [1] 0

``` r
write_tsv(captive_reps_178,file=file.path(outdir,'Table_captive_reps_178.txt'))
```

### Figure 5

Compare distribution of sulfatase genes in XB1A model to strains
isolated from captive apes

``` r
sulfa_outdir = file.path(dir,'sulfatase_circle')
dir.create(sulfa_outdir)

#filter to orthologous genes that are within the sulfatlas
sulfatase_genes <- annotation %>% 
  filter(str_detect(eggnog_best_OG_name,'COG3119@2'),
         Gene != 'pafA', #Gene doesn't align well
         Gene!= 'group_4698', #Gene doesn't align well
         No..isolates > 5
         )  #remove single gene group 
### make tree of sulfa genes
set.seed(131)


###subset from all protein seqs to only sulfatase genes in captive ape strains
#all_prot file not uploaded to github showing how sulftase.faa was generated
#all_prot_file = file.path(dir,'all_prot.faa') 
#all_prot = readAAStringSet(all_prot_file)
#names(all_prot) <- sapply(strsplit(names(all_prot)," "), `[`, 1)
#sulfatase_prot = all_prot[names(all_prot) %in% sulfatase_genes$Gene.ID] 
sulfatase_prot_file = file.path(sulfa_outdir,'sulfatase.faa')
#writeXStringSet(sulfatase_prot,sulfatase_prot_file)
#system(paste0('mafft --auto ',sulfatase_prot_file,' > ',sulfatase_prot_file,'.align'))
#system(paste0('fasttree ',sulfatase_prot_file,'.align > ',sulfatase_prot_file,'.tree'))
  
present_in_isolates <- function(list_of_genes,list_of_isolates,name) {
  #given list of genes and isolates, determine the total number of genes present in isolates
  isolate_gene_table <- annotation %>% 
    left_join(gene_table,by='Gene') %>%
    filter(Gene %in% list_of_genes) %>%
    select(Gene,list_of_isolates) %>%
    mutate(total = rowSums(select_(., "-Gene"))) %>%
    select(Gene,total) 
    
  isolate_gene_table = isolate_gene_table %>%
         mutate(binary = ifelse(total>0,name,NA)) %>%
         column_to_rownames('Gene')
  colnames(isolate_gene_table) <-c(paste0(name,'_total'),name)
  return(isolate_gene_table)
}

#determine which of the sulfatase_genes are present in XBA1 and captive ape lineages
mixedhost_cladeA =  metadata %>% filter(captive_clade=='mixedhost') %>% pull(isolate)
mixedhost_cladeA = present_in_isolates(sulfatase_genes$Gene,mixedhost_cladeA,'mixedhost_cladeA')
gorilla1_cladeA =  metadata %>% filter(captive_clade=='gorilla1') %>% pull(isolate)
gorilla1_cladeA = present_in_isolates(sulfatase_genes$Gene,gorilla1_cladeA,'gorilla1_cladeA')
gorilla2_cladeB =  metadata %>% filter(captive_clade=='gorilla2') %>% pull(isolate)
gorilla2_cladeB = present_in_isolates(sulfatase_genes$Gene,gorilla2_cladeB,'gorilla2_cladeB')
human_cladeA =  metadata %>% filter(clade=='cladeA',host=='human') %>% pull(isolate)
human_cladeA = present_in_isolates(sulfatase_genes$Gene,human_cladeA,'human_cladeA')
human_cladeB =  metadata %>% filter(clade=='cladeB',host=='human') %>% pull(isolate)
human_cladeB = present_in_isolates(sulfatase_genes$Gene,human_cladeB,'human_cladeB')
human_cladeC =  metadata %>% filter(clade=='cladeC',host=='human') %>% pull(isolate)
human_cladeC = present_in_isolates(sulfatase_genes$Gene,human_cladeC,'human_cladeC')

pres <- cbind(mixedhost_cladeA,gorilla1_cladeA,gorilla2_cladeB,human_cladeA,human_cladeB,human_cladeC) 
pres = pres %>% select(!contains('total')) 
#remove sulfastase group only present in outgroup Bxy strains
pres = pres[rowSums(is.na(pres)) != ncol(pres), ]
#total number of sulfatase groups
nrow(pres)
```

    ## [1] 91

``` r
#total number of sulfatase groups
core = pres[rowSums(!is.na(pres)) == ncol(pres), ]
nrow(core)
```

    ## [1] 14

``` r
#captive only number of sulfatase groups
captive_only = pres %>% filter(is.na(human_cladeC),is.na(human_cladeB),is.na(human_cladeA))
captive_only %>% group_by_all() %>% tally()
```

    ## # A tibble: 2 x 7
    ## # Groups:   mixedhost_cladeA, gorilla1_cladeA, gorilla2_cladeB, human_cladeA,
    ## #   human_cladeB [2]
    ##   mixedhost_cladeA gorilla1_cladeA gorilla2_cladeB human_cladeA human_cladeB
    ##   <chr>            <chr>           <chr>           <chr>        <chr>       
    ## 1 mixedhost_cladeA gorilla1_cladeA gorilla2_cladeB <NA>         <NA>        
    ## 2 mixedhost_cladeA <NA>            gorilla2_cladeB <NA>         <NA>        
    ## # … with 2 more variables: human_cladeC <chr>, n <int>

``` r
#human only number of sulfatase groups
human_only = pres %>% filter(is.na(mixedhost_cladeA),is.na(gorilla1_cladeA),is.na(gorilla2_cladeB))
human_only %>% group_by_all() %>% tally()
```

    ## # A tibble: 2 x 7
    ## # Groups:   mixedhost_cladeA, gorilla1_cladeA, gorilla2_cladeB, human_cladeA,
    ## #   human_cladeB [2]
    ##   mixedhost_cladeA gorilla1_cladeA gorilla2_cladeB human_cladeA human_cladeB
    ##   <chr>            <chr>           <chr>           <chr>        <chr>       
    ## 1 <NA>             <NA>            <NA>            human_cladeA <NA>        
    ## 2 <NA>             <NA>            <NA>            <NA>         <NA>        
    ## # … with 2 more variables: human_cladeC <chr>, n <int>

``` r
#mixedhost_cladeA number of sulfatase groups
pres %>% filter(!is.na(mixedhost_cladeA)) %>% nrow()
```

    ## [1] 68

``` r
#gorilla1_cladeA number of sulfatase groups
pres %>% filter(!is.na(gorilla1_cladeA)) %>% nrow()
```

    ## [1] 38

``` r
#gorilla2_cladeB number of sulfatase groups
pres %>% filter(!is.na(gorilla2_cladeB)) %>% nrow()
```

    ## [1] 51

``` r
#human_cladeA number of sulfatase groups
pres %>% filter(!is.na(human_cladeA)) %>% nrow()
```

    ## [1] 77

``` r
#human_cladeB number of sulfatase groups
pres %>% filter(!is.na(human_cladeB)) %>% nrow()
```

    ## [1] 38

``` r
#human_cladeC number of sulfatase groups
pres %>% filter(!is.na(human_cladeC)) %>% nrow()
```

    ## [1] 23

``` r
sulfadata  <- pres %>%
  rownames_to_column(var='Gene')%>%
  left_join(annotation,by='Gene') 

#reorder pres abs table
pres = pres[c('human_cladeC','human_cladeB','gorilla2_cladeB',
              'human_cladeA','mixedhost_cladeA','gorilla1_cladeA')]

#label tree tips with sulfatase gene ortholog
sulfa_tree <- read.tree(paste0(sulfatase_prot_file,'.tree')) #read in tree
get_group <- function(prot){sulfadata$Gene[sulfadata$Gene.ID == prot]}
sulfa_tree$tip.label <- as.character(lapply(sulfa_tree$tip.label,get_group))
sulfa_tree <- keep.tip(sulfa_tree,rownames(pres)) #remove tips not in dataframe

#format table for display with ggtree
#annotated function
func <- sulfadata %>% 
  dplyr::select(Gene,Annotation) %>% 
  column_to_rownames(var='Gene') %>% as.matrix() 
#sulfatase subfamily
subfamily <- sulfadata %>% 
  dplyr::select(Gene,sulfatlas_subfamily) %>% 
  filter(Gene %in% sulfa_tree$tip.label) %>% 
  mutate(sulfatase_subfamily = as.character(sulfatlas_subfamily)) %>%
  column_to_rownames(var='Gene') %>% as.matrix() 
colourCount = length(unique(subfamily[,'sulfatlas_subfamily'])) #increase color palette
subfamily_pal = colorRampPalette(brewer.pal(colourCount, "Paired"))(colourCount)

#identify sulfatase genes only present in captive ape clades
coregenes = rownames(pres)[rowSums(!is.na(pres)) == ncol(pres)]

#build sulfatase tree
sulfa_tree_gg = ggtree(sulfa_tree,layout='circular') +
  geom_point2(aes(subset=(label%in%coregenes)), 
              shape=21, size=1, fill='black') + 
  geom_treescale() 

#determine internal tree nodes for sulfatase subfamilies
got_subfam_node <- function(sub){
  subfam = sulfadata$Gene[sulfadata$sulfatlas_subfamily==sub]
  subfam = intersect(subfam,sulfa_tree$tip.label)
  subfam_node = getMRCA(sulfa_tree,subfam)
  return(subfam_node)
}
#got_subfam_node(15)

#number of orthologous gene groups per sulfatase subfamily
table(sulfadata$sulfatlas_subfamily)
```

    ## 
    ##  4  7  8  9 11 14 15 16 19 20 22 24 25 27 28 30 31 46 50 53 55 62 67 
    ##  7  3 12  2  2  3  6  3  1  3  2  1  1  2  2  3  1  1  1  1  1  1  2

``` r
nodes = lapply(unique(sulfadata$sulfatlas_subfamily),got_subfam_node)
names(nodes) = unique(sulfadata$sulfatlas_subfamily)
unlist(nodes) #subfamily & node
```

    ##  15  27   9  20  11  67   4  14   8  22   7  16  30  28 
    ## 153 114 108 149 118 105 143 133 163 160  97 138 112 180

``` r
#singleton subfamilies with 1 rep
singleton_S1 = names(table(sulfadata$sulfatlas_subfamily))[table(sulfadata$sulfatlas_subfamily)==1]
singleton_S1 = sulfadata %>% filter(sulfatlas_subfamily %in% singleton_S1) 
#groups with unassigned subfamilies
unassigned_S1 = sulfadata %>% filter(is.na(sulfatlas_subfamily))
nrow(unassigned_S1)
```

    ## [1] 30

``` r
#add sulfatase subfamilies to tree
sulfa_tree_gg_v2 = sulfa_tree_gg +
  geom_hilight(node=153, fill="darkgreen", alpha=.1) + 
  geom_cladelabel(node=153,label="S1_15") +  
  geom_hilight(node=114, fill="darkgreen", alpha=.1) + 
  geom_cladelabel(node=114,label="S1_27") +  
  geom_hilight(node=108, fill="blue", alpha=.1) +
  geom_cladelabel(node=108,label="S1_9") +  
  geom_hilight(node=149, fill="red", alpha=.1) +
  geom_cladelabel(node=149,label="S1_20") +  
  geom_hilight(node=118, fill="orange", alpha=.1) +
  geom_cladelabel(node=118,label="S1_11") +
  geom_hilight(node=105, fill="yellow", alpha=.1) + 
  geom_cladelabel(node=105,label="S1_67") +  
  geom_hilight(node=143, fill="purple", alpha=.1) +
  geom_cladelabel(node=143,label="S1_4") +  
  geom_hilight(node=133, fill="blue4", alpha=.1) +
  geom_cladelabel(node=133,label="S1_14") +  
  geom_hilight(node=163, fill="pink", alpha=.1) +
  geom_cladelabel(node=163,label="S1_8") +
  geom_hilight(node=160, fill="darkorange", alpha=.1) +
  geom_cladelabel(node=160,label="S1_22") +  
  geom_hilight(node=138, fill="green", alpha=.1) +
  geom_cladelabel(node=138,label="S1_16") +  
  geom_hilight(node=112, fill="red2", alpha=.1) +
  geom_cladelabel(node=112,label="S1_30") +  
  geom_hilight(node=180, fill="yellow4", alpha=.1) +
  geom_cladelabel(node=180,label="S1_28") + 
  geom_treescale() 

sulfa_tree_gg_v3 = sulfa_tree_gg_v2 %<+% sulfadata + xlim(NA, 5) +
    geom_tiplab(aes(label=sulfatlas_subfamily, subset=(label%in%singleton_S1$Gene)), parse=T) +
    geom_tippoint(aes(label=sulfatlas_subfamily, subset=(label%in%singleton_S1$Gene)),fill = 'blue',shape=23,parse=T)
sulfa_tree_gg_v3
```

![](figures_tables_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
color_vec <- recode(sort(colnames(pres)),
      "gorilla1_cladeA"="#1B9E77", #darkgreen
      "mixedhost_cladeA"="#7570B3", #purple
      "gorilla2_cladeB" = "#D95F02", #orange 
      "human_cladeA" = 'dodgerblue4',
      "human_cladeB"= "#E6AB02",
      "human_cladeC"='brown')

p <- gheatmap(sulfa_tree_gg_v3,pres,colnames_angle=90,hjust=1,width=.4)  + 
  scale_fill_manual(values=color_vec) 
p <- p + new_scale_fill()
(p2 <- gheatmap(p,func,colnames_angle=90,hjust=1,offset=1,width=.05) +
    scale_fill_manual(values = brewer.pal(15, "Paired")))
```

![](figures_tables_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->

``` r
ggsave(p2,file=file.path(outdir,'Figure5_sulfatase_circle.pdf'),height=6,width=6)
```

#### Determine AAI within orthologous genes

Shows leg work behind how AAI is calculated, requires all\_prot.faa
generated from processing\_pangenome final table with AAI among
sulfatase genes is uploaded

``` r
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

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   86.27   94.18   97.88   96.75   99.67  100.00

#### blastp of sulfatase genes to sulfatlas

requires download of sulfatlas as show in processing\_pangenome, blastp
results uploads

``` bash
#cd results/pangenome/Bacteroides_xylanisolvens
#blastp -query sulfatase_circle/sulfatase.faa -db sulfatlas/sulfatlas_v1.3.faa -qcov_hsp_perc 80 -outfmt "6 qseqid sseqid #salltitles pident evalue" -out sulfatase_circle/sulfatase_blastp.txt
```

``` r
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
```

    ## # A tibble: 6 x 32
    ##   Gene.ID   Bacteroides_ovat… Bacteroides_sp. Bacteroides_acid… Bacteroides_cac…
    ##   <chr>                 <dbl>           <dbl>             <dbl>            <dbl>
    ## 1 AAFILABK…             100               100              NA               NA  
    ## 2 AAFILABK…             100                NA              NA               NA  
    ## 3 AAFILABK…              99.4              NA              NA               NA  
    ## 4 AAFILABK…              98.4              NA              NA               NA  
    ## 5 AAFILABK…              89.8              NA              NA               NA  
    ## 6 BALBAENG…             100               100              92.5             92.7
    ## # … with 27 more variables: Bacteroides_xylanisolvens <dbl>,
    ## #   Bacteroides_thetaiotaomicron <dbl>, Bacteroides_finegoldii <dbl>,
    ## #   Bacteroides_caecimuris <dbl>, Bacteroides_faecichinchillae <dbl>,
    ## #   Bacteroides_faecis <dbl>, uncultured_Bacteroides <dbl>,
    ## #   Bacteroides_nordii <dbl>, Bacteroides_pyogenes <dbl>,
    ## #   Bacteroides_salyersiae <dbl>, gut_metagenome <dbl>,
    ## #   Parabacteroides_merdae <dbl>, Bacteroides_cellulosilyticus <dbl>,
    ## #   Bacteroides_dorei <dbl>, Bacteroides_fragilis <dbl>,
    ## #   Bacteroides_intestinalis <dbl>, Bacteroides_oleiciplenus <dbl>,
    ## #   Bacteroides_plebeius <dbl>, Bacteroides_stercorirosoris <dbl>,
    ## #   Bacteroides_vulgatus <dbl>, Phocaeicola_vulgatus <dbl>,
    ## #   Bacteroides_clarus <dbl>, Bacteroides_eggerthii <dbl>,
    ## #   Bacteroides_helcogenes <dbl>, Bacteroides_reticulotermitis <dbl>,
    ## #   Bacteroides_stercoris <dbl>, Bacteroides_uniformis <dbl>

``` r
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
```

    ##         Gene        Gene.ID                          Annotation
    ## 1 group_1490 KENJADKH_01220 N-acetylgalactosamine-6-O-sulfatase
    ## 2 group_5562 KENJADKH_04792                  Endo-4-O-sulfatase
    ## 3 group_2339 KENJADKH_02809  Delta 4,5-hexuronate-2-O-sulfatase
    ## 4     atsA_2 KENJADKH_00964 N-acetylgalactosamine-6-O-sulfatase
    ## 5  group_683 KENJADKH_00378 N-acetylgalactosamine-6-O-sulfatase
    ## 6 group_1598 KENJADKH_02788   N-acetylglucosamine-6-O-sulfatase
    ##   sulfatlas_subfamily mixedhost_cladeA gorilla1_cladeA gorilla2_cladeB
    ## 1                  15 mixedhost_cladeA gorilla1_cladeA gorilla2_cladeB
    ## 2                  27 mixedhost_cladeA gorilla1_cladeA gorilla2_cladeB
    ## 3                   9 mixedhost_cladeA gorilla1_cladeA gorilla2_cladeB
    ## 4                  20 mixedhost_cladeA gorilla1_cladeA gorilla2_cladeB
    ## 5                  15 mixedhost_cladeA gorilla1_cladeA gorilla2_cladeB
    ## 6                  11 mixedhost_cladeA gorilla1_cladeA gorilla2_cladeB
    ##   human_cladeA human_cladeB human_cladeC  ave_AAI Bacteroides_xylanisolvens
    ## 1 human_cladeA human_cladeB human_cladeC 91.36276                       100
    ## 2 human_cladeA human_cladeB human_cladeC 96.47488                       100
    ## 3 human_cladeA human_cladeB human_cladeC 94.03305                       100
    ## 4 human_cladeA human_cladeB human_cladeC 91.45549                       100
    ## 5 human_cladeA human_cladeB human_cladeC 86.92707                       100
    ## 6 human_cladeA human_cladeB human_cladeC 93.03145                       100
    ##   Bacteroides_ovatus Bacteroides_thetaiotaomicron Bacteroides_finegoldii
    ## 1             99.022                       89.041                 91.585
    ## 2             97.868                       91.535                 92.636
    ## 3             97.131                       94.467                 93.840
    ## 4             95.745                       87.835                 88.660
    ## 5             99.803                       87.598                 98.819
    ## 6             97.154                       92.220                 93.416
    ##   Bacteroides_caecimuris Bacteroides_faecis gut_metagenome
    ## 1                 96.673             88.650         86.301
    ## 2                 98.643             91.929         85.328
    ## 3                 97.336             92.828             NA
    ## 4                     NA             86.804             NA
    ## 5                 98.622                 NA             NA
    ## 6                     NA             91.841             NA
    ##   Parabacteroides_merdae Bt_sp
    ## 1                     NA   100
    ## 2                     NA   100
    ## 3                     NA   100
    ## 4                     NA   100
    ## 5                     NA   100
    ## 6                     NA   100

``` r
write_tsv(table_sulfadata_captive,file.path(outdir,'TableS5_sulfatase_blastp.txt'))
```

``` r
sessionInfo(package = NULL)
```

    ## R version 3.6.1 (2019-07-05)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /stor/system/opt/R/R-3.6.1/lib/R/lib/libRblas.so
    ## LAPACK: /stor/system/opt/R/R-3.6.1/lib/R/lib/libRlapack.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] sjmisc_2.8.6        PMCMRplus_1.9.0     Biostrings_2.54.0  
    ##  [4] XVector_0.26.0      IRanges_2.20.2      S4Vectors_0.24.4   
    ##  [7] BiocGenerics_0.32.0 RColorBrewer_1.1-2  ggstance_0.3.5     
    ## [10] ggnewscale_0.4.5    vegan_2.5-7         lattice_0.20-41    
    ## [13] permute_0.9-5       harrietr_0.2.3      cowplot_1.1.1      
    ## [16] seqinr_4.2-5        forcats_0.5.1       stringr_1.4.0      
    ## [19] dplyr_1.0.6         purrr_0.3.4         readr_1.4.0        
    ## [22] tidyr_1.1.3         tibble_3.1.1        ggplot2_3.3.3      
    ## [25] tidyverse_1.3.1     estimatr_0.30.2     ggtree_2.0.4       
    ## [28] treeio_1.10.0       ape_5.4-1          
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] TH.data_1.0-10      colorspace_2.0-0    rio_0.5.16         
    ##  [4] ellipsis_0.3.1      sjlabelled_1.1.7    estimability_1.3   
    ##  [7] fs_1.5.0            rstudioapi_0.13     farver_2.0.3       
    ## [10] fansi_0.4.1         mvtnorm_1.1-1       lubridate_1.7.10   
    ## [13] xml2_1.3.2          codetools_0.2-18    splines_3.6.1      
    ## [16] knitr_1.30          SuppDists_1.1-9.5   ade4_1.7-16        
    ## [19] Formula_1.2-4       jsonlite_1.7.2      broom_0.7.6        
    ## [22] cluster_2.1.0       Rmpfr_0.8-3         dbplyr_2.1.1       
    ## [25] BiocManager_1.30.10 compiler_3.6.1      httr_1.4.2         
    ## [28] rvcheck_0.1.8       emmeans_1.5.3       backports_1.2.1    
    ## [31] assertthat_0.2.1    Matrix_1.3-2        lazyeval_0.2.2     
    ## [34] cli_2.5.0           htmltools_0.5.1.1   prettyunits_1.1.1  
    ## [37] tools_3.6.1         gmp_0.6-1           coda_0.19-4        
    ## [40] gtable_0.3.0        glue_1.4.2          Rcpp_1.0.6         
    ## [43] carData_3.0-4       cellranger_1.1.0    vctrs_0.3.6        
    ## [46] nlme_3.1-151        insight_0.13.2      xfun_0.20          
    ## [49] openxlsx_4.2.3      rvest_1.0.0         lifecycle_1.0.0    
    ## [52] rstatix_0.6.0       MASS_7.3-53         zlibbioc_1.32.0    
    ## [55] zoo_1.8-8           scales_1.1.1        hms_1.0.0          
    ## [58] sandwich_3.0-0      curl_4.3            yaml_2.2.1         
    ## [61] memoise_1.1.0       stringi_1.5.3       tidytree_0.3.3     
    ## [64] zip_2.1.1           rlang_0.4.10        pkgconfig_2.0.3    
    ## [67] evaluate_0.14       labeling_0.4.2      tidyselect_1.1.0   
    ## [70] plyr_1.8.6          magrittr_2.0.1      R6_2.5.0           
    ## [73] generics_0.1.0      multcompView_0.1-8  multcomp_1.4-15    
    ## [76] BWStest_0.2.2       DBI_1.1.0           foreign_0.8-76     
    ## [79] pillar_1.6.0        haven_2.3.1         withr_2.3.0        
    ## [82] mgcv_1.8-33         abind_1.4-5         survival_3.2-7     
    ## [85] car_3.0-10          modelr_0.1.8        crayon_1.4.1       
    ## [88] utf8_1.1.4          rmarkdown_2.6       kSamples_1.2-9     
    ## [91] progress_1.2.2      grid_3.6.1          readxl_1.3.1       
    ## [94] data.table_1.13.6   reprex_2.0.0        digest_0.6.27      
    ## [97] xtable_1.8-4        munsell_0.5.0       viridisLite_0.3.0
