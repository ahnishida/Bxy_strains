identifying\_source\_genes
================

\-Run diamond blast of faa files of 18 Bxy strains in 9 PICs against nr
database

*./scripts/pangenome\_analyses/source\_genes/diamond\_PIC.sh*

\-Read in results determine appropiate pident and palign cutoff values
to find hits to other bacterial taxa that are comparably to hits to Bxy
strains and subset to hits to those above that threshold

\-Combine with taxonomic table from gtdb to assign higher taxonomic
ranks (Phyla,Order,Family) for hits to determine taxonomic distribution
of high similarity hits

\-Output with table with GeneID and taxonomic distribution to be used by
figures\_tables.Rmd

``` r
diamond_dir = 'results/pangenome/Bacteroides_xylanisolvens/diamond/'
all_files = list.files(file.path(diamond_dir))
(select_files = all_files[str_detect(all_files,'_blast.txt')])
```

    ##  [1] "GCA.000210075.1.ASM21007v1_blast.txt"                             
    ##  [2] "GCA.000273315.1.Bact.xyla.CL03T12C04.V1_blast.txt"                
    ##  [3] "GCA.003458755.1.ASM345875v1_blast.txt"                            
    ##  [4] "GCA.003464445.1.ASM346444v1_blast.txt"                            
    ##  [5] "GCA.003468875.1.ASM346887v1_blast.txt"                            
    ##  [6] "GCA.003474645.1.ASM347464v1_blast.txt"                            
    ##  [7] "GCA.004167295.1.ASM416729v1_blast.txt"                            
    ##  [8] "GCA.009101945.1.ASM910194v1_blast.txt"                            
    ##  [9] "GCA.009102085.1.ASM910208v1_blast.txt"                            
    ## [10] "GCA.009102105.1.ASM910210v1_blast.txt"                            
    ## [11] "GCA.009102805.1.ASM910280v1_blast.txt"                            
    ## [12] "GCA.015551805.1.ASM1555180v1_blast.txt"                           
    ## [13] "GCA.900107825.1.IMG.taxon.2623620516.annotated.assembly_blast.txt"
    ## [14] "GCA.900114865.1.IMG.taxon.2654588180.annotated.assembly_blast.txt"
    ## [15] "P14.E4_blast.txt"                                                 
    ## [16] "P21.11A_blast.txt"                                                
    ## [17] "P21.6E_blast.txt"

``` r
# defining the appropriate pident and palign cutoff
infile = file.path(diamond_dir,select_files[1])
df = read.csv(infile,sep='\t',header=F)
colnames(df) = c('qseqid', 'sseqid', 'pident', 'length', 'qlen','mismatch', 'gapopen',
               'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'salltitles')
df = df %>% mutate(Bxy_hit = if_else(grepl("xylan", salltitles),1,0))
Bxy_df = df %>% filter(Bxy_hit == 1) 
summary(Bxy_df$pident)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   27.70   97.90   99.10   96.88   99.70  100.00

``` r
Bxy_df = Bxy_df %>% filter(pident>.98)
summary(Bxy_df$length/Bxy_df$qlen)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ## 0.08543 1.00000 1.00000 0.97835 1.00000 2.33333

``` r
#apply 98% pident and 98% palign cutoff to hits
read_diamond <- function(infile){
  outfile = str_replace(infile,'_blast.txt','_blast98.txt')
  df = read.csv(infile,sep='\t',header=F)
  colnames(df) = c('qseqid', 'sseqid', 'pident', 'length', 'qlen','mismatch', 'gapopen',
               'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'salltitles')
  df = df %>% mutate(Bxy_hit = ifelse(str_contains('xylanisolvens',salltitles),1,0))
  
  df = df %>% filter(pident>98 & length/qlen > .98) 
  df = df %>% separate_rows(salltitles,sep='<>')
  df = df %>% mutate(taxonomy = stringr::str_extract(string = df$salltitles,
                              pattern = "(?<=\\[).*(?=\\])"),
                     taxonomy = sapply(str_split(taxonomy,'\\['), tail, 1),
                     taxonomy = str_replace(taxonomy,"Candidatus ",""),
                     taxonomy = str_replace(taxonomy,"uncultured ",""),
                     taxonomy = str_replace(taxonomy,"unclassified ",""),
                     palign = length/qlen,
                     Genus = sapply(str_split(taxonomy,'\\ '), head, 1),
                     Species = sapply(str_split(taxonomy,'\\ '), "[", 2),
                     Species = replace_na(Species,value='sp.'),
                     Genus_Species = paste0(Genus,'_',Species)) %>%
                     as.data.frame()
  df = df %>% group_by(qseqid,Genus_Species) %>% 
      slice_min(order_by = evalue, n=1,with_ties = F) 
  write_tsv(df,outfile)
  return(df)
}

#filter blast hits 
#read_diamond(file.path(diamond_dir,'P14.E4_blast.txt'))
#read_diamond(file.path(diamond_dir,'P21.6E_blast.txt'))
#read_diamond(file.path(diamond_dir,'P21.11A_blast.txt'))

#use bac120 taxonomy in GTDB to define higher taxonomic ranks for hits
bac = read_tsv('bin/bac120_taxonomy.tsv',col_names = F,col_types = cols()) %>%
  separate(X2,sep=';',into = c('Domain','Phylum','Class','Order','Family','Genus','Species'))
bac = bac %>% mutate(Phylum = str_replace(Phylum,'p__',''),
                     Order = str_replace(Order,'o__',''), 
                     Family = str_replace(Family,'f__',''), 
                     Genus = str_replace(Genus,'g__',''),
                     Species = str_replace(Species,'s__',''),
                     Species = str_replace(Species,' ','_')) %>%
                mutate(Phylum = recode(Phylum,
                                       Bacteroidota = 'Bacteroidetes',
                                       Firmicutes_A = 'Firmicutes',
                                       Actinobacteriota = 'Actinobacteria'))

#merge results for 3 representatives captive ape lineages
diamond = read_tsv(file.path(diamond_dir,'P14.E4_blast98.txt'),col_types = cols()) %>%
    add_row(read_tsv(file.path(diamond_dir,'P21.6E_blast98.txt'),col_types = cols())) %>%
    add_row(read_tsv(file.path(diamond_dir,'P21.11A_blast98.txt'),col_types = cols())) 

#merge with GTDB to get phylum, order and family designations
diamond_genus = diamond %>% 
  left_join(distinct(select(bac,Phylum,Order,Family,Genus)),by = 'Genus') %>%
  rename('GeneID' = 'qseqid') %>%
  select(GeneID,pident,palign,salltitles,Phylum,Order,Family,Genus,Species,Genus_Species) 

#assign taxonomy info for hist whose entry for Genus is not found in GTDB
diamond_notfound = diamond_genus %>% filter(is.na(Phylum))
sort(table(diamond_notfound$Genus))
```

    ## 
    ##       Dechloromonas    Flavobacteriales            Hallella    Erysipelotrichia 
    ##                   1                   1                   1                   2 
    ##               virus      Avibacteroides   Dysgonamonadaceae   Dysgonomonadaceae 
    ##                   2                   3                   3                   3 
    ##          Firmicutes       Lentisphaeria            Limisoma            Shigella 
    ##                   3                   3                   3                   3 
    ##      Proteobacteria       Rikenellaceae         Caccoplasma    Gallibacteroides 
    ##                   4                   4                   6                   6 
    ##       Solobacterium Alphaproteobacteria          CrAss-like       Clostridiales 
    ##                   6                   7                   8                   9 
    ##          Clostridia  Porphyromonadaceae        Microviridae         Podoviridae 
    ##                  10                  10                  12                  12 
    ##            organism          Tyzzerella          Myoviridae      Muribaculaceae 
    ##                  15                  15                  19                  21 
    ##       Bacteroidetes     Caedimonadaceae         Bacteroidia        Mediterranea 
    ##                  30                  32                  33                  33 
    ##   Paludibacteraceae       Bacteriophage    Oscillospiraceae     Lachnospiraceae 
    ##                  57                  64                  66                  91 
    ##           bacterium            Bacteria        Siphoviridae      Bacteroidaceae 
    ##                 139                 193                 547                 728 
    ##       Bacteroidales 
    ##                 920

``` r
reclassify <- function(df,taxanames,phylum_order){
  df_to_change = df %>% filter(Genus %in% taxanames)
  df = df %>% filter(!Genus %in% taxanames)
  df_to_change$Phylum <- phylum_order[[1]]
  df_to_change$Order <- phylum_order[[2]]
  df_to_change$Family <- NA
  df_to_change$Genus <- NA
  df_to_change$Species <- NA
  df=rbind(df,df_to_change)
  return(df)
  }

diamond_genus = reclassify(diamond_genus,
                           taxanames = c('Bacteroidia','Bacteroidetes'),
                           phylum_order = c('Bacteroidetes',NA))
diamond_genus = reclassify(diamond_genus,
                           taxanames = c('Flavobacteriales'),
                           phylum_order = c('Bacteroidetes','Flavobacteriales'))
diamond_genus = reclassify(diamond_genus,
                           taxanames = c('Bacteroidaceae','Bacteroidales','Muribaculaceae',
                                         'Paludibacteraceae','Porphyromonadaceae','Mediterranea',
                                         'Dysgonamonadaceae','Dysgonomonadaceae','Gallibacteroides',
                                         'Hallella','Avibacteroides','Limisoma','Rikenellaceae',
                                         'Caccoplasma'),
                           phylum_order = c('Bacteroidetes','Bacteroidales'))
diamond_genus = reclassify(diamond_genus,
                           taxanames = c('Firmicutes','Erysipelotrichia','Clostridia','Solobacterium',
                                         'Clostridiales','Lachnospiraceae','Oscillospiraceae',
                                         'Tyzzerella'),
                           phylum_order = c('Firmicutes',NA))
diamond_genus = reclassify(diamond_genus,
                           taxanames = c('Proteobacteria','Shigella','Dechloromonas',
                                         'Alphaproteobacteria','Caedimonadaceae'),
                           phylum_order = c('Proteobacteria',NA))
diamond_genus = reclassify(diamond_genus,
                           taxanames = c('bacterium','Bacteria','organism','Lentisphaeria'),
                           phylum_order = c(NA,NA))
diamond_genus = reclassify(diamond_genus,
                           taxanames = c('Microviridae','Podoviridae','Bacteriophage',
                                         'Siphoviridae','Myoviridae','CrAss-like','virus'),
                           phylum_order = c('Bacteriophage',NA))

diamond_notfound = diamond_genus %>% filter(is.na(Phylum))
sort(table(diamond_notfound$Genus))
```

    ## integer(0)

``` r
table(diamond_genus$Phylum)
```

    ## 
    ##    Actinobacteria     Bacteriophage     Bacteroidetes  Campylobacterota 
    ##                56               664             57309                13 
    ##    Fibrobacterota        Firmicutes    Proteobacteria Verrucomicrobiota 
    ##                 1               410               917                12

``` r
#for an individual GeneIDs, determine taxonomic distribution of hits to nr database
diamond_genus = diamond_genus  %>% group_by(GeneID) %>%
                     mutate(found_in_phage = if_else('Bacteriophage' %in% Phylum,1,0),
                     multiphyla =if_else(
                       length(setdiff(unique(Phylum),c('Bacteriophage',NA)))>1,1,0),
                     multiorder = if_else(length(setdiff(unique(Order),NA))>1|multiphyla==1,1,0),
                     multigenus = if_else(length(setdiff(unique(Genus),NA))>1|multiorder==1,1,0),
                     multispecies = if_else(
                       length(setdiff(Species,c('sp.',NA)))>1|multigenus==1,1,0),
                     singlespecies = if_else(length(setdiff(Species,c('sp.',NA)))==1,1,0))

diamond_genus = diamond_genus  %>% group_by(GeneID) %>%
                    mutate( 
                      #assign taxonomic categories
                      cat = ifelse(multiphyla==1,'multiphyla',
                            ifelse(multiorder==1,'multiorder',
                              ifelse(multigenus==1,'multigenus',
                                ifelse(multispecies==1,'multispecies',
                                  ifelse(singlespecies ==1,'singlespecies','unclassified'))))),
                      
                      #assign most conserved taxonomic rank
                      group = ifelse(multiphyla==1,'Bacteria',
                          ifelse(
                            multiorder==1,toString(setdiff(unique(Phylum),c(NA,'Bacteriophage'))),
                              ifelse(multigenus==1,toString(setdiff(unique(Order),NA)),
                                  ifelse(multispecies==1,toString(setdiff(unique(Genus),NA)),
                                      ifelse(singlespecies==1,
                                            paste(
                                                toString(setdiff(unique(Genus),NA)),
                                                toString(setdiff(unique(Species),c('sp.',NA))),
                                                sep='_'),
                                            'unclassified'))))),
                      #creates lists for individual taxonomic ranks
                      phylalist = toString(setdiff(unique(Phylum),NA)),
                      orderlist = toString(setdiff(unique(Order),NA)),
                      genuslist = toString(setdiff(unique(Genus),NA)),
                      specieslist = toString(setdiff(unique(Species),NA))
                      )
#some hits were unclassified ie hitting to bacteriophage or Bacteroides_sp. but no classified taxa
unclassified = diamond_genus %>% filter(group=='unclassified') 
unclassified = unclassified %>% mutate(group = if_else(!is.na(Genus),Genus_Species,Phylum))
classified = diamond_genus %>% filter(group!='unclassified')                                    
diamond_genus =  rbind(classified,unclassified) 

diamond_genus_short = diamond_genus %>%
  group_by(GeneID) %>%
    slice_max(order_by='pident',with_ties=F,n=1) %>% #get a single representative per GENEID
    select(-c(pident:singlespecies)) #remove columns specific to individual hits

write_tsv(diamond_genus_short,file.path(diamond_dir,'diamond_blast_taxonomy.txt'))
```
