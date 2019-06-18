#!/usr/bin/env Rscript

## PREP
# clean up
rm(list=ls())

# load functions
source("https://raw.githubusercontent.com/boopsboops/reference-libraries/master/scripts/funs.R")

# load reference lib
source("https://raw.githubusercontent.com/boopsboops/reference-libraries/master/scripts/references-load.R")
# loads objects: uk.species.table, uk.species.table.common, reflib.orig

# make a copy so don't have to keep reloading the reflib
reflib <- reflib.orig


## SUBSET BY MARKER

# list prefixes
prefixes.all <- c("coi.lerayxt.noprimers","coi.seamid.noprimers","coi.seashort.noprimers","coi.ward.noprimers","12s.miya.noprimers","12s.riaz.noprimers","12s.valentini.noprimers","12s.taberlet.noprimers","16s.berry.noprimers","cytb.minamoto.noprimers")

# subset by marker and regenerate list of dataframes
all.dbs <- mapply(function(x) subset_by_marker(x, df=reflib, thresh=0.5), prefixes.all, SIMPLIFY=FALSE, USE.NAMES=TRUE)


## OPTIONS FOR COMMON AND INTERSECTED ANALYSES - RUN ...

# get shared spp list
spp.list <- mapply(function(x) unique(reflib$sciNameValid[!is.na(reflib[[paste0("nucleotidesFrag.",x)]])]), prefixes.all)
spps.intersect <- Reduce(intersect, spp.list)
length(spps.intersect)# 221 shared spp
length(intersect(spps.intersect,uk.species.table.common$sciName))# 88 common shared spp

# !!!! subset only shared species spps.intersect !!!!
all.dbs <- mapply(function(x) x %>% filter(sciNameValid %in% spps.intersect), all.dbs, SIMPLIFY=FALSE, USE.NAMES=TRUE)

# !!!! subset common species !!!! # 176 common spp.
all.dbs <-mapply(function(x) x %>% filter(sciNameValid %in% uk.species.table.common$sciName), all.dbs, SIMPLIFY=FALSE, USE.NAMES=TRUE)

# !!! subset those those without 'mitochondri*' search term
#reflib %>% filter(source=="GENBANK") %>% filter(!grepl("mitochondr",notesGenBank)) 
all.dbs <-mapply(function(x) x %>% filter(source=="GENBANK" & grepl("mitochondr",notesGenBank) | source=="BOLD"), all.dbs, SIMPLIFY=FALSE, USE.NAMES=TRUE)


## NOW GENERATE STATS

# stats total
numSequences.df <- enframe(mapply(function(x) x %>% summarise(n=n()) %>% pull(n), all.dbs),name="marker",value="numSequences")

# median frag len
medianLength.df <- enframe(mapply(function(x,y) x %>% summarise(median=median(!!as.name(paste0("lengthFrag.",y)))) %>% pull(median), all.dbs, prefixes.all),name="marker",value="medianLength")

# num species
numSpecies.df <- enframe(mapply(function(x) x %>% summarise(n=n_distinct(sciNameValid)) %>% pull(n), all.dbs),name="marker",value="numSpecies")

# stats per species
meanIndPerSpecies.df <- enframe(mapply(function(x) x %>% group_by(sciNameValid) %>% summarise(n=n()) %>% summarise(mean=mean(n)) %>% pull(mean), all.dbs),name="marker",value="meanIndPerSpecies")
medianIndPerSpecies.df <- enframe(mapply(function(x) x %>% group_by(sciNameValid) %>% summarise(n=n()) %>% summarise(median=median(n)) %>% pull(median), all.dbs),name="marker",value="medianIndPerSpecies")

# collapse haplotypes species-wise (only duplicate haplotypes with the same species name are excluded)
all.dbs.collapsed <- mapply(function(x,y) bind_rows(mcmapply(FUN=function(z) hap_collapse_df(df=z,lengthcol=paste0("lengthFrag.",y),nuccol=paste0("nucleotidesFrag.",y)), split(x,x$sciNameValid), SIMPLIFY=FALSE, mc.cores=8)), all.dbs, prefixes.all, SIMPLIFY=FALSE, USE.NAMES=TRUE)

# calculate the names of seqs with identical haplotypes
all.dbs.collapsed.spp <- mapply(function(x,y) mcmapply(FUN=function(z) get_sames(df=x,ids="dbid",nucs=paste0("nucleotidesFrag.",y),sppVec="sciNameValid",query=z), x[[paste0("nucleotidesFrag.",y)]], mc.cores=8, SIMPLIFY=FALSE, USE.NAMES=FALSE), all.dbs.collapsed, prefixes.all, SIMPLIFY=FALSE, USE.NAMES=TRUE)

# run on the data
all.dbs.collapsed.haps <- mapply(function(x,y) x %>% mutate(nMatches=sapply(y, function(z) length(unique(z))), matchTax=sapply(y, function(w) paste(unique(w),collapse=" | "))), all.dbs.collapsed, all.dbs.collapsed.spp, SIMPLIFY=FALSE, USE.NAMES=TRUE)


# calc the overall rates
idRate.df <- enframe(mapply(function(x) x %>% group_by(sciNameValid) %>% summarise(perc=length(which(nMatches==1))/length(nMatches)) %>% summarise(nn=length(which(perc==1))/length(perc)) %>% pull(nn), all.dbs.collapsed.haps),name="marker",value="idRate")

# calc the rates per species
all.dbs.collapsed.haps.names <- mapply(function(x) x %>% group_by(sciNameValid) %>% summarise(perc=length(which(nMatches==1))/length(nMatches)), all.dbs.collapsed.haps, SIMPLIFY=FALSE)
all.dbs.collapsed.haps.names[[1]] %>% print(n=Inf) 

# make a joined df
list.df <- list(numSequences.df,
    medianLength.df,
    numSpecies.df,
    meanIndPerSpecies.df,
    medianIndPerSpecies.df,
    idRate.df)

# bind all
joined.df <- list.df %>% purrr::reduce(full_join, by="marker") %>% arrange(marker) 
joined.df %>% print()
dir.create("../temp/id-rates")
joined.df %>% write_csv(path="../temp/id-rates/id-results-common-mito.csv")
# change filename for other subsets
#id-results-all.csv
#id-results-intersect.csv
#id-results-intersect-common.csv
#id-results-common.csv
#id-results-all-mito.csv
#id-results-common-mito.csv
