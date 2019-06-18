#!/usr/bin/env Rscript
# Script to generate primer amplification success stats for metabarcoding markers using whole mitogenomes #
# Uses MFEprimer #

## load funs and reflibs
# clean
rm(list=ls())

# load functions
source("https://raw.githubusercontent.com/boopsboops/reference-libraries/master/scripts/funs.R")

# load reference lib
source("https://raw.githubusercontent.com/boopsboops/reference-libraries/master/scripts/references-load.R")
# loads objects: uk.species.table, uk.species.table.common, reflib.orig

# make a copy so don't have to keep reloading the reflib
reflib <- reflib.orig


## make primer fastas
# load up the primers
primers <- read_csv(file="../data/primers.csv")

# make names for the primers
primers %<>% mutate(name=paste(locus,reference,primer,sep="."))

# make fas for all primers and write out
primers.list <- lapply(split(primers,primers$name), function(x) tab2fas(df=x,seqcol="oligo",namecol="officialName"))

# remove orig Leray-XT, use modified
primers.list <- primers.list[-grep("orig", names(primers.list))]

# write out
# if required create a directory to run MFEprimer
dir.create(path="../temp/mitogenome")
mapply(function(x,y) write.FASTA(x, file=paste0("../temp/mitogenome/",y,".primers.fas")), primers.list, names(primers.list), USE.NAMES=FALSE)


## reference library filter
# keep only the mitogenomes 
reflib %<>% filter(length > 15000) 
plot(table(reflib$length))

# collapse haps by species 
reflib.haps <- bind_rows(mcmapply(hap_collapse_df, lengthcol="length",nuccol="nucleotides", split(reflib,reflib$sciNameValid), SIMPLIFY=FALSE, mc.cores=8))

# check numbers of indivs/spp
reflib.haps %>% summarise(n=n())
reflib.haps %>% distinct(sciNameValid)

# write out a nucleotides file
write.FASTA(tab2fas(df=reflib.haps,seqcol="nucleotides",namecol="dbid"), file="../temp/mitogenome/mitogenome.fasta")


### NOW RUN IN MFEPRIMER 'MFEprimer.sh' ###
### NOW RUN IN MFEPRIMER 'MFEprimer.sh' ###
### NOW RUN IN MFEPRIMER 'MFEprimer.sh' ###


# Process MFE results and estimate primer bias

# load prefixes
prefixes <- c(
    "12s.miya.universal",
    "12s.miya.elasmo",
    "12s.valentini.tele01",
    "12s.riaz.v5",
    "12s.taberlet.tele02",
    "12s.taberlet.elas02",
    "coi.leray.xt",
    "coi.sea.mid",
    "coi.sea.short",
    "coi.ward.fish",
    "16s.berry.fish",
    "cytb.minamoto.fish")

# load up all MFE data
mfe.list <- lapply(prefixes, function(x) read_tsv(file=paste0("../temp/mitogenome/",x,".primers.fas.results.tsv"),col_names=TRUE))
names(mfe.list) <- prefixes

# run function over all results files
mfe.annotated.all <- mapply(function(x,y) process_MFE(mfe=x,primers=primers,prefixes=y,references=reflib.haps, common=uk.species.table.common), mfe.list, prefixes, SIMPLIFY=FALSE)
# join
mfe.annotated.all <- bind_rows(mfe.annotated.all)
# Total 184 spp.

# run and summarise results per group
#fun to get prop
prop_fun <- function(x){length(which(x>0))/length(x)}

# run 
mfe.annotated.all %>% 
    group_by(marker,class) %>% 
    summarise(mean=mean(bestPPC), prop=prop_fun(bestPPC)) %>% 
    filter(class=="Actinopterygii" | class=="Elasmobranchii") %>% 
    arrange(class,marker) %>% write_csv(path="../temp/mitogenome/all.primer.results.csv") 

# for just common spp
mfe.annotated.all %>% 
    filter(common=="common") %>% 
    group_by(marker,class) %>% 
    summarise(mean=mean(bestPPC), prop=prop_fun(bestPPC)) %>% 
    filter(class=="Actinopterygii" | class=="Elasmobranchii") %>% 
    arrange(class,marker) %>% write_csv(path="../temp/mitogenome/common.primer.results.csv")
