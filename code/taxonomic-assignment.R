#!/usr/bin/env Rscript

# load libs
library("tidyverse")
library("magrittr")
library("lubridate")
library("stringr")
library("ape")
require("stringdist")
require("phangorn")
library("parallel")
#library("xtable")
library("digest")
library("rfishbase")
library("taxize")
# other sources
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/tab2fas.R")

# set marker and proj dir 
marker <- "miya"
marker <- "leray-xt"
marker <- "sea-short"
marker <- "sea-mid"
#
rel.path <- "../temp/fastq"
proj.path <- paste(rel.path,marker,sep="/") 


################ STEP 1 -- GLOBAL BLAST ####################
################ STEP 1 -- GLOBAL BLAST ####################
################ STEP 1 -- GLOBAL BLAST ####################
# rm(list=ls())
# first run HMM homology search in 'blast.sh'
# first run HMM homology search in 'blast.sh'
# first run HMM homology search in 'blast.sh'

# read asv table back in written table to get format correct
otu.tab <- as_tibble(read.table(file=paste0(proj.path,"/results/asv-table.tsv"), sep="\t", header=TRUE, as.is=TRUE, row.names=1, check.names=FALSE), rownames="asv")
dnas.curated <- read.FASTA(file=paste0(proj.path,"/results/asvs.fna"))

# load up data cleaned by the HMM
dnas.curated.clean <- read.FASTA(file=paste0(proj.path,"/results/asvs-clean.fasta"))
dnas.curated.dirty <- dnas.curated[which(!(names(dnas.curated) %in% names(dnas.curated.clean)))]
#write.FASTA(dnas.curated.dirty, file=paste0(proj.path,"/results/otus-dirty.fasta"))
# subset the clean sequences from the ASV table
otu.tab.clean <- otu.tab %>% filter(asv %in% names(dnas.curated.clean))
# write out
write_tsv(otu.tab.clean,path=paste0(proj.path,"/results/asv-table-clean.tsv"))

# get numbers of seqs lost during hmm search
total <- otu.tab %>% summarise_if(is.numeric, sum, na.rm=TRUE) %>% rowSums()
lost <- otu.tab %>% filter(!asv %in% names(dnas.curated.clean)) %>% summarise_if(is.numeric, sum, na.rm=TRUE) %>% rowSums()
retained <- otu.tab.clean %>% summarise_if(is.numeric, sum, na.rm=TRUE) %>% rowSums()
print(retained)
lost+retained == total

# prep for BLAST
# break up fasta for faster blast
div <- length(dnas.curated.clean)/8
chunk <- ceiling(div)
spl <- split(dnas.curated.clean, ceiling(seq_along(dnas.curated.clean)/chunk))
names(spl) <- paste0("blast",seq_along(spl),".fa")
mapply(function(x,i) write.FASTA(x, file=paste0(proj.path,"/results/",i)), spl, names(spl), USE.NAMES=FALSE)


# NOW run BLAST in `blast.sh` #
# NOW run BLAST in `blast.sh` #
# NOW run BLAST in `blast.sh` #


# read in and filter the BLAST results
blast.res <- read_tsv(file=paste0(proj.path,"/results/blast-results.tsv"),guess_max=100000)

# select "best" blast hit
blast.res.sorted <- blast.res %>% 
    group_by(qseqid) %>% 
    arrange(desc(bitscore),score,evalue,desc(pident),desc(nident),.by_group=TRUE) %>% 
    slice(1) %>%
    ungroup()

# remake the otus frame with DNAs
otus.df <- tibble(names=names(dnas.curated.clean), dnas=sapply(as.character(dnas.curated.clean), paste, collapse="")) %>% mutate(len=str_length(dnas))

# to make the annotated Blast-ASV table
# add read numbers to OTU table
# get number of reads per sample
otu.tab.clean %>% select_if(is.numeric) %>% summarise_all(list(sum))
# number reads per OTU
otu.tab.clean.sums <- otu.tab.clean %>% select_if(is.numeric) %>% mutate(sums=rowSums(.))
otu.tab.clean.sums %>% select(sums)

# make new table 
otu.reads <- tibble(qseqid=otu.tab.clean$asv, totalReads=otu.tab.clean.sums$sums)
otu.reads <- otu.reads %>% mutate(propTotal=totalReads/sum(totalReads))

# change var names of otu list and merge with blast results and reads tables
otus.df.curated <- otus.df %>% rename(qseqid=names, asvLength=len)
asv.results <- left_join(left_join(otus.df.curated, blast.res.sorted, by="qseqid"), otu.reads, by="qseqid")

# clean up taxids
asv.results %<>% mutate(staxids=str_replace_all(staxids, ";.*", ""), sscinames=str_replace_all(sscinames, ";.*", ""), scomnames=str_replace_all(scomnames, ";.*", ""))

# add taxonomy
ids <- unique(asv.results$staxids)

# run taxize classification
# be sure to add "ENTREZ_KEY=ncbi.key" to ~/.Renviron
# takes about a min for 1,000 (may need to repeat a few times)
taxize::classification(ids[1], db='ncbi', check=FALSE)
ght <- mclapply(ids, taxize::classification, db='ncbi', check=FALSE, mc.cores=1)
# check - should all be class classification
table(sapply(ght,class))

# join and make mini frames
hj <- do.call(rbind, ght)

# filter ranks
kin <- hj %>% filter(rank=="kingdom")
phy <- hj %>% filter(rank=="phylum")
cla <- hj %>% filter(rank=="class")
ord <- hj %>% filter(rank=="order")
fam <- hj %>% filter(rank=="family")
gen <- hj %>% filter(rank=="genus")

# add to a table
asv.results %<>% mutate(
    kingdom=kin$name[match(asv.results$staxids, kin$query)], 
    phylum=phy$name[match(asv.results$staxids, phy$query)], 
    class=cla$name[match(asv.results$staxids, cla$query)], 
    order=ord$name[match(asv.results$staxids, ord$query)], 
    family=fam$name[match(asv.results$staxids, fam$query)], 
    genus=gen$name[match(asv.results$staxids, gen$query)])

# rearrange
asv.results %<>% 
    select(qseqid,sseqid,staxids,sacc,sskingdoms,kingdom,phylum,class,order,family,genus,sscinames,scomnames,asvLength,evalue,score,bitscore,length,pident,nident,totalReads,propTotal,dnas) %>% 
    mutate(marker=marker) %>% 
    arrange(desc(totalReads))

# write out to view
write_tsv(asv.results, path=paste0(proj.path,"/results/asv-results-sorted.tsv"))
#asv.results <- read_tsv(file=paste0(proj.path,"/results/asv-results-sorted.tsv"))


################ STEP 2 -- BLAST LOCAL FISH LIBRARY ####################
################ STEP 2 -- BLAST LOCAL FISH LIBRARY ####################
################ STEP 2 -- BLAST LOCAL FISH LIBRARY ####################
# rm(list=ls())
# import reference library
# import reflib and exclusions file
# load reference lib and clean
source("https://raw.githubusercontent.com/boopsboops/reference-libraries/master/scripts/references-load.R")
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/tab2fas.R")
ref.lib <- reflib.orig

# set marker and proj dir 
marker <- "miya"
marker <- "leray-xt"
marker <- "sea-short"
marker <- "sea-mid"
#
rel.path <- "../temp/fastq"
proj.path <- paste(rel.path,marker,sep="/") 

# create fasta and write out the reference library
ref.fas <- tab2fas(ref.lib, seqcol="nucleotides", namecol="dbid")
write.FASTA(ref.fas, file=paste0(rel.path,"/references.fasta"))

# read in 1st step blast
asv.results <- read_tsv(file=paste0(proj.path,"/results/asv-results-sorted.tsv"))

# load fishbase taxonomy and filter out fish reads
data(fishbase)
head(fishbase)
asv.results %<>% filter(class=="Actinopteri" | class=="Cephalaspidomorphi" | class=="Chondrichthyes" | order %in% fishbase$Order | family %in% fishbase$Family | genus %in% fishbase$Genus) %>% 
    filter(kingdom=="Metazoa")
print(asv.results)

# write out
write.FASTA(tab2fas(df=asv.results, seqcol="dnas", namecol="qseqid"), file=paste0(proj.path,"/results/fishqueries.fasta"))


# now run local fish only BLAST from `blast.sh` #
# now run local fish only BLAST from `blast.sh` #
# now run local fish only BLAST from `blast.sh` #

# load blast result
local.db.blast <- read_tsv(paste0(proj.path,"/results/fish-blast-result.tsv"))

# chose "best" hit based on bitscore
# also add scinames
local.db.blast.sorted <- local.db.blast %>%
    group_by(qseqid) %>%
    arrange(desc(bitscoreLocal),.by_group=TRUE) %>%
    filter(bitscoreLocal==max(bitscoreLocal)) %>%
    mutate(sciNameValidLocal=ref.lib$sciNameValid[match(sseqidLocal,ref.lib$dbid)]) %>%
    arrange(sciNameValidLocal,.by_group=TRUE) %>%
    mutate(sciNameValidLocal=paste(unique(sciNameValidLocal),collapse="; ")) %>%
    slice(1) %>% 
    ungroup()
local.db.blast.sorted %>% arrange(desc(bitscoreLocal)) %>% print(n=Inf)

# write out for later 
local.db.blast.sorted %>% write_csv(paste0(proj.path,"/results/fish-blast-result-sorted.csv"))


################ STEP 3 -- EPA ASSIGNMENT ####################
################ STEP 3 -- EPA ASSIGNMENT ####################
################ STEP 3 -- EPA ASSIGNMENT ####################
# rm(list=ls())
# join the reflib sequences and the query sequences
# list prefixes
prefix <- "12s.miya.noprimers"
prefix <- "coi.lerayxt.noprimers"
prefix <- "coi.seamid.noprimers"
prefix <- "coi.seashort.noprimers"

# also marker 
marker <- "miya"
marker <- "leray-xt"
marker <- "sea-mid"
marker <- "sea-short"

#
rel.path <- "../temp/fastq"
proj.path <- paste(rel.path,marker,sep="/") 


# reflib
source("https://raw.githubusercontent.com/boopsboops/reference-libraries/master/scripts/references-load.R")
source("https://raw.githubusercontent.com/boopsboops/reference-libraries/master/scripts/funs.R")
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/tab2fas.R")
# change reflib
ref.lib <- reflib.orig

# subset the marker from the reflib - using 50% length cutoff
reflib.sub <- subset_by_marker(prefix=prefix,df=ref.lib,thresh=0.5)

# collapse haps
reflib.red <- bind_rows(mcmapply(FUN=function(x) hap_collapse_df(df=x,lengthcol=paste0("lengthFrag.",prefix),nuccol=paste0("nucleotidesFrag.",prefix)), split(reflib.sub,reflib.sub$sciNameValid), SIMPLIFY=FALSE, mc.cores=8))

# load up the queries
fish.queries <- read.FASTA(file=paste0(proj.path,"/results/fishqueries.fasta"))

# concatenate with the reference and align
combined.fas <- c(fish.queries,tab2fas(reflib.red, seqcol=paste0("nucleotidesFrag.",prefix),namecol="dbid"))

# remove Ns
combined.fas <- rm_ns(bin=combined.fas)

# align
combined.fas.aligned <- mafft(combined.fas,path="mafft",method="retree 1")

# make a subdir
dir.create(paste0(proj.path,"/results/epa"))
#write.FASTA(combined.fas.aligned, file=paste0(proj.path,"/results/epa/combined.aligned.fasta"))

# make a reference and query alignment
references.aligned <- combined.fas.aligned[grep("asv[0-9]+",rownames(combined.fas.aligned),invert=TRUE),]
queries.aligned <- combined.fas.aligned[grep("asv[0-9]+",rownames(combined.fas.aligned)),]

write.FASTA(references.aligned, file=paste0(proj.path,"/results/epa/references.aligned.fasta"))
write.FASTA(queries.aligned, file=paste0(proj.path,"/results/epa/queries.aligned.fasta"))

# generate taxonomy file
reflib.red %>% mutate(sciNameValid=str_replace_all(sciNameValid," ", "_")) %>% 
    mutate(taxonomy=paste(subphylum,class,order,family,genus,sciNameValid,sep=";")) %>% 
    select(dbid,taxonomy) %>% 
    write_tsv(path=paste0(proj.path,"/results/epa/references.taxonomy.tsv"),col_names=FALSE)


# NOW run EPA and GAPPA from the 'epa.sh' script #
# NOW run EPA and GAPPA from the 'epa.sh' script #
# NOW run EPA and GAPPA from the 'epa.sh' script #

# load up epa results
epa.results <- read_tsv(file=paste0(proj.path,"/results/epa/per_query.tsv"))
epa.results %<>% rename(qseqid=name)
# LWR: likelihood weight that was assigned to this exact taxonomic path
# fract: LWR divided by the global total likelihood weight
# aLWR: accumulated likelihood weights that were assigned either to this taxonomic path or any taxonomic path below this
# afract: aLWR divided by the global total likelihood weight
# taxopath: the taxonomic path

# best the best id
epa.results.best <- epa.results %>% 
    group_by(qseqid) %>%
    arrange(desc(LWR),.by_group=TRUE) %>% #print(n=300)
    slice(1) %>% 
    mutate(epaID=sapply(str_split(taxopath,";"),last), epaID=str_replace_all(epaID,"_"," ")) %>% 
    ungroup()

# get species IDs
epa.results.species <- epa.results %>% 
    group_by(qseqid) %>% 
    filter(grepl("_",taxopath)) %>% 
    arrange(taxopath,.by_group=TRUE) %>% 
    mutate(epaBestSpp=sapply(str_split(taxopath,";"),last), epaBestSpp=str_replace_all(epaBestSpp,"_"," ")) %>%  
    mutate(epaAllSpp=paste(unique(epaBestSpp),collapse="; ")) %>% 
    arrange(desc(LWR),.by_group=TRUE) %>% 
    slice(1) %>% 
    ungroup() %>%
    select(qseqid,epaBestSpp,epaAllSpp)

# join and tidy
epa.results.filtered <- left_join(epa.results.best,epa.results.species) %>% 
    select(qseqid,LWR,epaID,epaBestSpp,epaAllSpp) %>% 
    rename(epaIdLWR=LWR)

# write out 
epa.results.filtered %>% write_csv(path=paste0(proj.path,"/results/epa/epa-results-filtered.csv"))


################ STEP 4 -- SWARM CLUSTERING ####################
################ STEP 4 -- SWARM CLUSTERING ####################
################ STEP 4 -- SWARM CLUSTERING ####################
# rm(list=ls())
# also marker 
marker <- "leray-xt"
marker <- "sea-mid"
marker <- "sea-short"
marker <- "miya"
#
rel.path <- "../temp/fastq"
proj.path <- paste(rel.path,marker,sep="/") 

# make a file for swarm input
fish.queries <- read.FASTA(file=paste0(proj.path,"/results/fishqueries.fasta"))
asv.results <- read_tsv(file=paste0(proj.path,"/results/asv-results-sorted.tsv"))
# make names
asv.results %<>% mutate(names=paste0(qseqid,";size=",totalReads))
names(fish.queries) <- asv.results$names[match(names(fish.queries),asv.results$qseqid)]
# write out
dir.create(paste0(proj.path,"/results/swarm"))
write.FASTA(fish.queries, paste0(proj.path,"/results/swarm/fishqueries-clean-abundances.fasta"))

# NOW run SWARM from the 'swarm.sh' script #
# NOW run SWARM from the 'swarm.sh' script #
# NOW run SWARM from the 'swarm.sh' script #


# load up swarm results - ignore errors
swarm.tab <- read_table2(file=paste0(proj.path,"/results/swarm/swarm.tsv"),col_names=paste0("col",seq_along(1:1000)))

swarm.tab %<>% 
    gather(key=column,value=qseqid,-col1,na.rm=TRUE) %>% 
    select(-column) %>% rename(swarm=col1) %>% 
    arrange(swarm)# %>% print(n=Inf)

swarm.tab %>% write_csv(path=paste0(proj.path,"/results/swarm/swarm-tabular.csv"))


################ STEP 5 -- JOIN THE RESULTS ####################
################ STEP 5 -- JOIN THE RESULTS ####################
################ STEP 5 -- JOIN THE RESULTS ####################
# rm(list=ls())
# also marker 
marker <- "leray-xt"
marker <- "sea-mid"
marker <- "sea-short"
marker <- "miya"
#
rel.path <- "../temp/fastq"
proj.path <- paste(rel.path,marker,sep="/") 


# load global blast, local blast, epa, swarm results
data(fishbase)
fish.otus <- read_tsv(file=paste0(proj.path,"/results/asv-results-sorted.tsv")) %>% filter(class=="Actinopteri" | class=="Cephalaspidomorphi" | class=="Chondrichthyes" | order %in% fishbase$Order | family %in% fishbase$Family | genus %in% fishbase$Genus) %>% 
   filter(kingdom=="Metazoa") %>% 
   select(-marker)
local.db.blast.sorted <- read_csv(paste0(proj.path,"/results/fish-blast-result-sorted.csv"))
epa.results <- read_csv(file=paste0(proj.path,"/results/epa/epa-results-filtered.csv"))
swarm.results <- read_csv(file=paste0(proj.path,"/results/swarm/swarm-tabular.csv"))

# combine
taxonomy.results <- purrr::reduce(list(fish.otus,local.db.blast.sorted,epa.results,swarm.results),left_join, by="qseqid")

# Mode FUN
mode_avg <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}
minLen <- ceiling(mode_avg(taxonomy.results$asvLength)*0.9)
# plot(sort(taxonomy.results$bitscoreLocal))
# plot(sort(taxonomy.results$pidentLocal))

# process
# two step
# either (a) or (b) must be satisfied
# (a) epaAssign = highest likelihood EPA id and blast id the same, EPA likelihood > 0.9, and pident > 0.9
# (b) blastAssign = the best species-level EPA id (not necesarily best overall id) and blast id the same, and pident >0.97 and lengthLocal >0.9 of modal ASV length 
taxonomy.results %<>% 
    mutate(epaAssign=if_else(epaID==sciNameValidLocal & epaIdLWR>=0.9 & pidentLocal>=90,TRUE,FALSE)) %>% 
    mutate(blastAssign=if_else(str_detect(sciNameValidLocal, epaBestSpp) & lengthLocal>=minLen & pidentLocal>=97,TRUE,FALSE)) %>% 
    mutate(assigned=if_else(epaAssign==TRUE | blastAssign==TRUE,TRUE,FALSE)) %>% 
    mutate(marker=marker)
#taxonomy.results %>% write_csv(path=paste0(proj.path,"/results/swarm/sanity-check.csv"))

# take a look
taxonomy.results %>% 
    select(qseqid,totalReads,asvLength,lengthLocal,pidentLocal,nidentLocal,bitscoreLocal,sciNameValidLocal,epaIdLWR,epaID,epaBestSpp,epaAllSpp,epaAssign,blastAssign,assigned) %>% #,matchNames,assigned 
    arrange(desc(assigned),desc(epaAssign),desc(blastAssign),desc(bitscoreLocal)) %>% 
    #write_csv(path=paste0(proj.path,"/results/sanity-check.csv"))
    print(n=500)

# save
taxonomy.results %>% write_csv(path=paste0(proj.path,"/results/taxonomic-assignment.csv"))
# read back in 
#taxonomy.results <- read_csv(file=paste0(proj.path,"/results/taxonomic-assignment.csv"))


################ STEP 6 -- PER SAMPLE STATS ####################
################ STEP 6 -- PER SAMPLE STATS ####################
################ STEP 6 -- PER SAMPLE STATS ####################
# rm(list=ls())
# also marker 
marker <- "leray-xt"
marker <- "sea-mid"
marker <- "sea-short"
marker <- "miya"
#
rel.path <- "../temp/fastq"
proj.path <- paste(rel.path,marker,sep="/") 


# load tables
plates <- read_csv(file="../data/sample-plates.csv")
taxonomy.results <- read_csv(file=paste0(proj.path,"/results/taxonomic-assignment.csv"))
otu.tab <- read.table(file=paste0(proj.path,"/results/asv-table.tsv"), sep="\t", header=TRUE, as.is=TRUE, row.names=1, check.names=FALSE)
curated.table.fish <- as_tibble(otu.tab, rownames="asv")

# filter, and replace the names in the otu table
taxonomy.results.assigned <- taxonomy.results %>% filter(assigned==TRUE)
curated.table.fish %<>% filter(asv %in% taxonomy.results.assigned$qseqid) %>% 
    mutate(asv=taxonomy.results.assigned$sciNameValidLocal[match(asv,taxonomy.results.assigned$qseqid)])

# collapse by species
curated.table.fish %<>% group_by(asv) %>% summarise_all(list(sum))

# get the plates codes
# filter by marker
plates %<>% filter(primerSet==marker) %>% mutate(eventMD5=str_trunc(sapply(eventID, digest, algo="md5", serialize=FALSE, USE.NAMES=FALSE), width=12, side="right", ellipsis=""))

# process the table 
curated.fish.col <- curated.table.fish %>% 
    # reshape into column format
    gather(key=sample,value=nreads,-asv) %>% 
    # make a hash and a rep-PCR column
    mutate(hash=str_split_fixed(sample,"\\.",3)[,1], rep=paste(str_split_fixed(sample,"\\.",3)[,2],str_split_fixed(sample,"\\.",3)[,3],sep=".")) %>% 
    # remove zero reads and rename species column
    filter(nreads > 0) %>% rename(species=asv) %>% 
    # add site names -- RUN first half of the prep-barcodes.R script to get the table 'plates'
    mutate(eventID=plates$eventID[match(hash,plates$eventMD5)]) %>% 
    # clean 
    select(-hash,-sample) %>% 
    # add marker
    mutate(marker=marker) %>% 
    # select cols
    select(marker,species,nreads,rep,eventID)

# write out the final table
curated.fish.col %>% write_csv(path=paste0(proj.path,"/results/fish-final-table.csv"))


################ STEP 6 -- JOIN ALL RESULTS ####################
################ STEP 6 -- JOIN ALL RESULTS ####################
################ STEP 6 -- JOIN ALL RESULTS ####################

# run all steps for all markers before running this step

# taxon assignments
bind_rows(read_csv(file=paste0(rel.path,"/sea-short/results/taxonomic-assignment.csv"),col_types=cols(sseqidLocal=col_character())),
    read_csv(file=paste0(rel.path,"/sea-mid/results/taxonomic-assignment.csv"),col_types=cols(sseqidLocal=col_character())),
    read_csv(file=paste0(rel.path,"/miya/results/taxonomic-assignment.csv"),col_types=cols(sseqidLocal=col_character())),
    read_csv(file=paste0(rel.path,"/leray-xt/results/taxonomic-assignment.csv"),col_types=cols(sseqidLocal=col_character()))) %>% 
    write_csv(path="../data/fish-assignments.csv")

# otu tables
bind_rows(read_csv(file=paste0(rel.path,"/sea-short/results/fish-final-table.csv")),
    read_csv(file=paste0(rel.path,"/sea-mid/results/fish-final-table.csv")),
    read_csv(file=paste0(rel.path,"/miya/results/fish-final-table.csv")),
    read_csv(file=paste0(rel.path,"/leray-xt/results/fish-final-table.csv"))) %>% 
    write_csv(path="../data/results-by-marker.csv")

# otu asv blast table
bind_rows(read_tsv(file=paste0(rel.path,"/sea-short/results/asv-results-sorted.tsv")),
    read_tsv(file=paste0(rel.path,"/sea-mid/results/asv-results-sorted.tsv")),
    read_tsv(file=paste0(rel.path,"/miya/results/asv-results-sorted.tsv")),
    read_tsv(file=paste0(rel.path,"/leray-xt/results/asv-results-sorted.tsv"))) %>% 
    write_csv(path="../data/blast-results-all.csv")
