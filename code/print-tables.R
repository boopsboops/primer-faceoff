#!/usr/bin/env Rscript
# summary results
library("xtable")
library("digest")
rm(list=ls())

# load functions
source("https://raw.githubusercontent.com/boopsboops/reference-libraries/master/scripts/funs.R")

# load reference lib
source("https://raw.githubusercontent.com/boopsboops/reference-libraries/master/scripts/references-load.R")
# loads objects: uk.species.table, uk.species.table.orig, uk.species.table.common, reflib.orig, 


###### reference database stats #######
###### reference database stats #######
###### reference database stats #######

# load up the ID rates data
id.results <- read_csv(file="../data/in-silico-results.csv")
# all names 
uk.species.table.orig %>% summarise(n=n()) %>% pull(n) # 4897 names
# all valid
uk.species.table %>% summarise(n=n()) %>% pull(n) # 531 valid
# common and valid
uk.species.table.common %>% summarise(n=n()) %>% pull(n) # 176 common valid

# numbers of species/seqs in reflib
reflib.orig %>% summarise(n=n()) %>% pull(n)# 43,366 sequences for all
reflib.orig %>% distinct(sciNameValid) %>% summarise(n=n()) %>% pull(n) # 491 species total for all markers
reflib.orig %>% filter(sciNameValid %in% uk.species.table.common$sciName) %>% summarise(n=n()) %>% pull(n) # 25,799 sequences for common 
reflib.orig %>% filter(sciNameValid %in% uk.species.table.common$sciName) %>% distinct(sciNameValid) %>% summarise(n=n()) %>% pull(n) # 172 common spp 
common.got <- reflib.orig %>% filter(sciNameValid %in% uk.species.table.common$sciName) %>% distinct(sciNameValid) %>% pull(sciNameValid)
setdiff(uk.species.table.common$sciName, common.got)# those common species with no data

# add the reflib proportions
id.results %<>% mutate(nSpp=if_else(subset=="all",(nSpp/531),(nSpp/176))) %>% rename(libProp=nSpp)

# print
id.results %>% arrange(subset,locus) %>%
    mutate(libProp=round(libProp,digits=2),meanPerSpp=round(meanPerSpp,digits=1),idRate=round(idRate,digits=2),idRateShared=round(idRateShared,digits=2),primerBiasActin=round(primerBiasActin,digits=1),primerBiasElasmo=round(primerBiasElasmo,digits=1)) %>% 
    mutate_all(.funs=function(x){prettyNum(x,big.mark=",")}) %>% 
    mutate(meanPerSpp=paste0(meanPerSpp," (",medianPerSpp,")")) %>% 
    mutate(idRate=paste0(idRate," (",idRateShared,")")) %>% 
    select(locus,dataset,subset,totalSeqs,libProp,meanPerSpp,idRate,primerBiasActin,primerBiasElasmo) %>% 
    xtable(caption="blahhh", digits=c(0,0,0,0,0,0,0,0,0,0)) %>% 
    print.xtable(include.rownames=FALSE,booktabs=TRUE,sanitize.text.function=identity,caption.placement="top",size="scriptsize")



###### primers table #######
###### primers table #######
###### primers table #######
# read in
primers <- read_csv(file="../data/primers.csv")

# print out
primers %>% mutate(name=paste(reference,primer,sep=".")) %>% 
    select(name,locus,officialName,oligo,lengthFrag,refTex) %>%
    xtable(caption="blahhh", digits=c(0,0,0,0,0,0,0)) %>% 
    print.xtable(include.rownames=FALSE, sanitize.text.function=identity)


###### sequencing summary #######
###### sequencing summary #######
###### sequencing summary #######

# add these results into '../data/sequencing-summary.csv'

# to get the rawReads, aquatheatReads, primerReads, barcodeReads, trimReads, filterReads, mergeReads, chimReads, homolReads,
# go to '../temp/fastq/$MARKER/logs'

# to get the minPerSampleReads, maxPerSampleReads, meanPerSampleReads, sdPerSampleReads
# go to 'make-plots.R' and run from ### PLOT OF READ DEPTHS PER SAMPLE 

# ASV global blast taxonomy results
# load
blast.results.all <- read_csv(file="../data/blast-results-all.csv")
# check total reads
blast.results.all %>% group_by(marker) %>% summarise(reads=sum(totalReads))
# by superkingdom
blast.results.all %>% group_by(marker,sskingdoms) %>% summarise(reads=sum(totalReads)) %>% mutate(propTot=reads/sum(reads)) %>% arrange(sskingdoms,marker)
# by kingdom
blast.results.all %>% group_by(marker,kingdom) %>% summarise(reads=sum(totalReads)) %>% mutate(propTot=reads/sum(reads)) %>% arrange(kingdom,marker)
# by phylum
blast.results.all %>% group_by(marker,phylum) %>% summarise(reads=sum(totalReads)) %>% mutate(propTot=reads/sum(reads)) %>% arrange(phylum,marker) %>% print(n=Inf)
# by fish
data(fishbase)
blast.results.all %>% 
    mutate(isFish=if_else(kingdom=="Metazoa" & class=="Actinopteri" | kingdom=="Metazoa" & class=="Cephalaspidomorphi" | kingdom=="Metazoa" & class=="Chondrichthyes" | kingdom=="Metazoa" & order %in% fishbase$Order | kingdom=="Metazoa" & family %in% fishbase$Family | kingdom=="Metazoa" & genus %in% fishbase$Genus,"Fish", "Not fish",missing="Not fish")) %>% 
    group_by(marker,isFish) %>% 
    summarise(reads=sum(totalReads)) %>% 
    mutate(propTot=reads/sum(reads)) %>% 
    arrange(isFish,marker)

#blast.results.all %>% 
#    mutate(isFish=if_else(class=="Actinopteri" | class=="Cephalaspidomorphi" | class=="Chondrichthyes" | order %in% fishbase$Order | family %in% fishbase$Family | genus %in% fishbase$Genus,"Fish","Not fish",missing="Not fish")) %>% 
#    filter(marker=="leray-xt",isFish=="Fish") %>% print(n=Inf)

# ADD the fish assignments
tab <- read_csv(file="../data/fish-assignments.csv")

# ASVs, swarm otus and reads
tab %>% group_by(marker) %>% summarise(ASVs=n_distinct(qseqid),swarmOTUs=n_distinct(swarm),reads=sum(totalReads))
# assigned
tab %>% filter(assigned==TRUE) %>% group_by(marker) %>% summarise(species=n_distinct(sciNameValidLocal),reads=sum(totalReads))
# unassigned fish reads
tab %>% filter(assigned==FALSE) %>% group_by(marker) %>% summarise(unassignedASVs=n_distinct(qseqid),swarmOTUs=n_distinct(swarm),reads=sum(totalReads))
# unassigned OTHER reads
#tab %>% filter(is.na(assigned)) %>% group_by(marker) %>% summarise(unassignedASVs=n_distinct(qseqid),reads=sum(totalReads))
# open and save into "../data/sequencing-summary.csv"

# to PRINT - load up sequencing summary
reads.results <- read_csv(file="../data/sequencing-summary.csv")

# substract the aquatheat reads and remove unwanted cols
# then rename for the table
reads.results.reform <- reads.results %>% mutate(totalReads=totalReads-aquatheatReads,primerReads=primerReads-aquatheatReads,barcodeReads=barcodeReads-aquatheatReads,trimReads=trimReads-aquatheatReads) %>%
    select(-aquatheatReads,-minPerSampleReads,-maxPerSampleReads,-fishASVs,-fishSwarms,-fishSpeciesID,-fishUnassignedASVs,-fishUnassignedSwarms) %>% 
    rename(Total=totalReads,`With primer`=primerReads,`With barcode`=barcodeReads,Trimmed=trimReads,
        `Filtered`=filterReads,Merged=mergeReads,`Cleaned chimaeric`=chimReads,`Homologous`=homolReads,
        `Mean per replicate`=meanPerSampleReads,`(SD)`=sdPerSampleReads,
        Bacteria=bacteriaReads,Eukaryota=eukaryoteReads,Metazoa=metazoaReads,Chordata=chordateReads,
        `Fish (Putative)`=fishReads,`Fish (assigned)`=fishSpeciesReads,`Fish (unassigned)`=fishUnassignedReads)

# print out in xtable after transposing
reads.results.reform %>%
    gather(key=Statistic, value=nReads, 2:ncol(reads.results.reform), factor_key=TRUE) %>%
    spread(key="marker", value=nReads) %>%
    select(Statistic,`leray-xt`,`sea-mid`,`sea-short`,miya) %>%
    mutate_all(.funs=function(x){prettyNum(x,big.mark=",")}) %>% 
    xtable(caption="blahhh", digits=c(0,0,0,0,0,0)) %>% 
    print.xtable(include.rownames=FALSE,booktabs=TRUE,sanitize.text.function=identity,caption.placement="top",size="scriptsize")


# make a quick perc summary of the global blast results
reads.results %>% 
    group_by(marker) %>% 
    summarise(bact=bacteriaReads/homolReads, bactExclChords=bacteriaReads/(homolReads-chordateReads), euk=eukaryoteReads/homolReads, eukExclChords=(eukaryoteReads-chordateReads)/(homolReads-chordateReads),
        meta=metazoaReads/homolReads, metaExclChords=(metazoaReads-chordateReads)/(homolReads-chordateReads), chord=chordateReads/homolReads, fish=fishReads/homolReads) %>%
    group_by(marker) %>%
    summarise_all(function(x){x*100})

# make a quick perc summary of the stringent ids
reads.results %>% 
    group_by(marker) %>% 
    summarise(assigned=fishSpeciesReads/fishReads,unassigned=fishUnassignedReads/fishReads) %>%
    group_by(marker) %>%
    summarise_all(function(x){x*100})

# summary of otus and asvs
reads.results %>% 
    group_by(marker) %>% 
    summarise(fishSpeciesID=fishSpeciesID,fishASVs=fishASVs,fishUnassignedASVs=fishUnassignedASVs,fishSwarms=fishSwarms,fishUnassignedSwarms=fishUnassignedSwarms)


## get ratios of the number species per marker  
results.all <- read_csv(file="../data/results-by-marker.csv")

results.all %>% group_by(marker,eventID) %>% 
    summarise(spp=length(unique(species))) %>% 
    arrange(eventID,desc(spp)) %>% 
    group_by(eventID) %>% 
    summarise(prop=first(spp)/nth(spp,2))


## supporting table - unassigned Miya reads ##
## supporting table - unassigned Miya reads ##
## supporting table - unassigned Miya reads ##


# set marker and proj dir 
marker <- "miya"
marker <- "leray-xt"
marker <- "sea-short"
marker <- "sea-mid"
#
rel.path <- "../temp/fastq"
proj.path <- paste(rel.path,marker,sep="/") 

# load
taxonomy.results <- read_csv(file="../data/fish-assignments.csv", col_names=TRUE)
plates <- read_csv(file="../data/sample-plates.csv")
otu.tab <- read.table(file=paste0(proj.path,"/results/asv-table-clean.tsv"), sep="\t", header=TRUE, as.is=TRUE, row.names=1, check.names=FALSE)
curated.table.fish <- as_tibble(otu.tab, rownames="asv")

# subset unassigned
taxonomy.results.unassigned <- taxonomy.results %>% filter(marker==as.character(!!marker) & assigned==FALSE) #%>% select(marker) %>% print(n=Inf)

# get plates
plates %<>% filter(primerSet==marker) %>% mutate(eventMD5=str_trunc(sapply(eventID, digest, algo="md5", serialize=FALSE, USE.NAMES=FALSE), width=12, side="right", ellipsis=""))

# filter process the table 
unassigned.by.location <- curated.table.fish %>% 
    filter(asv %in% taxonomy.results.unassigned$qseqid) %>% 
    gather(key=sample,value=nreads,-asv) %>% filter(nreads>0) %>% 
    mutate(hash=str_split_fixed(sample,"\\.",3)[,1]) %>%
    mutate(eventID=plates$eventID[match(hash,plates$eventMD5)]) %>% 
    group_by(asv) %>% 
    summarise(location=paste(unique(eventID),collapse=", ")) %>% 
    rename(qseqid=asv) %>%
    arrange(qseqid) #%>% print(n=Inf)

# reduce
unassigned.reduced <- taxonomy.results.unassigned %>% left_join(unassigned.by.location, by="qseqid") %>% select(qseqid,sscinames,pident,nident,totalReads,location,swarm,epaID)

# process and tidy
unassigned.reduced %<>% 
    group_by(swarm) %>% 
    mutate(blastGlobalIDs=paste(unique(sscinames),collapse=", "), allLocation=paste(unique(location),collapse=", "), nASVs=length(unique(qseqid)), allEPAs=paste(unique(epaID),collapse=", "), swarmReads=sum(totalReads),probable=NA,comment=NA) %>% 
    ungroup() %>%
    select(swarm,nASVs,swarmReads,blastGlobalIDs,allEPAs,probable,comment,allLocation) %>% 
    arrange(desc(swarmReads)) %>% 
    distinct() %>%
    mutate(allLocation=sapply(sapply(str_split(allLocation,", "),unique),paste,collapse=", ")) 

# print 
unassigned.reduced %>% 
    mutate(allLocation=str_replace_all(allLocation,"EA-2016.10.26-ESK-FYKE","Esk-Fyke"),
    allLocation=str_replace_all(allLocation,"EA-2016.10.26-ESK-SEI1","Esk-Seine"),
    allLocation=str_replace_all(allLocation,"EA-2016.10.31-TEES-SEI1","Tees"),
    allLocation=str_replace_all(allLocation,"MBA-2016.11.04-WHIT-02","Whitsand Bay"),
    allLocation=str_replace_all(allLocation,"PISC-2016.11.24-SOT-01","Test"),
    blastGlobalIDs=str_replace_all(blastGlobalIDs,"^","\\\\emph{"),
    blastGlobalIDs=str_replace_all(blastGlobalIDs,"$","}"),
    swarm=str_replace_all(swarm,"swarm","otu"),
        swarmReads=prettyNum(swarmReads,big.mark=",")) %>%
    rename(OTU=swarm,`Number ASVs`=nASVs,`Total reads`=swarmReads,`GenBank BLAST match`=blastGlobalIDs,`EPA identification`=allEPAs,`Probable species`=probable,Comment=comment,Locations=allLocation) %>% 
    xtable(caption="ff", digits=c(0,0,0,0,0,0,0,0,0)) %>% 
    print.xtable(include.rownames=FALSE,booktabs=TRUE,sanitize.text.function=identity)

# unassigned
taxonomy.results %>% filter(marker==as.character(!!marker) & assigned==FALSE) %>% summarise(tot=sum(totalReads))
# get number of Hippoglossoides and Sygnathus
taxonomy.results %>% filter(marker==as.character(!!marker) & assigned==TRUE & sciNameValidLocal=="Hippoglossoides platessoides") %>% summarise(tot=sum(totalReads))
taxonomy.results %>% filter(marker==as.character(!!marker) & assigned==TRUE & sciNameValidLocal=="Syngnathus typhle") %>% summarise(tot=sum(totalReads))


## supporting table - all assigned reads and trad surveys

# load the data
edna <- read_csv(file="../data/results-by-marker.csv")
trad <- read_csv(file="../data/traditional-survey.csv")

# format trad
trad %<>% rename(nreads=individuals) %>% mutate(marker="trad")

# join frames
edna.trad <- full_join(edna,trad)

# clean data
edna.trad %<>% group_by(species,marker,eventID) %>% 
    summarise(nreads=sum(nreads)) %>% 
    ungroup() %>% 
    mutate(nreads=prettyNum(nreads,big.mark=",")) %>%
    spread(key=marker,value=nreads) %>% 
    mutate(eventID=as_factor(eventID)) %>% 
    select(eventID,species,trad,miya,`leray-xt`,`sea-mid`,`sea-short`) %>%
    rename(Species=species,Traditional=trad,`12S MiFish-U`=miya,`COI Leray-XT`=`leray-xt`,`COI SeaDNA-mid`=`sea-mid`,`COI SeaDNA-short`=`sea-short`)
print(edna.trad)

# split dataframe
edna.trad.list <- split(edna.trad, edna.trad$eventID)

# reorder
edna.trad.list <- edna.trad.list[c(2,4,5,3,1)]
names(edna.trad.list)

# function to create tables
tex_tab <- function(df,name){
df %<>% select(-eventID) %>% mutate(Species=str_replace_all(Species,"^","\\\\emph{"),Species=str_replace_all(Species,"$","}"))
x.tab <- xtable(df, caption=paste("Metabarcoding and traditional fish survey results for the",name,"site survey. Values correspond to the number of reads identified to species for the molecular markers, and the number of individuals caught on the traditional surveys. Species separated by semicolon are those for which matches were ambiguous."), 
    digits=0, align=c("l","l","r","r","r","r","r"))
print.xtable(x.tab,include.rownames=FALSE,booktabs=TRUE,sanitize.text.function=identity,caption.placement="top",size="scriptsize")
}

# print the tables
mapply(function(x,y) tex_tab(df=x,name=y), edna.trad.list, names(edna.trad.list))


## eDNA congruence
# load the data
edna <- read_csv(file="../data/results-by-marker.csv")

# collapse reps
edna %<>% group_by(marker,eventID,species) %>% summarise_at(vars(nreads),sum) %>% ungroup() %>% arrange(marker,eventID,desc(nreads)) 

# make a new line for ambig spp
edna %>% filter(grepl("; ",species)) %>% group_by(eventID) %>% distinct(species, .keep_all=TRUE) %>% arrange(eventID)
edna %<>% mutate(species2=str_split_fixed(species,"; ",2)[,2], species=str_split_fixed(species,"; ",2)[,1])

# this breaks the ambiguous spp into two rows and rejoins
sp2 <- edna %>% filter(species2!="") %>% select(-species) %>% rename(species=species2)
sp1 <- edna %>% select(-species2)
edna <- bind_rows(sp2,sp1)

# remove dups
edna %<>% group_by(species,eventID) %>% select(-nreads) %>% distinct()

# eDNA congrence 
edna %>% group_by(species,eventID) %>% summarise(timesDetected=n()) %>% arrange(desc(timesDetected),eventID) %>% print(n=Inf)
edna %>% group_by(species,eventID) %>% 
    summarise(timesDetected=n()) %>% 
    group_by(eventID) %>% 
    count(timesDetected) %>% 
    mutate(totSpp=sum(n),prop=n/totSpp) %>%
    arrange(timesDetected)
