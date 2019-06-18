#!/usr/bin/env Rscript
rm(list=ls())
# load libs
library("tidyverse")
library("magrittr")
library("lubridate")
library("stringr")
library("ape")
require("stringdist")
require("phangorn")
library("dada2")# used dada2 1.10.1
library("parallel")
library("xtable")

# load funs
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/tab2fas.R")
source("https://raw.githubusercontent.com/boopsboops/reference-libraries/master/scripts/funs.R")
source("https://raw.githubusercontent.com/boopsboops/reference-libraries/master/scripts/references-load.R")
reflib <- reflib.orig


# set marker and proj dir 
marker <- "miya"
marker <- "leray-xt"
marker <- "sea-short"
marker <- "sea-mid"
#
rel.path <- "../temp/fastq"
proj.path <- paste(rel.path,marker,sep="/") 

# file paths
path.sense <- paste0(proj.path,"/sense/trimmed")
path.antisense <- paste0(proj.path,"/antisense/trimmed")
files.sense.R1 <- sort(list.files(path.sense,".R1.fastq.gz", full.names=TRUE))
files.sense.R2 <- sort(list.files(path.sense,".R2.fastq.gz", full.names=TRUE))
files.antisense.R1 <- sort(list.files(path.antisense,".R1.fastq.gz", full.names=TRUE))
files.antisense.R2 <- sort(list.files(path.antisense,".R2.fastq.gz", full.names=TRUE))

# remove the aquatheatr sequences
files.sense.R1 <- files.sense.R1[grep("AquaTheat",files.sense.R1,invert=TRUE)]
files.sense.R2 <- files.sense.R2[grep("AquaTheat",files.sense.R2,invert=TRUE)]
files.antisense.R1 <- files.antisense.R1[grep("AquaTheat",files.antisense.R1,invert=TRUE)]
files.antisense.R2 <- files.antisense.R2[grep("AquaTheat",files.antisense.R2,invert=TRUE)]


# quality plots
plotQualityProfile(files.sense.R1[1:9])# just 9
plotQualityProfile(files.sense.R2[1:9])# just 9
plotQualityProfile(files.antisense.R1[1:9])# just 9
plotQualityProfile(files.antisense.R2[1:9])# just 9


# make file paths for filtered fastq
files.sense.filt.R1 <- str_replace_all(files.sense.R1, "trimmed", "filtered")
files.sense.filt.R2 <- str_replace_all(files.sense.R2, "trimmed", "filtered")
files.antisense.filt.R1 <- str_replace_all(files.antisense.R1, "trimmed", "filtered")
files.antisense.filt.R2 <- str_replace_all(files.antisense.R2, "trimmed", "filtered")


# trucLens
trucVal <- c(105,105)# MIYA - truncLen 105 gives at least 29 bp overlap for the longest amplicons (e.g. Raja clavata @ 181 bp), and 40 bp for the regular 170 bp
trucVal <- c(170,170)# Leray-XT - truncLen 170 gives about 28 bp overlap
trucVal <- c(40,40)# SEA-SHORT - truncLen 40 gives about 25 bp overlap
trucVal <- c(80,80)# SEA-MID - truncLen 80 gives about 30 bp overlap

# quality trim Ns and truncate
filterAndTrim(fwd=files.sense.R1, filt=files.sense.filt.R1, rev=files.sense.R2, filt.rev=files.sense.filt.R2, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, truncLen=trucVal, multithread=TRUE, verbose=TRUE)
filterAndTrim(fwd=files.antisense.R1, filt=files.antisense.filt.R1, rev=files.antisense.R2, filt.rev=files.antisense.filt.R2, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, truncLen=trucVal, multithread=TRUE, verbose=TRUE)

# learn errors (can restart from here)
set.seed(42)
files.sense.filt.R1.errs <- learnErrors(files.sense.filt.R1, multithread=TRUE, randomize=TRUE, nbases=1e+08)
set.seed(42)
files.sense.filt.R2.errs <- learnErrors(files.sense.filt.R2, multithread=TRUE, randomize=TRUE, nbases=1e+08)
set.seed(42)
files.antisense.filt.R1.errs <- learnErrors(files.antisense.filt.R1, multithread=TRUE, randomize=TRUE, nbases=1e+08)
set.seed(42)
files.antisense.filt.R2.errs <- learnErrors(files.antisense.filt.R2, multithread=TRUE, randomize=TRUE, nbases=1e+08)

# plot the errors
#plotErrors(dq=files.sense.filt.R1.errs, nominalQ=TRUE)
#plotErrors(dq=files.sense.filt.R2.errs, nominalQ=TRUE)
#plotErrors(dq=files.antisense.filt.R1.errs, nominalQ=TRUE)
#plotErrors(dq=files.antisense.filt.R2.errs, nominalQ=TRUE)

# derep the sequences
files.sense.filt.R1.derep <- derepFastq(files.sense.filt.R1, verbose=TRUE)
files.sense.filt.R2.derep <- derepFastq(files.sense.filt.R2, verbose=TRUE)
files.antisense.filt.R1.derep <- derepFastq(files.antisense.filt.R1, verbose=TRUE)
files.antisense.filt.R2.derep <- derepFastq(files.antisense.filt.R2, verbose=TRUE)

# make sample names
names(files.sense.filt.R1.derep) <- str_replace_all(names(files.sense.filt.R1.derep), ".fastq.gz", "")
names(files.sense.filt.R2.derep) <- str_replace_all(names(files.sense.filt.R2.derep), ".fastq.gz", "")
names(files.antisense.filt.R1.derep) <- str_replace_all(names(files.antisense.filt.R1.derep), ".fastq.gz", "")
names(files.antisense.filt.R2.derep) <- str_replace_all(names(files.antisense.filt.R2.derep), ".fastq.gz", "")

# make some fish priors
prefix <- "coi.lerayxt.noprimers"
prefix <- "coi.seamid.noprimers"
prefix <- "coi.seashort.noprimers"
prefix <- "12s.miya.noprimers"

# subset the marker from the reflib - using 50% length cutoff
reflib.sub <- subset_by_marker(prefix=prefix,df=reflib,thresh=0.5)
# collapse haps
reflib.red <- bind_rows(mcmapply(FUN=function(x) hap_collapse_df(df=x,lengthcol=paste0("lengthFrag.",prefix),nuccol=paste0("nucleotidesFrag.",prefix)), split(reflib.sub,reflib.sub$sciNameValid), SIMPLIFY=FALSE, mc.cores=8))
# pull out seqs
fish.priors <- reflib.red %>% pull(paste0("nucleotidesFrag.",prefix))

# fun to revcomp
library("seqinr")
flip <- function(x){
    revcomp <- c2s(rev(comp(s2c(x),ambiguous=TRUE)))
    return(revcomp)}

# get revcomps and trim
fish.priors.revcomp <- mapply(flip,fish.priors,USE.NAMES=FALSE)
fish.priors <- unique(str_trunc(fish.priors, width=trucVal[1], side="right", ellipsis=""))
fish.priors.revcomp <- unique(str_trunc(fish.priors.revcomp, width=trucVal[1], side="right", ellipsis=""))


# run dada denoising
files.sense.filt.R1.dada <- dada(files.sense.filt.R1.derep, err=files.sense.filt.R1.errs, multithread=TRUE, pool=FALSE, priors=fish.priors)
files.sense.filt.R2.dada <- dada(files.sense.filt.R2.derep, err=files.sense.filt.R2.errs, multithread=TRUE, pool=FALSE, priors=fish.priors.revcomp)
files.antisense.filt.R1.dada <- dada(files.antisense.filt.R1.derep, err=files.antisense.filt.R1.errs, multithread=TRUE, pool=FALSE, priors=fish.priors.revcomp)
files.antisense.filt.R2.dada <- dada(files.antisense.filt.R2.derep, err=files.antisense.filt.R2.errs, multithread=TRUE, pool=FALSE, priors=fish.priors)

# look at dada objects
files.sense.filt.R1.dada[[1]]
files.sense.filt.R2.dada[[1]]
files.antisense.filt.R1.dada[[1]]
files.antisense.filt.R2.dada[[1]]

# merge the R1 and R2
sense.merged <- mergePairs(files.sense.filt.R1.dada, files.sense.filt.R1.derep, files.sense.filt.R2.dada, files.sense.filt.R2.derep, verbose=TRUE)
antisense.merged <- mergePairs(files.antisense.filt.R1.dada, files.antisense.filt.R1.derep, files.antisense.filt.R2.dada, files.antisense.filt.R2.derep, verbose=TRUE)

head(sense.merged[[1]])
head(antisense.merged[[1]])

# make an OTU table
sense.seqtab <- makeSequenceTable(sense.merged)
antisense.seqtab <- makeSequenceTable(antisense.merged)
dim(sense.seqtab)
dim(antisense.seqtab)

# reverse comp the antisense
colnames(antisense.seqtab) <- dada2:::rc(colnames(antisense.seqtab))

# fix the names before merging
rownames(sense.seqtab) <- apply(str_split_fixed(rownames(sense.seqtab),"\\.",5)[,2:4],1,paste,collapse=".")
rownames(antisense.seqtab) <- apply(str_split_fixed(rownames(antisense.seqtab),"\\.",5)[,2:4],1,paste,collapse=".")

# merge the tables
merged.seqtab <- mergeSequenceTables(table1=sense.seqtab, table2=antisense.seqtab, repeats="sum")
dim(merged.seqtab)

# remove chimaeras
merged.seqtab.nochim <- removeBimeraDenovo(merged.seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(merged.seqtab.nochim)

# stats to copy into 'logs/stats.csv'
table(nchar(getSequences(merged.seqtab.nochim)))
sum(merged.seqtab)
sum(merged.seqtab.nochim)
sum(merged.seqtab.nochim)/sum(merged.seqtab)

# make df and fasta for IDs
otus.df <- tibble(names=paste0("asv",seq_along(colnames(merged.seqtab.nochim))), dnas=colnames(merged.seqtab.nochim)) %>% mutate(len=str_length(dnas))

# write out
write.FASTA(tab2fas(df=otus.df, seqcol="dnas", namecol="names"), file=paste0(proj.path,"/results/asvs.fna"))

# save the OTU table as df
colnames(merged.seqtab.nochim) <- paste0("asv",seq_along(colnames(merged.seqtab.nochim)))
write_tsv(as_tibble(t(merged.seqtab.nochim), rownames="asv"), path=paste0(proj.path,"/results/asv-table.tsv"))
