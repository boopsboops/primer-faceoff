#!/usr/bin/env Rscript

# load libs
library("tidyverse")
library("magrittr")
library("lubridate")
library("digest")
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/tab2fas.R")

# set the name of the folder that the file will be added to 
marker <- "miya"
trlens <- c(21,27)
#
marker <- "leray-xt"
trlens <- c(26,26)
#
marker <- "sea-short"
trlens <- c(20,20)
#
marker <- "sea-mid"
trlens <- c(20,21)

# set proj
project <- "primer-faceoff"

# load up the data
plates <- read_csv(file="../data/sample-plates.csv")

# filter by marker
plates %<>% filter(primerSet==marker)

# create the barcode tags
plates %<>% mutate(barcodesFwd=str_replace_all(oligoFwd,"N",""), 
    barcodesFwd=str_trunc(barcodesFwd, width=10, side="right", ellipsis=""),
    barcodesRev=str_replace_all(oligoRev,"N",""),
    barcodesRev=str_trunc(barcodesRev, width=10, side="right", ellipsis=""),
    primerFwd=str_trunc(oligoFwd, width=trlens[1], side="left", ellipsis=""),
    primerRev=str_trunc(oligoRev, width=trlens[2], side="left", ellipsis=""),
    labelFwd=str_trunc(oligoFwd, width=unique(str_length(oligoFwd)-trlens[1]), side="right", ellipsis=""),
    labelRev=str_trunc(oligoRev, width=unique(str_length(oligoRev)-trlens[2]), side="right", ellipsis="")
    )

# create a unique md5 hash for each sample
plates %<>% mutate(eventMD5=str_trunc(sapply(eventID, digest, algo="md5", serialize=FALSE, USE.NAMES=FALSE), width=12, side="right", ellipsis=""))

# create the labels
plates %<>% mutate(senseLabel=paste0(labelFwd, ".", eventMD5, ".", rep, ".", pcr), antisenseLabel=paste0(labelRev, ".", eventMD5, ".", rep, ".", pcr)) 


# make fasta 
bcF <- tab2fas(df=plates, seqcol="barcodesFwd", namecol="senseLabel")
bcR <- tab2fas(df=plates, seqcol="barcodesRev", namecol="antisenseLabel")

# write out new dirs if needed to run first time
# dir.create(path=paste0("../temp/fastq/", marker),recursive=TRUE)

# write out in folder for running analyses
write.FASTA(bcF, file=paste0("../temp/fastq/", marker, "/barcodes-sense.fas"))
write.FASTA(bcR, file=paste0("../temp/fastq/", marker,"/barcodes-antisense.fas"))
