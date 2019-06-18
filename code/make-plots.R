#!/usr/bin/env Rscript

# load libs
library("tidyverse")
library("usedist")
library("RColorBrewer")
library("magrittr")
library("VennDiagram")
library("ggthemes")
library("digest")

# load data
results.all <- read_csv(file="../data/results-by-marker.csv")
dir.create("../temp/figs")


# plot boxplot of num species per marker and site event
p1 <- results.all %>% group_by(marker,eventID) %>% 
    mutate(totSpp=length(unique(species))) %>% 
    group_by(rep,add=TRUE) %>%
    summarise(spp=length(unique(species)),reads=sum(nreads),totSpp=unique(totSpp)) %>% 
    ungroup() %>% 
    mutate(marker=str_replace_all(marker,"miya","Miya-12S"),
        marker=str_replace_all(marker,"sea-short","COI-fish-short"),
        marker=str_replace_all(marker,"sea-mid","COI-fish-mid"),
        marker=str_replace_all(marker,"leray-xt","COI-Leray-XT"),
        eventID=str_replace_all(eventID,"EA-2016.10.26-ESK-FYKE","River Esk (fyke)"),
        eventID=str_replace_all(eventID,"EA-2016.10.26-ESK-SEI1","River Esk (seine)"),
        eventID=str_replace_all(eventID,"EA-2016.10.31-TEES-SEI1","River Tees"),
        eventID=str_replace_all(eventID,"MBA-2016.11.04-WHIT-02","Whitsand Bay"),
        eventID=str_replace_all(eventID,"PISC-2016.11.24-SOT-01","River Test")) %>% 
    mutate(marker=fct_relevel(marker,"Miya-12S","COI-fish-short","COI-fish-mid","COI-Leray-XT"),
        event=fct_relevel(eventID, "River Tees","River Esk (fyke)","River Esk (seine)","Whitsand Bay","River Test")) %>% 
    ggplot(aes(x=marker,y=spp)) + 
        geom_boxplot(aes(fill=marker)) + 
        #geom_jitter(width=0.1,alpha=0.2) +
        geom_point(aes(x=marker,y=totSpp),shape=4,size=2,alpha=0.1) +
        facet_wrap(.~event) +
        theme_bw() + scale_color_ptol() + scale_fill_ptol() +
        theme(axis.text.x=element_blank(),
            axis.text.y=element_text(size=10),
            axis.title.x=element_text(size=14),
            axis.title.y=element_text(size=14),
            strip.text.x=element_text(size=11),
            strip.text.y=element_text(size=11),
            strip.background=element_blank(),
            #panel.grid.major=element_blank(),
            #panel.grid.minor=element_blank(),
            legend.text=element_text(size=10),
            plot.background=element_blank(),
            axis.ticks.x = element_blank(),
            legend.title=element_text(size=12)) +
        labs(fill="Primer set",x="Primer set",y="Number species detected")
plot(p1)
ggsave(filename="../temp/figs/marker-by-species.pdf",plot=p1,width=10,height=7)


##### Function to make a heatmap of PCR amplification sucess for top n species

# choose number of top n
top.val <- 10

# sort by top n total reads per marker-eventID combination
# exclude leray-xt
# filter(marker!="leray-xt") %>% 
top.n <- results.all %>% group_by(marker,eventID,species) %>% summarise_at(vars(nreads),sum) %>% ungroup() %>% arrange(marker,eventID,desc(nreads)) %>% group_by(marker,eventID) %>% slice(1:top.val) %>% ungroup() %>% arrange(marker,eventID,desc(nreads))
top.n %>% print(n=Inf)

# look at congruence with markers
ev.mar <- results.all %>% group_by(marker,eventID,species) %>% summarise_at(vars(nreads),sum) %>% ungroup() %>% mutate(species=str_split_fixed(species,"; ",2)[,1])
ev.mar %>% distinct(species) %>% arrange(species) %>% print(n=Inf)
ev.mar %<>% group_by(eventID,species) %>% count(species) %>% group_by(eventID) %>% count(n) %>% group_by(eventID) %>% mutate(tot=sum(nn), prop=nn/tot) 
ev.mar %>% print()
ev.mar %>% filter(n==3) %>% print()
ev.mar %>% filter(n==2) %>% print()
ev.mar %>% filter(n==1) %>% print()

# split by top n species by marker/event into a list
top.split <- split(top.n,list(top.n$marker,top.n$eventID))

# extract the species
top.spp <- lapply(top.split,function(x) x$species)

# again remove leray from full data
#results.sub <- results.all %>% filter(marker!="leray-xt")

# split the full data by marker-event and check names
results.split <- split(results.all ,list(results.all$marker,results.all$eventID))
names(top.spp)==names(results.split)

# fun to extract the top n species species from full data 
subset_species <- function(df,subsett){
    res <- df %>% filter(species %in% subsett)
    return(res)
}

# subset the full data by top n species
setted <- mapply(function(x,y) subset_species(x,y), results.split, top.spp, SIMPLIFY=FALSE)

# tabulate all the species occurances
setted.mat <- lapply(setted,function(x) as.matrix(with(x, table(rep,species))))

# take a look at transposed matrices
lapply(setted.mat, t)

# function to calculate number of shared species between two vectors (i.e. reps)
similarity <- function(x,y){
    res <- length(which(x==1 & y==1))/top.val
    return(res)
}

# run this fun over each matrix using custion distance matrix function
setted.sim <- lapply(setted.mat, function(x) as.matrix(dist_make(x, similarity)))

# create the diagonals
diag.vals <- lapply(setted.mat, function(k) apply(k,1,function(x){length(which(x==1))/top.val}))
names(setted.sim)==names(diag.vals)

# add the diagonals to the matrices
setted.sim.diags <- mapply(function(x,y){
  diag(x) <- y
  return(x)},
  setted.sim, diag.vals, SIMPLIFY=FALSE)
 
# convert to tibble 
sims.tibs <- lapply(setted.sim.diags, function(x) as_tibble(rownames_to_column(as.data.frame(x),var="rep1")))

# convert to long format
sims.tibs.g <- lapply(sims.tibs,function(x) gather(x,key=rep2,value=similarity,-rep1))

# add marker and event labels
sims.tibs.complete <- mapply(function(x,y,z) mutate(x,marker=y,event=z), sims.tibs.g, str_split_fixed(names(sims.tibs.g), "\\.", 2)[,1], str_split_fixed(names(sims.tibs.g), "\\.", 2)[,2],SIMPLIFY=FALSE)

# join into one dataframe
sims.df <- bind_rows(sims.tibs.complete)

# relab everything
sims.df %<>% mutate(rep1=str_replace_all(rep1,"\\.pcr",""),
    rep2=str_replace_all(rep2,"\\.pcr",""),
    marker=str_replace_all(marker,"miya","Miya-12S"),
    marker=str_replace_all(marker,"sea-short","COI-fish-short"),
    marker=str_replace_all(marker,"sea-mid","COI-fish-mid"),
    marker=str_replace_all(marker,"leray-xt","COI-Leray-XT"),
    event=str_replace_all(event,"EA-2016.10.26-ESK-FYKE","River Esk (fyke)"),
    event=str_replace_all(event,"EA-2016.10.26-ESK-SEI1","River Esk (seine)"),
    event=str_replace_all(event,"EA-2016.10.31-TEES-SEI1","River Tees"),
    event=str_replace_all(event,"MBA-2016.11.04-WHIT-02","Whitsand Bay"),
    event=str_replace_all(event,"PISC-2016.11.24-SOT-01","River Test"))

# relevel the factors
sims.df %<>% mutate(marker=fct_relevel(marker,"Miya-12S","COI-fish-short","COI-fish-mid","COI-Leray-XT"),
    event=fct_relevel(event, "River Tees","River Esk (fyke)","River Esk (seine)","Whitsand Bay","River Test"))

# get average repeatability
sims.df %>% filter(!rep1==rep2) %>% group_by(marker) %>% summarise(avg=mean(similarity))

#plot the heatmap
p2 <- ggplot(data=sims.df, aes(x=rep1, y=rep2)) + 
    geom_tile(aes(fill=similarity)) + 
    scale_fill_gradientn(colours=brewer.pal(n=9, name="Blues"), limits=c(0,1)) + 
    facet_grid(rows=vars(marker), cols=vars(event)) +
    labs(fill="Proportion of top n species amplified",x="Site location",y="Marker fragment") +
    theme_bw() +
    theme(axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        strip.text.x=element_text(size=12),
        strip.text.y=element_text(size=14),
        strip.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.text=element_text(size=12),
        plot.background=element_blank(),
        legend.title=element_text(size=12)) 
plot(p2)
# write out
ggsave(filename="../temp/figs/heatmap.pdf",plot=p2,width=15,height=10)


### PLOT OF READ DEPTHS PER SAMPLE

# load the plates
plates <- read_csv(file="../data/sample-plates.csv")
plates %<>% mutate(eventMD5=str_trunc(sapply(eventID, digest, algo="md5", serialize=FALSE, USE.NAMES=FALSE), width=12, side="right", ellipsis=""))

# list the markers and paths
markers.all <- c("miya","leray-xt","sea-short","sea-mid")
rel.path <- "../temp/fastq"

# function to load and rowsum all reads 
load_otus <- function(marker){
    path <- paste(rel.path,marker,"results","asv-table-clean.tsv",sep="/")
    tib <- as_tibble(t(read.table(file=path, sep="\t", header=TRUE, as.is=TRUE, row.names=1, check.names=FALSE)), rownames="rep")
    tib %<>% mutate(marker=marker)
    tib %<>% mutate_if(is.integer, as.numeric) %>% mutate(sum=pmap_dbl(select(.,-c(rep,marker)), sum)) %>% select(rep,sum,marker)
    return(tib)
}

# run the function
reps.by.marker <- bind_rows(mapply(load_otus,markers.all,SIMPLIFY=FALSE))

# break up the columns
reps.by.marker %<>% separate(rep,into=c("hash","rep","pcr"))

# get the eventID numbers
reps.by.marker %<>% mutate(eventID=plates$eventID[match(hash,plates$eventMD5)])

# get stats per marker
reps.by.marker %>% group_by(marker) %>% summarise(tot=sum(sum),min=min(sum),max=max(sum),mean=mean(sum),sd=sd(sum))

# plot boxplot
library(scales)
p1 <- reps.by.marker %>% 
    mutate(marker=str_replace_all(marker,"miya","Miya-12S"),
        marker=str_replace_all(marker,"sea-short","COI-fish-short"),
        marker=str_replace_all(marker,"sea-mid","COI-fish-mid"),
        marker=str_replace_all(marker,"leray-xt","COI-Leray-XT"),
        eventID=str_replace_all(eventID,"EA-2016.10.26-ESK-FYKE","River Esk (fyke)"),
        eventID=str_replace_all(eventID,"EA-2016.10.26-ESK-SEI1","River Esk (seine)"),
        eventID=str_replace_all(eventID,"EA-2016.10.31-TEES-SEI1","River Tees"),
        eventID=str_replace_all(eventID,"MBA-2016.11.04-WHIT-02","Whitsand Bay"),
        eventID=str_replace_all(eventID,"PISC-2016.11.24-SOT-01","River Test")) %>% 
    mutate(marker=fct_relevel(marker,"Miya-12S","COI-fish-short","COI-fish-mid","COI-Leray-XT"),
        event=fct_relevel(eventID, "River Tees","River Esk (fyke)","River Esk (seine)","Whitsand Bay","River Test")) %>% 
    group_by(marker,eventID) %>% 
    ggplot(aes(eventID,sum,fill=eventID)) + 
    #geom_bar(stat="identity") + 
    geom_boxplot() + 
    facet_wrap(.~marker, scales="free") +# or "fixed"
    scale_y_continuous(labels = comma) +
    theme_bw() + scale_color_ptol() + scale_fill_ptol() +
            theme(axis.text.x=element_blank(),
            axis.text.y=element_text(size=10),
            axis.title.x=element_text(size=14),
            axis.title.y=element_text(size=14),
            strip.text.x=element_text(size=11),
            strip.text.y=element_text(size=11),
            strip.background=element_blank(),
            #panel.grid.major=element_blank(),
            #panel.grid.minor=element_blank(),
            legend.text=element_text(size=10),
            plot.background=element_blank(),
            axis.ticks.x = element_blank(),
            legend.title=element_text(size=12)) +
        labs(fill="Location",x="Location",y="Total number reads")
plot(p1)
ggsave(filename="../temp/figs/marker-by-reads.pdf",plot=p1,width=10,height=7)


#### Compare the traditional surveys (Venn diagram)
library("cowplot")

# reload data
results.all <- read_csv(file="../data/results-by-marker.csv")

# load survey data
trad.all <- read_csv(file="../data/traditional-survey.csv")

# collapse by rep
results.all %<>% group_by(marker,eventID,species) %>% summarise_at(vars(nreads),sum) %>% ungroup() %>% arrange(marker,eventID,desc(nreads)) 

# check for ambiguous IDs and get names
results.all %>% filter(grepl("; ",species)) %>% group_by(eventID) %>% distinct(species, .keep_all=TRUE) %>% arrange(eventID)
dups <- results.all %>% filter(grepl("; ",species)) %>% group_by(eventID) %>% distinct(species) %>% pull(species) %>% str_split("; ") %>% unlist() %>% unique()

# create a new col for species 2
results.all %<>% mutate(species2=str_split_fixed(species,"; ",2)[,2], species=str_split_fixed(species,"; ",2)[,1])

# load uk species 
source("https://raw.githubusercontent.com/boopsboops/reference-libraries/master/scripts/references-load.R")

# check for spelling errors etc
setdiff(unique(trad.all$species),uk.species.table.common$sciName)

# add family names 
results.all %<>% mutate(family=uk.species.table.common$family[match(species,uk.species.table.common$sciName)])

# look at fams to make list
unique(results.all$family)

# fw species
spp.fw <- c(
    "Cottus gobio",
    "Gasterosteus aculeatus",
    "Cyprinidae",
    "Petromyzontidae",
    "Nemacheilidae",
    "Percidae",
    "Salmonidae")

# remove fw species from results
results.all %<>% filter(!family %in% spp.fw) %>% filter(!species %in% spp.fw)

# going to assume that for the duplicate spp., that if not present in trad survey then only the common one was found
# if both found, then also both in the mol survey

setdiff(dups,trad.all$species)
results.all %<>% mutate(species=str_replace_all(species, paste(setdiff(dups,trad.all$species),collapse="|"), ""), 
    species2=str_replace_all(species2, paste(setdiff(dups,trad.all$species),collapse="|"), ""))


# trad species not in each reflib
reflib <- reflib.orig
prefix <- "coi.lerayxt.noprimers"
prefix <- "coi.seamid.noprimers"
prefix <- "coi.seashort.noprimers"
prefix <- "12s.miya.noprimers"
reflib.sub <- reflib %>% filter(!is.na(!!as.name(paste0("nucleotidesFrag.",prefix))))
setdiff(unique(trad.all$species),unique(reflib.sub$sciNameValid))
(length(setdiff(unique(trad.all$species),unique(reflib.sub$sciNameValid))))/(length(unique(trad.all$species)))


#### run fun

# split the table into marker-event
results.all.split <- split(results.all,list(results.all$marker,results.all$eventID))

# get a name vector
results.all.split.names <- str_split_fixed(names(results.all.split),"\\.",2)[,2]

# get trad species
trad.all.split <- lapply(results.all.split.names, function(x) trad.all$species[trad.all$eventID==x])

# get lengths and intersection of lists
inters.three <- mapply(function(x,y) c(length(unique(c(x$species,x$species2))), length(unique(y)), length(intersect(c(x$species,x$species2),y))), results.all.split, trad.all.split, SIMPLIFY=FALSE)

# reorder for plot
inters.three <- inters.three[c(10,2,6,14,18,12,4,8,16,20,11,3,7,15,19,9,1,5,13,17)]

# calculate means
inters.list <- mapply(function(x) x[3]/x[2], inters.three, USE.NAMES=TRUE)
inters.df <- gather(map_dfr(inters.list, `[`), key="markerLocation", value="mean") %>% mutate(marker=str_split_fixed(markerLocation,"\\.",2)[,1], location=str_split_fixed(markerLocation,"\\.",2)[,2]) %>% select(-markerLocation)
inters.df %>% arrange(marker,mean) %>% print()
inters.df %>% group_by(marker) %>% summarise(avg=mean(mean))

# venn drawing fun
save_venn <- function(vct){
    p.venn <- draw.pairwise.venn(vct[1], vct[2], vct[3], 
        fill=brewer.pal(n=3,name="Set1")[c(1,2)],
        alpha=c(0.6,0.6),lty=c("blank","blank"),
        scaled=TRUE,
        euler.d=TRUE,
        margin=c(0.05,0.05),
        ind=FALSE,
        sep.dist=0.001,
        ext.text=FALSE,
        fontfamily=rep("sans",3),
        cex=rep(1.5,3))
    return(p.venn)
}

venn.all <- lapply(inters.three,save_venn)
venn.all.grob <- lapply(venn.all,grobTree)
pg <- plot_grid(plotlist=venn.all.grob, nrow=4, ncol=5)#, labels=names(venn.all.grob)
save_plot("../temp/figs/venn.pdf", pg, nrow=4, ncol=5)
