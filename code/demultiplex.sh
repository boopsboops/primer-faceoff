#!/usr/bin/env sh

# cd to directory
cd ../temp/fastq

# set params #

# MIYA
MARKER="miya"
FWD="GTCGGTAAAACTCGTGCCAGC"
REV="CATAGTGGGGTATCTAATCCCAGTTTG"
RAWR1="SeaDNA-Fish5_S2_L001_R1_001.fastq.gz"
RAWR2="SeaDNA-Fish5_S2_L001_R2_001.fastq.gz"
MINLEN="21"
TRUCVAL="105"

# LERAY-XT
MARKER="leray-xt"
FWD="GGWACWRGWTGRACWNTNTAYCCYCC"
REV="TANACYTCNGGRTGNCCRAARAAYCA"
RAWR1="SFI1_real_R1_001.fastq.gz"
RAWR2="SFI1_real_R2_001.fastq.gz"
MINLEN="26"
TRUCVAL="170"

# SEA-SHORT
MARKER="sea-short"
FWD="GGAGGCTTTGGMAAYTGRYT"
REV="GGGGGAAGAARYCARAARCT"
RAWR1="SeaDNA-Fish3_S3_L001_R1_001.fastq.gz"
RAWR2="SeaDNA-Fish3_S3_L001_R2_001.fastq.gz"
MINLEN="20"
TRUCVAL="40"

# SEA-MID
MARKER="sea-mid"
FWD="GGAGGCTTTGGMAAYTGRYT"
REV="TAGAGGRGGGTARACWGTYCA"
RAWR1="SeaDNA-Fish4_S1_L001_R1_001.fastq.gz"
RAWR2="SeaDNA-Fish4_S1_L001_R2_001.fastq.gz"
MINLEN="20"
TRUCVAL="80"

# make a grepable degenerate primer
FWDGR="$(echo "$FWD" | sed -e 's/W/[A\|T]/g' -e 's/R/[A\|G]/g' -e 's/Y/[C\|T]/g' -e 's/M/[A\|C]/g' -e 's/N/[A\|C\|G\|T]/g' -e 's/$/\|\$/g')"
REVGR="$(echo "$REV" | sed -e 's/W/[A\|T]/g' -e 's/R/[A\|G]/g' -e 's/Y/[C\|T]/g' -e 's/M/[A\|C]/g' -e 's/N/[A\|C\|G\|T]/g' -e 's/$/\|\$/g')"


# make dirs for logs and sense/antisense, results etc
mkdir -p "$MARKER"/logs "$MARKER"/sense/dmplx "$MARKER"/sense/trimmed "$MARKER"/antisense/dmplx "$MARKER"/antisense/trimmed "$MARKER"/results "$MARKER"/trash

# make a stats file
# most of the stats can be found in cutadapt logs
echo -e "stat,reads\npf,\naquatheat,\nprimer,\nbarcode,\ntrim,\nfilter,\nmerge,\nchim,\nhomol," > "$MARKER"/logs/stats.csv


# demultiplex by orientation (need to wait for fwd to finish before doing rev)
# "--pair-adapters" means that the primers must be present in both forward and reverse to be RETAINED

# sense
cutadapt --error-rate 0.15 --overlap "$MINLEN" --pair-adapters --action=none -g senseF="$FWD" -G senseR="$REV" --untrimmed-output "$MARKER"/trash/untrimmed.R1.fastq.gz --untrimmed-paired-output "$MARKER"/trash/untrimmed.R2.fastq.gz -o "$MARKER"/sense/R1.fastq.gz -p "$MARKER"/sense/R2.fastq.gz "$RAWR1" "$RAWR2" > "$MARKER"/logs/cutadapt.sense.log

# antisense
cutadapt --error-rate 0.15 --overlap "$MINLEN" --pair-adapters --action=none -g antisenseF="$REV" -G antisenseR="$FWD" --untrimmed-output "$MARKER"/trash/noprimer.R1.fastq.gz --untrimmed-paired-output "$MARKER"/trash/noprimer.R2.fastq.gz -o "$MARKER"/antisense/R1.fastq.gz -p "$MARKER"/antisense/R2.fastq.gz "$MARKER"/trash/untrimmed.R1.fastq.gz "$MARKER"/trash/untrimmed.R2.fastq.gz > "$MARKER"/logs/cutadapt.antisense.log

# check output and trash for primers 
gzip -cd "$MARKER"/sense/R1.fastq.gz | sed -n '2~4p' | head -n 40 | GREP_COLOR="1;44" grep --color=always -E "$FWDGR" | GREP_COLOR="1;41" grep --color=always -E "$REVGR"
gzip -cd "$MARKER"/trash/noprimer.R1.fastq.gz | sed -n '2~4p' | head -n 40 | GREP_COLOR="1;44" grep --color=always -E "$FWDGR" | GREP_COLOR="1;41" grep --color=always -E "$REVGR"
gzip -cd "$MARKER"/trash/noprimer.R2.fastq.gz | sed -n '2~4p' | head -n 40 | GREP_COLOR="1;44" grep --color=always -E "$FWDGR" | GREP_COLOR="1;41" grep --color=always -E "$REVGR"

# Demultiplex by barcode SENSE
# First generate barcodes files with 'prep-barcodes.R'
cutadapt --no-indels --error-rate 0.1 --overlap 10 --action=none -g file:"$MARKER"/barcodes-sense.fas -o "$MARKER"/sense/dmplx/{name}.R1.fastq.gz -p "$MARKER"/sense/dmplx/{name}.R2.fastq.gz "$MARKER"/sense/R1.fastq.gz "$MARKER"/sense/R2.fastq.gz > "$MARKER"/logs/cutadapt.dmplx.barcodes.sense.log
mv "$MARKER"/sense/dmplx/unknown.R1.fastq.gz "$MARKER"/trash/sense.unknown.R1.fastq.gz
mv "$MARKER"/sense/dmplx/unknown.R2.fastq.gz "$MARKER"/trash/sense.unknown.R2.fastq.gz

# Demultiplex by barcode ANTISENSE
cutadapt --no-indels --error-rate 0.1 --overlap 10 --action=none -g file:"$MARKER"/barcodes-antisense.fas -o "$MARKER"/antisense/dmplx/{name}.R1.fastq.gz -p "$MARKER"/antisense/dmplx/{name}.R2.fastq.gz "$MARKER"/antisense/R1.fastq.gz "$MARKER"/antisense/R2.fastq.gz > "$MARKER"/logs/cutadapt.dmplx.barcodes.antisense.log
mv "$MARKER"/antisense/dmplx/unknown.R1.fastq.gz "$MARKER"/trash/antisense.unknown.R1.fastq.gz
mv "$MARKER"/antisense/dmplx/unknown.R2.fastq.gz "$MARKER"/trash/antisense.unknown.R2.fastq.gz

# check the kept and discarded sequences
cat "$MARKER"/sense/dmplx/*.R1.fastq.gz | gzip -cd | sed -n '2~4p' | GREP_COLOR="1;44" grep --color=always -E "$FWDGR" | GREP_COLOR="1;41" grep --color=always -E "$REVGR"
gzip -cd "$MARKER"/trash/sense.unknown.R1.fastq.gz | sed -n '2~4p' | head -n 40 | GREP_COLOR="1;44" grep --color=always -E "$FWDGR" | GREP_COLOR="1;41" grep --color=always -E "$REVGR"
gzip -cd "$MARKER"/trash/sense.unknown.R2.fastq.gz | sed -n '2~4p' | head -n 40 | GREP_COLOR="1;44" grep --color=always -E "$FWDGR" | GREP_COLOR="1;41" grep --color=always -E "$REVGR"


# trim SENSE with cutadapt
sense="$(ls "$MARKER"/sense/dmplx/*.fastq.gz | sed --expression='s/\.R.\.fastq\.gz//g' | sed --expression='s/unknown//g' | sed --expression="s/$MARKER\/sense\/dmplx\///g" | uniq)"
# check
for f in $sense; do echo "$f"; done
# now run
for i in $sense; do
cutadapt -n 5 --error-rate 0.15 --minimum-length "$TRUCVAL" -g "$FWD" -G "$REV" -o "$MARKER"/sense/trimmed/"$i".R1.fastq.gz -p "$MARKER"/sense/trimmed/"$i".R2.fastq.gz --discard-untrimmed "$MARKER"/sense/dmplx/"$i".R1.fastq.gz "$MARKER"/sense/dmplx/"$i".R2.fastq.gz >> "$MARKER"/logs/cutadapt.sense.trimming.log &
done

# trim ANTISENSE with cutadapt
antisense="$(ls "$MARKER"/antisense/dmplx/*.fastq.gz | sed --expression='s/\.R.\.fastq\.gz//g' | sed --expression='s/unknown//g' | sed --expression="s/$MARKER\/antisense\/dmplx\///g" | uniq)"
# check
for f in $antisense; do echo "$f"; done
# now run
for i in $antisense; do
cutadapt -n 5 --error-rate 0.15 --minimum-length "$TRUCVAL" -g "$REV" -G "$FWD" -o "$MARKER"/antisense/trimmed/"$i".R1.fastq.gz -p "$MARKER"/antisense/trimmed/"$i".R2.fastq.gz --discard-untrimmed "$MARKER"/antisense/dmplx/"$i".R1.fastq.gz "$MARKER"/antisense/dmplx/"$i".R2.fastq.gz >> "$MARKER"/logs/cutadapt.antisense.trimming.log &
done

# check output for primers (should be none)
cat "$MARKER"/sense/trimmed/*.R1.fastq.gz | gzip -cd | sed -n '2~4p' | GREP_COLOR="1;44" grep --color=always -E "$FWDGR" | GREP_COLOR="1;41" grep --color=always -E "$REVGR"
cat "$MARKER"/antisense/trimmed/*.R1.fastq.gz | gzip -cd | sed -n '2~4p' | GREP_COLOR="1;44" grep --color=always -E "$FWDGR" | GREP_COLOR="1;41" grep --color=always -E "$REVGR"


# join and get lengths and counts 
cat "$MARKER"/sense/trimmed/*.R1.fastq.gz | fqtools lengthtab
cat "$MARKER"/sense/trimmed/*.R2.fastq.gz | fqtools lengthtab

# get counts
cat "$MARKER"/sense/dmplx/*.R1.fastq.gz | fqtools count
cat "$MARKER"/antisense/dmplx/*.R1.fastq.gz | fqtools count

cat "$MARKER"/sense/trimmed/*.R1.fastq.gz | fqtools count
cat "$MARKER"/antisense/trimmed/*.R1.fastq.gz | fqtools count

# count number of aquatheatr reads
cat "$MARKER"/sense/dmplx/*AquaTheat*R1.fastq.gz | fqtools count
cat "$MARKER"/antisense/dmplx/*AquaTheat*R1.fastq.gz | fqtools count

####################################
# NOW MOVE TO DADA2 IN R 
####################################
####################################

# counts after filter in dada2
cat "$MARKER"/sense/filtered/*.R1.fastq.gz | fqtools count
cat "$MARKER"/antisense/filtered/*.R1.fastq.gz | fqtools count
