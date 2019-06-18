#!/usr/bin/env/ sh
 
# blast ALL of genbank locally for ALL sequences #
# blast ALL of genbank locally for ALL sequences #
# blast ALL of genbank locally for ALL sequences #
# blast ALL of genbank locally for ALL sequences #
# blast ALL of genbank locally for ALL sequences #


### First download and assemble the BLAST nt database (only need to do this step once)
# cd to directory
cd ../temp

# README - blast search carried out July 30th 2018
# download the v5 archives from ncbi
wget -v "ftp://ftp.ncbi.nlm.nih.gov/blast/db/v5/nt_v5.??.tar.gz" -P blastdbv5

# untar them into a single folder
mkdir -p blastdbv5/nt
cd blastdbv5
for i in *.tar.gz; do tar xvfz $i -C nt/; done

# move the taxonomy db
cp nt/taxdb.bti ../fastq/taxdb.bti
cp nt/taxdb.btd ../fastq/taxdb.btd
cd ../fastq
###


# do sequence homology search
# set marker
MARKER="miya"
MARKER="leray-xt"
MARKER="sea-short"
MARKER="sea-mid"
#
HMM="12s.miya"
HMM="coi.lerayxt"
HMM="coi.seashort"
HMM="coi.seamid"

# subset just the fragment using an HMM
hmmsearch -E 0.01 --incE 0.01 ../../../reference-libraries/hmms/"$HMM".noprimers.hmm "$MARKER"/results/asvs.fna | grep ">> asv[0-9+]" | sed -e 's/>> //g' -e 's/[[:space:]]//g' -e 's/$/$/' | sort | uniq > "$MARKER"/results/hmm-out.txt
grep -A 1 -f "$MARKER"/results/hmm-out.txt "$MARKER"/results/asvs.fna | grep -v "-" > "$MARKER"/results/asvs-clean.fasta
rm "$MARKER"/results/hmm-out.txt

# now break up the 'otus-clean.fasta' in 'taxonomic-assignment.R'
# now break up the 'otus-clean.fasta' in 'taxonomic-assignment.R'
# now break up the 'otus-clean.fasta' in 'taxonomic-assignment.R'

# date/time of file modification
#stat -c %y blast1.fa

# run a bulk blast search
date +"%H:%M"
for i in "$MARKER"/results/*.fa; do blastn -db ../blastdbv5/nt/nt_v5 -query "$i" -num_threads 1 -task blastn -evalue 1000 -word_size 11 -max_target_seqs 500 -outfmt "6 qseqid sseqid sacc score bitscore evalue length pident nident staxids sscinames scomnames sskingdoms" -out "$i".out & done

# make a header file
echo -e "qseqid\tsseqid\tsacc\tscore\tbitscore\tevalue\tlength\tpident\tnident\tstaxids\tsscinames\tscomnames\tsskingdoms" > "$MARKER"/results/headers
cat "$MARKER"/results/headers "$MARKER"/results/*.fa.out > "$MARKER"/results/blast-results.tsv
rm "$MARKER"/results/headers



# blast reference database locally for only fish reads #
# blast reference database locally for only fish reads #
# blast reference database locally for only fish reads #
# (first pull out fish sequences using 'taxonomic-assignment.R')

# set marker
MARKER="miya"
MARKER="leray-xt"
MARKER="sea-short"
MARKER="sea-mid"

# make blast db (only need to do this step once)
makeblastdb -in references.fasta -parse_seqids -dbtype nucl -blastdb_version 5

# blast the blast db
# get better hits with smaller word size
blastn -task blastn -num_threads 4 -evalue 1000 -word_size 7 -max_target_seqs 500 -db references.fasta -outfmt "6 qseqid sseqid evalue length pident nident score bitscore" -out "$MARKER"/results/fish-blast.out -query "$MARKER"/results/fishqueries.fasta

# join the header
echo -e "qseqid\tsseqidLocal\tevalueLocal\tlengthLocal\tpidentLocal\tnidentLocal\tscoreLocal\tbitscoreLocal" > "$MARKER"/results/headers
cat "$MARKER"/results/headers "$MARKER"/results/fish-blast.out > "$MARKER"/results/fish-blast-result.tsv
rm "$MARKER"/results/fish-blast.out
rm "$MARKER"/results/headers

# now go back to 'taxonomic-assignment.R' to process
