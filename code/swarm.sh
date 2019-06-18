#!/usr/bin/env/ sh

# install 'swarm' from https://github.com/torognes/swarm

# set marker
MARKER="miya"
MARKER="leray-xt"
MARKER="sea-short"
MARKER="sea-mid"

# run swarm
swarm -t 4 -d 1 -f -z -o ../temp/fastq/"$MARKER"/results/swarm/swarm.out ../temp/fastq/"$MARKER"/results/swarm/fishqueries-clean-abundances.fasta

# clean swarm output
sed -e 's/;size=[0-9]*//g' ../temp/fastq/"$MARKER"/results/swarm/swarm.out | nl -w 1 | sed -e 's/^/swarm/g' -e 's/ /\t/g' > ../temp/fastq/"$MARKER"/results/swarm/swarm.tsv
