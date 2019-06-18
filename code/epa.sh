#!/usr/bin/env/ sh

# install 'epa-ng' and 'gappa' from github

# set marker
MARKER="miya"
MARKER="leray-xt"
MARKER="sea-short"
MARKER="sea-mid"

# make RAxML trees for all markers
# run these in different terminals (leray-xt takes about 1.5 h)
# cd to wd (RAxML can't do relative file paths)
cd ../temp/fastq/"$MARKER"/results/epa/
raxmlHPC-AVX -s references.aligned.fasta -f d -m GTRCATX -n raxml -p 42
mv RAxML_bestTree.raxml references.tree.nwk
rm *.raxml

# get model params from raxml
raxmlHPC-AVX -f e -s references.aligned.fasta -t references.tree.nwk -n txt -m GTRGAMMAX
mv RAxML_info.txt raxml.model.params
mv RAxML_result.txt references.tree.nwk
rm *.txt
rm *.reduced

# run EPA
epa-ng --ref-msa references.aligned.fasta --tree references.tree.nwk --query queries.aligned.fasta --outdir . --model raxml.model.params --redo --preserve-rooting off

# run gappa
gappa analyze assign --jplace-path epa_result.jplace --taxon-file references.taxonomy.tsv
