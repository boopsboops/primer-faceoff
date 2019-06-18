#!/usr/bin/env sh
# MFE primer checks for primer specificity against a user supplied database

## installing MFEprimer
# install dependancy
sudo apt install python-psutil
sudo pip install psutil --upgrade
# cd and download git repo
cd Software
git clone https://github.com/quwubin/MFEprimer.git
cd MFEprimer
# add to path
export PATH=~/Software/MFEprimer:$PATH
#echo 'export PATH=~/Software/MFEprimer:$PATH' >> ~/.bashrc

# first make mitogenome fasta and primer fastas with 'primer-rates.R' #
# first make mitogenome fasta and primer fastas with 'primer-rates.R' #
# first make mitogenome fasta and primer fastas with 'primer-rates.R' #


# make the MFEprimer db (k=9)
IndexDb.sh ../temp/mitogenome/mitogenome.fasta 9

# get list of fas primers
fas="$(ls ../temp/mitogenome/*.fas)"
for f in $fas; do echo "$f"; done

for i in $fas; do
	MFEprimer.py -k 9 --ppc 0 --tab --size_start 50 --size_stop 800 -i "$i" -d ../temp/mitogenome/mitogenome.fasta > "$i".results.tsv &
done
