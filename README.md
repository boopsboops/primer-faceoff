# Non-specific amplification compromises environmental DNA metabarcoding with COI

Final data and code for the primer comparison paper published in _Methods in Ecology and Evolution_:

Collins, R.A., Bakker, J., Wangensteen, O.S., Soto, A.Z., Corrigan, L., Sims, D.W., Genner, M.J., Mariani, S. (2019) Non-specific amplification compromises environmental DNA metabarcoding with COI. _Methods in Ecology and Evolution_. [https://doi.org/10.1111/2041‐210X.13276](https://doi.org/10.1111/2041‐210X.13276).


### Contents

* **`code/`** - R and shell scripts to run analyses.
    - `blast.sh` - run local blast searches
    - `dada2.R` - run the dada2 workflow to denoise the sequence data
    - `demultiplex.sh` - demultiplex with cutadapt
    - `epa.sh` - run the Evolutionary Placement Algorithm
    - `get-data.sh` - download fastq files from Dryad
    - `id-rates.R` - calculate marker discriminatory power
    - `make-plots.R` - generate the plots from the manuscript
    - `MFEprimer.sh` - run _in silico_ PCR with MFE primer
    - `prep-barcodes.R` - generate oligo tag files for demultiplexing
    - `primer-rates.R` - process _in silico_ PCR results
    - `print-tables.R` - generate the tables from the manuscript
    - `swarm.sh` - run OTU clustering with Swarm
    - `taxonomic-assignment.R` - process blast and EPA assignments

* **`data/`** - Data generated from the analyses.
    - `blast-results-all.csv` - top blast hits for each ASV
    - `fish-assignments.csv` - species assigned to each fish ASV by blast/EPA
    - `in-silico-results.csv` - _in silico_ results (Table 2)
    - `primers.csv` - PCR primers and citations (Table 1)
    - `results-by-marker.csv` - number of reads per event/species/marker/sample
    - `sample-plates.csv` - PCR primers and oligo tags for demultiplexing
    - `sequencing-summary.csv` - bioinformatic processing results (Table 3)
    - `traditional-survey.csv` - traditional survey results

* **`temp/`** - Temporary file directory that is not committed to the repository, but needs to be created locally first, to run the scripts and see the output. Ignored by git.