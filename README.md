# primer-faceoff

Final data and code for the primer-faceoff paper published in _???_:

Collins, R. A., Bakker, J., Wangensteen, O. S., Soto, A. Z., Corrigan, L., Sims, D. W., Genner, M. J., Mariani, S. (XXX) Non-specific amplification compromises environmental DNA metabarcoding with COI. _XXX_. [https://doi.org/XXX](https://doi.org/XXX).

### Contents

* **`code/`** - R and shell scripts to run analyses.
    - `blast.sh` - run local blast searches
    - `dada2.R` - run the dada2 workflow to denoise the sequence data
    - `demultiplex.sh` - demultiplex with cutadapt
    - `epa.sh` - run the Evolutionary Placement Algorithm
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