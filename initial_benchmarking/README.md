# Initial Benchmarking

There are already many benchmarkings concerning differential expression in scRNA-seq data. However, they are missing at least of one the following points.

1. Many benchmarkings do not test most recent methods (like scDD, D3E, BPSC).
2. They test differential means but not other kinds of differential distributions (as elaborated in Korthauer (2016) for example).
3. They test one group vs another without introducing a batch effect. However, batch effects in single cell data are substantial according to Hicks (2015) and Tung (2016).

All details are given [here](http://b210-research.dkfz.de/computational-genome-biology/scRNAseq/initial_benchmarking/).

## Simulation

The simulation is in *R/generate_datasets.R*.
The functions `PoiBeta` and `Generator` do all the work.

## Reproduce the 4 Scenarios

There are 4 python scripts for reproducing each of the 4 scenarios.
Execute them in the subdirectory *initial_benchmarking*:

    ./strongDD_strongBatch.py
    ./weakDD_strongBatch.py
    ./strongDD_weakBatch.py
    ./weakDD_weakBatch.py

For each scenario a subdirectory with 3 simulations is created.
Then for each simulation the reports are created (this will run for several hours).

The parameters for generating the datasets are written into a csv file.
In the subdirectories of each of the 3 simulations you will find the 3 gitbooks about *Data Quality*, *Benchmarking*, and *Data Characterization*.
The actual `sce` object of the dataset can be found in *intermediate_data*, the results of the benchmarking in *results*.

## Analyse Results

In two reports results are analysed (*benchmarking_results*) and a scDD-edgeR ensemble is derived (*consensus_approach*).
Render both reports with:

    Rscript -e 'rmarkdown::render("benchmarking_results.Rmd")'
    Rscript -e 'rmarkdown::render("consensus_approach.Rmd")'

They will create html reports and a few *rds* files.

**Packages**

Before you start, make sure you have all the `R` packages installed (it would be annoying to find out the script halted because a package is not installed):
`scran`, `scater`, `data.table`, `ggplot2`, `limma`, `fitdistrplus`, `vcd`, `FAdist`, `edgeR`, `BPSC`, `scDD`, `SummarizedExperiment`, `MAST`, `bookdown`.
