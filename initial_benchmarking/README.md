# Initial Benchmarking

Please find the code on [github](https://github.com/mRcSchwering/scRNAseq_batchEffects/tree/master/batchEffect_explo).
Clone the repository.
The code for this experiment is in *initial_benchmarking*.

    git clone https://github.com/mRcSchwering/scRNAseq_batchEffects
    cd scRNAseq_batchEffects/initial_benchmarking/


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

**Packages**

Before you start, make sure you have all the `R` packages installed (it would be annoying to find out the script halted because a package is not installed):
`scran`, `scater`, `data.table`, `ggplot2`, `limma`, `fitdistrplus`, `vcd`, `FAdist`, `edgeR`, `BPSC`, `scDD`, `SummarizedExperiment`, `MAST`, `bookdown`.
