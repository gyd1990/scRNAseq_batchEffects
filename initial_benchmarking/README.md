# Initial Benchmarking

## Simulation

The simulation is in *R/generate_datasets.R*.
The functions `PoiBeta` and `Generator` do all the work.

## Reproduce the 4 Scenarios

There are 4 python scripts for reproducing each of the 4 scenarios.
Execute them in this directory. 
E.g.

    ./strongDD_strongBatch.py

to create the scenario with a strong differential expression and a strong batch effect (this can take a while).
A directory will be created with 3 subdirectories for each of the 3 simulated datasets.
The parameters for generating the datasets are written into a csv file.
In the subdirectories of each simulation you will find the 3 gitbooks about *Data Quality*, *Benchmarking*, and *Data Characterization*.
The actual `sce` dataset can be found in *intermediate_data*, the results of the benchmarking in *results*.


