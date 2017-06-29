# Batch Effect Exploration

Please find the code on [github](https://github.com/mRcSchwering/scRNAseq_batchEffects/tree/master/batchEffect_explo).
Clone the repository.
The code for this experiment is in *batchEffect_explo*.

    git clone https://github.com/mRcSchwering/scRNAseq_batchEffects
    cd scRNAseq_batchEffects/batchEffect_explo/


## Reproduce Simulations

For both grids (in subdirectories *grid1/*, *grid2/*) the python script
*run_gridSearch.py* will create all simulations (768 for each grid).

    cd grid1/
    ./run_gridSearch.py
    cd ../grid2/
    ./run_gridSearch.py
    cd ..

This will run for several hours and take up around 80GB.
For each parameter combination a directory is created
(*n500_b2.0_g1.0_c0.5/* for *n = 500*, *b = 2.0*, *g = 1.0*, *c = 0.5*)
with 3 subdirectories for each simulation.

They contain html reports for quality control (*quali.html*), differential expression calling
(*call.html*), result comparison (*compare.html*), and batch effect estimation
(*estimate.html*).
The actual `sce` object of the dataset can be found in *intermediate_data*, result objects are in *results*.

## Regress FDR Control

After the above step, results are gathered and pre-processed,
and different models are fit to predict the loss of FDR control from the
parameters estimated form the datasets.
This is done for scDD using pooled cells and the Kolmogorov-Smirnov test,
and for edgeR using summed up batches.

    Rscript -e 'rmarkdown::render("pre.Rmd")'
    Rscript -e 'rmarkdown::render("model_scDD.Rmd")'
    Rscript -e 'rmarkdown::render("model_edgeR.Rmd")'

Some of the models require hyper parameter tuning,
which is why the respective scripts might run over 1 hour.
Each of them will create a html report and several  *rds* files.

**Packages**

Before you start, make sure you have all the `R` packages installed:
`scran`, `scater`, `data.table`, `ggplot2`, `edgeR`, `DESeq2`, `scDD`, `SummarizedExperiment`, `rmarkdown`,
`kBET`, `limma`, `earth`, `mgcv`, `xgboost`, `scatterplot3d`, `corrplot`, `highcharter`.
