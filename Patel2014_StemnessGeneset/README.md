# Produce Stemness Geneset from Patel 2014

Here, I re-analyze a part of
[Patel et al](http://science.sciencemag.org/content/344/6190/1396.long).
This is a great publication about glioblastoma intratumoral heterogeneity.
One part of the work includes a scRNA-seq experiment.
In a comparative transcriptomics experiment a geneset indicating *stemness* is derived.
It is interesting to see how the results change if current state-of-the-art scRNA-seq analysis methods are used.
All details are given [here](http://b210-research.dkfz.de/computational-genome-biology/scRNAseq/stemnessGeneset/).

## Quality Control and Normalization

Quality control and normalization are done using the `scater` framework.
A standard quality control workflow is used and `scran` library size normalization is done.
Furthermore, the batch effect is removed using a linear mixed model (this can actually take a while to compute).

    Rscript -e 'rmarkdown::render("QC.Rmd")'
    Rscript -e 'rmarkdown::render("Norm.Rmd")'

This will create html reports about quality control and normalization,
as well as some intermediate data.

## Stemness Geneset

First, a stemness geneset is derived using scDD and edgeR (scDD runs for a while).
Then this geneset is compared to the one derived by Patel *et al*.
GSEAs for both genesets are performed and compared as well.

    Rscript -e 'rmarkdown::render("benchmarking_results.Rmd")'
    Rscript -e 'rmarkdown::render("consensus_approach.Rmd")'

This will create 2 html reports, some intermediate data, and some results.

**Packages**

Before you start, make sure you have all the `R` packages installed:
`scran`, `scater`, `data.table`, `ggplot2`, `limma`, `MultiAssayExperiment`, `SummarizedExperiment`, 
`edgeR`, `scDD`, `DT`, `Lattirl`, `ggthemes`.
