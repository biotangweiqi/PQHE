# PQHE
**Pooled QTL Heritability Estimator (PQHE)** is a method to estimate QTL heritability using pooled sequencing data obtained under different experimental designs. This method was implemented via R programming. The R script of PQHE can be used at command line from terminal.

# Introduction
Bulked segregant analysis (BSA) is originally proposed for rapid identification of molecular markers associated with a trait of interest. In recent years, next generation sequencing (NGS) technologies have been applied to BSA. By referring to a known genome sequence, a huge number of markers (mainly single nucleotide polymorphisms, SNPs) can be detected and mapped simultaneously by NGS. This greatly enhances the power of BSA. Such NGS-assisted BSA (abbr. BSA-seq) has proven to be an efficient and cost-effective method for quick mapping of single major genes as well as quantitative trait loci (QTLs). QTL mapping by BSA-seq, or termed ‘pooled QTL mapping’, is very fascinating because it requires to genotype (sequence) only a pair of DNA pools in-stead of many individuals. 

Although pooled QTL mapping has the significant advantage of being simple and fast, it also has an obvious weakness. In conventional QTL mapping studies, both QTL location and QTL heritability (the proportion of phenotypic variation contributed by a QTL) can be estimated. Up to now, however, QTL heritability estimation in pooled QTL mapping has not been realized.

To solve this weakness, we propose a method for pooled QTL heritability estimation. To make this method widely useful, this method was implemented by R language and the R script is provided in this project (PQHE).

# Download and installation
R script (pooled_QTL_heritability_estimator.R) can be directly downloaded and used in batch mode of R.

# Required R package
In the R script (pooled_QTL_heritability_estimator.R), R package [rootSolve](https://cran.r-project.org/web/packages/rootSolve) is used for solving nonlinear equations, thus rootSolve need be installed in R.

    install.packages("rootSolve")

# Usage
At command line interface, enter command like this:

    Rscript pooled_QTL_heritability_estimator.R example_conf.txt
