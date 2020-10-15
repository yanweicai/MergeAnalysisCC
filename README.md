Merge analysis tests whether the QTL signal identified by haplotype association can be explained by the pattern of alleles at a known variant by comparing the fit of the founder allele haplotype model with a simpler “merged” model in which haplotypes carrying the same variant allele are forced to have the same haplotype effect, or equivalently, haplotypes with the same allele are “merged”. 

This repository is specifically designed to run merge analysis in the Colloabrative Cross mice, using the ISBdb database (http://github.com/danoreper/ISVdb)

The example here use a toy example and run by
```R
Rscript MergeExample.r
```

Input variables | Discription
------------ | -------------
genodb       |   ISVdb database path, prefer to be on clusters
CHR          | Chomesome location of intrested peaks
peakI        |Peak interval, a two element vector showing the start and end of the QTL in MB. eg.c(39,41)
tmp.file.path | Path of output folder
