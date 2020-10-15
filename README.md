# Merge Analysis Introduction
Merge analysis tests whether the QTL signal identified by haplotype association can be explained by the pattern of alleles at a known variant by comparing the fit of the founder allele haplotype model with a simpler “merged” model in which haplotypes carrying the same variant allele are forced to have the same haplotype effect, or equivalently, haplotypes with the same allele are “merged”. 

# Collabrative Cross and ISVdb

For example, consider a variant with k=3 alleles and strain distribution pattern 00111222, meaning that allele 0 is shared by the first 2 founders (129S1, AJ), allele 1 shared by the next 3 (B6, NOD, NZO) and allele 2 is in the last 3 founders (CAST, PWK, WSB). In the merged model for this variant, the QTL term in Eq 3 would have only 3 effects, $\text{QTL}_{sm}=\beta_0 \textit{x}_0+\beta_1 \textit{x}_1+\beta_2 \textit{x}_2$, where, for example, $\textit{x}_1=\textit{x}_{\text{B6},sm} + \textit{x}_{\text{NOD},sm} + \textit{x}_{\text{NZO},sm}$

# Code 
MergeExample.r
