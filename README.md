# ABC Regulatory Density Permutation Analysis

R script for permutation testing of regulatory density differences between female and male sex-stratified obesity variants using ABC enhancer-gene predictions.

## Method
Compares the ratio of target genes to credible SNPs (regulatory density) across 22 tissue categories. Statistical significance assessed via 10,000 permutations with FDR correction (Benjamini-Hochberg).

## Requirements
- Base R (no packages required)

## Input
CSV files with ABC enhancer-gene mappings containing:
- `gene_name`: Target gene
- `variant_id`: SNP identifier  
- `Tissue`: Tissue category
- `sex`: "Female" or "Male"

## Usage
```r
# Edit file paths in lines 9-10, then run:
source("abc_regulatory_density_permutation.R")
