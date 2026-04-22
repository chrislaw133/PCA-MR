# PCA-MR: PCA-Based Mendelian Randomization

A package for performing two-sample mendelian randomization with summary statistics.

### Install from source

``` r
install.packages("devtools")
devtools::install_github("chrislaw133/PCA-MR")
```

***Important:*** Make sure your alleles are aligned, and that your vectors are numeric and in the same order as the LD matrix!

### Read in plink .zst square matrix (SNP correlation matrix)

``` r
library(data.table)
ld <- as.matrix(fread(cmd = "unzstd -c plink.zst"))
```

### Read in vectors 

``` r
vector <- as.numeric(readLines("path/to/vector.txt")
```

## Usage

### PCA-MR-IVW

``` r
library(PCAMR)

fit_ivw <- pcamr_ivw(
  bx = bx,
  by = by,
  bxse = bxse,
  byse = byse,
  ld = ld,
  model = "random" #"fixed" or "random", "random" is recommended
)
```

### PCA-MR-ML

``` r
library(PCAMR)

fit_ml <- pcamr_ml(
  bx = bx,
  by = by,
  bxse = bxse,
  byse = byse,
  ld = ld,
  alpha = 0.05,
  model = "random" #"fixed" or "random", "random" is recommended
)
```

### PCA-MR-Egger

``` r
library(PCAMR)

fit_egger <- pcamr_egger(
  bx = bx,
  by = by,
  bxse = bxse,
  byse = byse,
  ld = ld
)
```

##References

##Yavorska, O. O., & Burgess, S. (2017). MendelianRandomization: An R package for performing Mendelian randomization analyses using summarized data. International Journal of Epidemiology, 46(6), 1734–1739. https://doi.org/10.1093/ije/dyx034

##Sun, J., Dong, Q., Wei, J., Gao, Y., Yu, Z., Hu, X., & Zhang, Y. (2025). ti-scMR: Trajectory-inference-based dynamic single-cell Mendelian randomization identifies causal genes underlying phenotypic differences. NAR Genomics and Bioinformatics, 7(3), lqaf082. https://doi.org/10.1093/nargab/lqaf082





