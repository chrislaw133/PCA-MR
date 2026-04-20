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





