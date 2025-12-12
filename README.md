# PCA-MR: PCA-Based Mendelian Randomization

A package for performing two-sample mendeleian randomization with summary statistics.

### Install from source

``` r
install.packages("devtools")
devtools::install_github("chrislaw133/PCA-MR")
```

***Important:*** Make sure your alleles are aligned, and that your vectors are numeric and in the same order as the LD matrix!

### Read in plink .zst square matrix

``` r
library(data.table)
ld <- as.matrix(fread(cmd = "unzstd -c plink.zst"))
```

### Read in vectors 

``` r
vector <- as.numeric(readLines("path/to/vector.txt")
```

### Usage

``` r
pca_mr(bx, by, sey, ld)
```





