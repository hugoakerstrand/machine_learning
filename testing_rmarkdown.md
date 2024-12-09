---
title: "Testing R Markdown :O"
output: 
  html_document: 
    keep_md: true
date: "2024-12-07"
---


``` r
library(flowCore)
library(mclust)
```

```
## Warning: package 'mclust' was built under R version 4.4.2
```

```
## Package 'mclust' version 6.1.1
## Type 'citation("mclust")' for citing this R package in publications.
```

``` r
data("GvHD")
clustered <- Mclust(GvHD.pos)
```
