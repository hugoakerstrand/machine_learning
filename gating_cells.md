---
title: "Machine learning and facs gating: Cells (FSC-H vs SSC-H)"
author: "Hugo Akerstrand"
date: "2024-12-09"
output: 
  html_document: 
    keep_md: true
---

# Welcome!

This is a document that aims to explore different ways of gating flow cytometry
data in R, primarily using a tidymodels approach.

To do this, I am using the data set `GvHD` from `flowCore`. Note that I am not 
compensating or cleaning the data, since this is not meant to represent a complete
flow cytometry analysis pipeline.



Start by loading the data and make it into a tibble with the relevant list columns:


``` r
gvhd_tibble <- tibble(
  exprs = purrr::map(GvHD, ~ exprs(.x)),                       # This contains detector information
  keywords = purrr::map(GvHD, ~ keyword(.x)),                  # This contains meta data
  exprs_tibble = purrr::map(exprs, function(.x) as_tibble(.x)) # This is for plotting
)

head(gvhd_tibble)
```

```
## # A tibble: 6 × 3
##   exprs              keywords           exprs_tibble         
##   <list>             <list>             <list>               
## 1 <dbl [3,420 × 8]>  <named list [170]> <tibble [3,420 × 8]> 
## 2 <dbl [3,405 × 8]>  <named list [170]> <tibble [3,405 × 8]> 
## 3 <dbl [3,435 × 8]>  <named list [170]> <tibble [3,435 × 8]> 
## 4 <dbl [8,550 × 8]>  <named list [170]> <tibble [8,550 × 8]> 
## 5 <dbl [10,410 × 8]> <named list [170]> <tibble [10,410 × 8]>
## 6 <dbl [3,750 × 8]>  <named list [170]> <tibble [3,750 × 8]>
```

We will gate the 'cell gate' using `FSC-H` and `SSC-H` columns. First inspect the data:


``` r
plot_list <- gvhd_tibble$exprs_tibble |> purrr::map(function(.x) {
      ggplot(data= .x, aes(x = `FSC-H`, y = `SSC-H`)) +
        geom_hex()
    }
  )

wrap_plots(plot_list[1:4])
```

![](gating_cells_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

We clean up the data by filtering out FSC-H < 100:


``` r
gvhd_tibble <- gvhd_tibble |> 
  mutate(exprs_tibble = purrr::map(gvhd_tibble$exprs_tibble, ~ filter(.x, `FSC-H` > 99)))

plot_list <- gvhd_tibble$exprs_tibble |> purrr::map(function(.x) {
      ggplot(data= .x, aes(x = `FSC-H`, y = `SSC-H`)) +
        geom_hex() +
    scale_x_continuous(limits = c(0,1000))
    }
  )

wrap_plots(plot_list[1:4])
```

```
## Warning: Removed 12 rows containing non-finite outside the scale range
## (`stat_binhex()`).
```

```
## Warning: Removed 9 rows containing non-finite outside the scale range
## (`stat_binhex()`).
```

```
## Warning: Removed 12 rows containing non-finite outside the scale range
## (`stat_binhex()`).
```

```
## Warning: Removed 18 rows containing non-finite outside the scale range
## (`stat_binhex()`).
```

![](gating_cells_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

Use the cleaned data for making training and testing sets


``` r
gvhd_tibble <- gvhd_tibble |> mutate(
  initial_split = purrr::map(exprs_tibble, ~ initial_split(.x, strata = `FSC-H`, prop = 0.8)),
  training = purrr::map(initial_split, ~ training(.x)),
  test = purrr::map(initial_split, ~ testing(.x))
)

gvhd_tibble
```

```
## # A tibble: 35 × 6
##    exprs    keywords     exprs_tibble initial_split       training test    
##    <list>   <list>       <list>       <list>              <list>   <list>  
##  1 <dbl[…]> <named list> <tibble>     <split [2195/551]>  <tibble> <tibble>
##  2 <dbl[…]> <named list> <tibble>     <split [1874/471]>  <tibble> <tibble>
##  3 <dbl[…]> <named list> <tibble>     <split [2310/579]>  <tibble> <tibble>
##  4 <dbl[…]> <named list> <tibble>     <split [5258/1316]> <tibble> <tibble>
##  5 <dbl[…]> <named list> <tibble>     <split [6850/1715]> <tibble> <tibble>
##  6 <dbl[…]> <named list> <tibble>     <split [2678/671]>  <tibble> <tibble>
##  7 <dbl[…]> <named list> <tibble>     <split [9838/2461]> <tibble> <tibble>
##  8 <dbl[…]> <named list> <tibble>     <split [1615/407]>  <tibble> <tibble>
##  9 <dbl[…]> <named list> <tibble>     <split [9261/2317]> <tibble> <tibble>
## 10 <dbl[…]> <named list> <tibble>     <split [7984/1998]> <tibble> <tibble>
## # ℹ 25 more rows
```

Now we have a look at the recorded marker expression


``` r
plot_list <- gvhd_tibble$exprs_tibble |> purrr::map(function(.x) {
  .x |> 
    pivot_longer(where(is.numeric)) |>
    ggplot(aes(x = log10(value + 1), y = name)) +
      ggridges::geom_density_ridges() 
    }
  )

wrap_plots(plot_list[1:4])
```

```
## Picking joint bandwidth of 0.109
```

```
## Picking joint bandwidth of 0.0974
```

```
## Picking joint bandwidth of 0.0457
```

```
## Picking joint bandwidth of 0.044
```

![](gating_cells_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

It looks like FL2-A is a clearly separated population, which we can use for 
semi-supervised clustering.

We will try three different methods for gating the FSC-H vs SSC-H:
1) `Mclust`
2) `dbscan`
