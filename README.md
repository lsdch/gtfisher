# Groupwise Truncated Fisher (GTF)

Implements the Groupwise Truncated Fisher (GTF) score to aggregate
groups of p-values from independent tests into a single p-value per group, using an approach derived from Fisher's method.
It is designed as an alternative to Fisher's method, that does not suffer as much
from the non uniformity of p-values under the null hypothesis, 
and is expected to work well as long as there is not an enrichment in very low p-values under the null.
The API is designed to integrate with the package `dplyr`, and optionally the package `progressr`

## Installation

```R
devtools::install_github("https://github.com/lsdch/gtfisher")

# optional 
install.packages("progressr")
```
Computation of scores is fast, but can take a few minutes on large datasets, and takes longer when the simulated sample size is chosen large (see `n` parameter in `gtf_predict()`). 
This packages implements progress signaling using [progressr](https://github.com/HenrikBengtsson/progressr), 
which gives freedom to the user on how progress should be reported. 


## Usage 

```R
library(gtfisher)
library(progressr) # optional

# no progress bar
group_pvals = gtf_predict(pvals, group, p)

# with progress bar
handlers("progress") # pick your favorite one, see `progressr` vignette
with_progress({
  group_pvals = gtf_predict(pvals, group, p)
})
```
