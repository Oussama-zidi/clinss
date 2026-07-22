# clinss

<!-- badges: start -->
[![R CMD check](https://github.com/YOURUSERNAME/clinss/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/YOURUSERNAME/clinss/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Sample size and power calculations for clinical trials, in base R.

An open-source alternative to PASS, organized the way trial statisticians
think: by endpoint (continuous, binary, survival, ...), design (parallel,
paired, crossover, cluster-randomized, ...), and hypothesis (superiority,
non-inferiority, superiority by a margin, equivalence).

## Design principles

1. **Base R only.** No runtime dependencies beyond `stats`, to maximize the
   R Validation Hub / `riskmetric` risk score and to keep the package easy
   to qualify in regulated environments. Every formula is implemented from
   scratch.
2. **Primary references, not software.** Every procedure cites the
   published formula it implements (Julious; Chow, Shao & Wang; Machin et
   al.; ...).
3. **Exact power, then search.** Each procedure has a `power_*()` engine
   (exact where feasible, e.g., noncentral t) and an `ss_*()` wrapper that
   finds the smallest integer sample size reaching the target power,
   seeded by the closed-form normal approximation.
4. **Rich, self-documenting results.** Every `ss_*()` function returns a
   `clinss_result` object with sample sizes, achieved power, assumptions,
   method, and reference, plus `print()`, `summary()`, `as.data.frame()`,
   and `report()` (protocol-ready text) methods.
5. **Validation is a feature.** Every procedure is tested against an
   independent implementation and/or published closed-form values. See
   `inst/validation/`.

## Installation

```r
# From the built source tarball:
install.packages("clinss_0.2.0.tar.gz", repos = NULL, type = "source")
```

## Quick start

```r
library(clinss)

# Superiority: detect a 5-unit difference (SD 10), 90% power, alpha 0.05
res <- ss_two_means(effect = 5, sd = 10, power = 0.90, alpha = 0.05,
                    sides = 2, dropout = 0.10)
summary(res)
report(res)

# Non-inferiority
ss_two_means(effect = 0, sd = 20, margin = 10, alpha = 0.025,
             power = 0.90, hypothesis = "noninferiority")

# Equivalence (TOST), 2:1 allocation
ss_two_means(effect = 0, sd = 10, margin = 5, alpha = 0.05, power = 0.80,
             hypothesis = "equivalence", ratio = 2)

# Binary endpoint: superiority, 60% vs 40% response
ss_two_proportions(p1 = 0.60, p2 = 0.40, power = 0.90, alpha = 0.05,
                   sides = 2)

# Binary non-inferiority (Farrington-Manning score test)
ss_two_proportions(p1 = 0.85, p2 = 0.85, margin = 0.10, alpha = 0.025,
                   power = 0.90, hypothesis = "noninferiority")

# One proportion, exact binomial test
ss_one_proportion(p = 0.35, p0 = 0.20, power = 0.80, alpha = 0.05,
                  sides = 1)

# Sensitivity table over a range of standard deviations
do.call(rbind, lapply(c(8, 10, 12), function(s)
  as.data.frame(ss_two_means(effect = 5, sd = s, power = 0.90,
                             alpha = 0.05, sides = 2))))
```

## Package layout

```
R/
  core-validation.R          input checks shared by all procedures
  core-quantiles.R           z_alpha(), z_beta(), dropout, allocation
  core-result.R              clinss_result class + print/summary/as.data.frame
  core-clinss-result-class.R clinss_result help page
  core-report.R              report() protocol text
  continuous-two-means.R     two independent means (all 4 hypotheses)
  continuous-one-mean.R      one mean and paired means
  binary-two-proportions.R   two proportions (pooled z / Farrington-Manning)
  binary-one-proportion.R    one proportion (exact binomial / z)

.github/
  workflows/
    R-CMD-check.yaml         CI: build + check on Windows and Ubuntu
```

## Development workflow

```sh
# In R (after every edit to a source file):
roxygen2::roxygenise()      # regenerate man/ and NAMESPACE

# In a terminal:
R CMD build clinss
R CMD check clinss_0.2.0.tar.gz

# Then push:
git add .
git commit -m "Add ..."
git push
```

## Roadmap

| Version | Scope |
|---------|-------|
| 0.1 | Core utilities + continuous endpoints (one, two, paired means), all four hypothesis types |
| 0.2 | Binary endpoints (two proportions: chi-square/z, Farrington-Manning NI, Fisher exact, one proportion) |
| 0.3 | Survival (Schoenfeld and Freedman log-rank, exponential, accrual/follow-up) |
| 0.4 | ANCOVA, regression (linear, logistic, Cox, Poisson) |
| 0.5 | Cluster-randomized (design effect, unequal clusters, stepped wedge) and crossover |
| 0.6 | Group sequential (Lan-DeMets spending functions, O'Brien-Fleming, Pocock) |
| 0.7 | Plotting (power curves, sensitivity analyses) |
| 1.0 | Stable API, full validation dossier, CRAN submission |
