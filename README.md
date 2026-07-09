# AncCond: The Ancestral Condition Test

This repository contains a modern, dependency-light implementation of the **Ancestral
Condition (AncCond) test** together with all of the code, data, and results needed to
reproduce the power analyses, adversarial simulations, and empirical squamate example
from the accompanying manuscript.

The AncCond test asks a focused evolutionary question:

> **Do transitions in a binary discrete trait tend to occur when a co-evolving
> continuous trait sits at an extreme (high or low) value?**

A classic motivating example is the evolution of viviparity (live birth) in reptiles:
does the switch from egg-laying to live birth tend to arise in lineages that have
already evolved a particular body size? AncCond reconstructs the continuous trait along
the tree, locates where the discrete trait changes state, and asks whether those change
points fall at unusually high or low values of the continuous trait relative to a null
model in which the discrete trait evolves independently.

---

## Table of contents

- [Overview](#overview)
- [How the method works](#how-the-method-works)
- [Repository structure](#repository-structure)
- [Installation and requirements](#installation-and-requirements)
- [Quick start](#quick-start)
- [The `AncCond()` function](#the-anccond-function)
- [Reproducible analyses](#reproducible-analyses)
- [Data](#data)
- [Figures](#figures)
- [Results](#results)
- [Working directories and paths](#working-directories-and-paths)
- [Citation and authors](#citation-and-authors)
- [Contact](#contact)

---

## Overview

The AncCond test was originally introduced by the Blackmon lab and distributed in the
[`evobiR`](https://rdrr.io/github/coleoguy/evobir/) R package. This repository provides a
**rewritten, self-contained implementation** (`R/AncCond.R`) that:

- depends only on `ape` and `phytools` (the latter used solely for `fastAnc`);
- fits the discrete-trait model internally with a fast, log-rescaled pruning algorithm
  rather than calling out to a general Mk fitter;
- tests **both** directions of change (gains, `0 → 1`, and losses, `1 → 0`) in a single run;
- supports one-sided and two-sided hypotheses, an optional user-supplied rate matrix, and
  a "pruned" reconstruction mode that estimates the continuous trait using only lineages
  still in the ancestral state.

The repository is organized as a **reproducibility project**: the `R/` module is the
method itself, and the `analyses/`, `figures/`, `data/`, and `results/` directories let a
reader regenerate every number and figure in the manuscript.

## How the method works

For a rooted phylogeny with a continuous trait `x` and a binary trait `y` at the tips,
`AncCond()` proceeds in five conceptual steps:

1. **Continuous-trait reconstruction.** Ancestral values of `x` at internal nodes are
   estimated with maximum likelihood (`phytools::fastAnc`). In `prune = TRUE` mode the
   reconstruction uses only the tips that remain in state `0`, so that selection acting
   after the discrete transition does not bias the ancestral estimates.
2. **Discrete-trait model.** A two-state Markov (Mk) model is fitted to `y`. The code
   compares an all-rates-different (ARD) model against a unidirectional model (`q10 = 0`)
   and keeps the AIC-preferred fit (conservatively favouring the simpler model). A fixed
   rate matrix `Q` may be supplied to skip model selection entirely.
3. **Observed test statistic.** Discrete-trait histories are sampled by stochastic
   mapping. On every branch that changes state, the transition point is placed by drawing
   from a truncated exponential and the continuous trait is linearly interpolated between
   the parent and child reconstructed values. The statistic `T` is the mean interpolated
   value at transition points, computed separately for `0 → 1` and `1 → 0` changes.
4. **Null distribution.** Holding the continuous trait fixed, many discrete-trait
   datasets are simulated under the fitted Mk model, and `T` is recomputed for each,
   yielding a null distribution for gains and for losses.
5. **P-values.** The observed `T` is compared to the null with a continuity correction,
   producing upper- and lower-tail p-values, a hypothesis-adjusted p-value per direction,
   and a Bonferroni-combined p-value across both directions.

## Repository structure

```
anc.cond/
├── R/
│   ├── AncCond.R            # Core method: AncCond() + print method + internal helpers
│   ├── AncCond.Rd           # R help page for AncCond()
│
├── analyses/
│   ├── 01_sim_power_analysis.R      # Power & false-positive rate (symmetric Mk)
│   ├── 02_sim_OU_adversarial.R      # Adversarial "causation reversal" (state-dependent OU)
│   ├── 03_sim_asymm_comparision.R   # Power & FPR under asymmetric transition rates
│   ├── 04_empirical_squamates.R     # Empirical example: body size vs. viviparity
|
├── data/
│   ├── output.nex                   # 100 posterior squamate trees (Tonini et al. 2016)
│   ├── squamate_traits.csv          # Body size (log SVL) + reproductive mode per species
│
├── figures/
│   ├── *.pdf                        # Manuscript figures (see Figures section)
│   ├── figure_scripts/              # Scripts that build each figure from result CSVs
│
├── results/
│   ├── 01_sim_power_*.csv           # Symmetric power/FPR raw + summary tables
│   ├── 01_sim_asymm_power_*.csv     # Asymmetric power/FPR raw + summary tables
│   ├── 02_sim_OU_adversarial_*      # Adversarial FPR results (.csv and .rds)
│   ├── empirical_squamate_results.csv  # Per-tree AncCond results for squamates
│
└── README.md                        # This file
```

## Installation and requirements

The code targets **R ≥ 4.0**. There is no `DESCRIPTION`/`NAMESPACE`; this is a script
and data repository rather than an installable package, so you use it by cloning the repo
and `source()`-ing the function.

```bash
git clone https://github.com/coleoguy/anc.cond.git
cd anc.cond
```

Install the R packages used across the method, the simulations, and the figures:

```r
# Core method (required for AncCond itself)
install.packages(c("ape", "phytools"))

# Additional packages used by the analysis scripts
install.packages(c("phylolm", "nlme", "R.utils", "parallel"))

# Additional packages used by the figure scripts
install.packages(c("ggplot2", "dplyr", "tidyr"))
```

| Package     | Used by                          | Purpose                                             |
| ----------- | -------------------------------- | --------------------------------------------------- |
| `ape`       | method + all analyses            | Core phylogenetic data structures                   |
| `phytools`  | method + all analyses            | `fastAnc`, `fastBM`, `pbtree`, `sim.history`        |
| `phylolm`   | analyses 01, 03                  | `phyloglm` comparison method                        |
| `nlme`      | analyses 01, 03                  | `gls` / PGLS comparison method                      |
| `parallel`  | analyses 01–04                   | Cross-platform socket clusters                      |
| `R.utils`   | analysis 02                      | `withTimeout` guard for long replicates             |
| `ggplot2`   | figure scripts                   | Plotting                                            |
| `dplyr`     | figure scripts                   | Data wrangling                                      |
| `tidyr`     | figure scripts                   | Reshaping (`pivot_longer`/`pivot_wider`)            |

## Quick start

```r
library(ape)
library(phytools)
source("R/AncCond.R")

# Simulate a toy dataset
set.seed(1)
tree <- pbtree(n = 100, scale = 1)
x <- fastBM(tree)                                   # continuous trait
y <- rbinom(100, 1, 0.3); names(y) <- tree$tip.label # binary trait

# Run the test (small settings for a quick demo)
res <- AncCond(tree, x, y, n.maps = 50, n.sims = 500)

res                 # pretty-printed summary (print.anccond)
res$p.combined      # Bonferroni p-value across both directions
res$p.01            # p-value for gains (0 -> 1)
res$p.10            # p-value for losses (1 -> 0)
```

For a directional hypothesis (e.g. "gains happen at high continuous values"):

```r
res <- AncCond(tree, x, y, hypothesis = "higher")
res$p.upper.01
```

## The `AncCond()` function

```r
AncCond(tree, x, y,
        n.maps     = 100,
        n.sims     = 1000,
        prune      = FALSE,
        Q          = NULL,
        hypothesis = "none")
```

### Arguments

| Argument     | Type            | Default  | Description                                                                                                   |
| ------------ | --------------- | -------- | ------------------------------------------------------------------------------------------------------------- |
| `tree`       | `phylo`         | —        | A rooted, bifurcating phylogeny with branch lengths.                                                          |
| `x`          | named numeric   | —        | Continuous trait values at the tips; names must match `tree$tip.label`.                                       |
| `y`          | named integer   | —        | Binary discrete trait (`0`/`1`) at the tips; names must match `tree$tip.label`.                               |
| `n.maps`     | integer         | `100`    | Stochastic maps drawn per test-statistic computation.                                                         |
| `n.sims`     | integer         | `1000`   | Null-model simulations used to build the reference distribution.                                              |
| `prune`      | logical         | `FALSE`  | If `TRUE`, reconstruct the continuous trait using only state-`0` lineages (guards against post-transition selection). |
| `Q`          | 2×2 matrix      | `NULL`   | Optional fixed rate matrix (rows/cols labelled `"0"`,`"1"`). If supplied, model selection is skipped.         |
| `hypothesis` | character       | `"none"` | `"none"` (two-sided), `"higher"` (transitions at high values), or `"lower"` (transitions at low values).      |

> **Note on `Q`.** The function expects a 2×2 matrix. The `.Rd` example shows the legacy
> `Q = c(0,2,1,0)` vector form; when calling the current `AncCond()` directly, pass a
> matrix, e.g. `matrix(c(-q01, q01, q10, -q10), 2, 2, byrow = TRUE)`, or leave `Q = NULL`
> to let the function fit the Mk model.

### Return value

An object of class `"anccond"` (a list) with, among others:

| Element                      | Meaning                                                                    |
| ---------------------------- | -------------------------------------------------------------------------- |
| `T.obs.01`, `T.obs.10`       | Observed mean continuous value at `0→1` and `1→0` transitions.             |
| `p.01`, `p.10`               | Hypothesis-adjusted p-values for gains and losses.                         |
| `p.combined`                 | Bonferroni-corrected p-value across both directions: `min(1, 2·min(p.01, p.10))`. |
| `p.upper.01`, `p.lower.01`   | One-sided tail probabilities for `0→1` (`P(T_null ≥ T_obs)`, `P(T_null ≤ T_obs)`). |
| `p.upper.10`, `p.lower.10`   | One-sided tail probabilities for `1→0`.                                    |
| `T.null.01`, `T.null.10`     | Null distributions (length `n.sims`) for each direction.                   |
| `Q`                          | The fitted (or supplied) rate matrix.                                      |
| `model`                      | `"ARD"`, `"unidirectional"`, or `"user-supplied"`.                         |
| `hypothesis`                 | The hypothesis direction used.                                             |

A `print` method (`print.anccond`) renders a compact, human-readable summary of the model,
observed statistics, and p-values.

### Internal helpers

`R/AncCond.R` is fully self-contained. The exported `AncCond()` and `print.anccond()` are
built on internal functions that a curious reader may find useful:

- `.tp2()` — analytical 2-state transition probability matrix.
- `.pruning2()` — Felsenstein pruning with log-rescaling for numerical stability.
- `.fitMk2()` — fits ARD vs. unidirectional Mk models and returns the AIC-preferred fit.
- `.simMk2()` — vectorized batch simulation of discrete tip states under Mk.
- `.anccond.T()` — samples stochastic maps and computes the transition-point statistics.
- `.map.pruned()` — maps pruned-tree ancestral estimates back onto the full tree.

## Reproducible analyses

All four analysis scripts live in `analyses/`, use socket clusters for parallelism
(portable across macOS/Linux/Windows), and set seeds for reproducibility. Each `source()`s
the method file and writes tabular results. See
[Working directories and paths](#working-directories-and-paths) for the exact setup the
scripts assume.

### `01_sim_power_analysis.R` — power and false-positive rates

Simulates a symmetric Mk discrete trait alongside a Brownian continuous trait, then injects
a controlled association by scaling branch lengths in the extreme quartiles of the
continuous trait by a factor `s`. Compares **AncCond** against **`phyloglm`** and **PGLS**.

- Grid: `n.tips ∈ {50, 100, 200, 500}` × `q ∈ {0.1, 0.5, 1.0, 2.0}` × `s ∈ {1, 2, 3, 5, 10}`
  (`s = 1` is the null / false-positive case), with 200 replicates per cell.
- Datasets require ≥ 10% of tips in each state; over-imbalanced draws are rejected and logged.
- **Outputs:** `01_sim_power_raw.csv`, `01_sim_power_results.csv`, `01_sim_discard_log.csv`.

```bash
Rscript 01_sim_power_analysis.R
```

### `02_sim_OU_adversarial.R` — adversarial "causation reversal" test

A stress test of specificity. Here the **discrete trait drives the continuous trait** (the
opposite of AncCond's assumption): the continuous trait evolves under a two-optimum
Ornstein–Uhlenbeck process whose optimum switches with the current discrete state. A
well-behaved test should *not* raise many false positives.

- 200 tips, 1000 replicates; sweeps OU pull-back `alpha ∈ {0.5, 2, 8}` and optimum gap
  `delta.theta ∈ {2, 4, 8}` (`sigma = 1`).
- Runs both the standard and pruned variants of AncCond, with per-replicate timeouts.
- **Outputs:** `02_sim_OU_adversarial_results.rds`, `02_sim_OU_adversarial_results.csv`,
  `02_sim_OU_adversarial_raw.csv`.

### `03_sim_asymm_comparision.R` — asymmetric transition rates

Mirrors script 01 but under **asymmetric** discrete evolution: the reverse rate is held low
(`q10 = 0.1`) while the forward rate varies (`q01 ∈ {0.1, 0.5, 1.0, 2.0}`), with the root
fixed to state `0`. This probes how directionality in the discrete trait affects power and
false-positive control.

- Grid: `n.tips ∈ {50, 100, 200, 500}` × `q01 ∈ {0.1, 0.5, 1.0, 2.0}` × `s ∈ {1, 2, 3, 5, 10}`,
  200 replicates per cell; results reshaped to long form with per-method significance counts.
- **Outputs:** `01_sim_asymm_power_raw.csv`, `01_sim_asymm_power_results.csv`.

### `04_empirical_squamates.R` — body size and viviparity in squamates

Applies AncCond to a real question: **is the evolution of viviparity (live birth) in
squamate reptiles associated with body size?** The continuous trait is `log` maximum
snout–vent length (SVL); the discrete trait is oviparous (`0`, ancestral) vs. viviparous
(`1`, derived).

- Reads the 100 posterior trees in `../data/output.nex` and the traits in
  `../data/squamate_traits.csv`, deduplicates species (mean SVL), and matches the shared
  taxa. Runs AncCond on every posterior tree in parallel (`n.maps = 100`, `n.sims = 500`).
- **Outputs:** `../results/empirical_squamate_results.csv` (per-tree statistics and
  p-values) and `../figures/figure6_empirical.pdf` (null distributions with observed
  values overlaid).
- In the committed results, none of the 100 posterior trees yield `p.combined < 0.05`,
  i.e. this analysis does not detect a significant association between body size and
  reproductive-mode transitions across the posterior. (Consult the manuscript for
  interpretation.)

> **Heads-up:** `04_empirical_squamates.R` sets `n.cores <- 28`. Lower this to match your
> machine (e.g. `parallel::detectCores() - 2`) before running.

## Data

| File                        | Description                                                                                                  |
| --------------------------- | ------------------------------------------------------------------------------------------------------------ |
| `data/output.nex`           | 100 posterior squamate phylogenies subsampled/pruned from [VertLife](http://vertlife.org), derived from Tonini et al. (2016). |
| `data/squamate_traits.csv`  | One row per species record with columns `species`, `log_svl_mm` (log maximum SVL in mm), and `repro_mode` (`Oviparous`/`Viviparous`). ~5,955 records; body-size/reproductive-mode data from SquamBase (Meiri 2024). |


## Figures

Manuscript-ready figures are stored as PDFs in `figures/`, and the scripts that build them
from the result CSVs live in `figures/figure_scripts/`:

| Figure script                          | Reads                                          | Produces                            |
| -------------------------------------- | ---------------------------------------------- | ----------------------------------- |
| `01_sim_power_visualization.R`         | `01_sim_power_results.csv`                     | Power curves (AncCond vs. phyloglm) |
| `01_head_to_head_visualization.R`      | `01_sim_power_results.csv`                     | Head-to-head power scatter          |
| `01_sim_asymm_visualization.R`         | `01_sim_power_results.csv`, `01_sim_asymm_power_results.csv` | Symmetric vs. asymmetric power |
| `02_sim_adversarial_visualization.R`   | `02_sim_OU_adversarial_results.csv`            | Adversarial false-positive rates    |
| `supp_symmetry_check.R`                | `01_sim_power_results.csv`                     | `0→1` vs. `1→0` symmetry check (CCC) |

Rendered outputs include `Method Performance.pdf`, `Asym power comparision.pdf`,
`OU model.pdf`, `Figure 1.pdf`, `phylo_vs_AncCond.pdf`, `AncCond_Symmetry_Check.pdf`, and
`figure6_empirical.pdf`.

## Results

The `results/` directory holds the tabular outputs consumed by the figure scripts:

| File                                   | Contents                                                                 |
| -------------------------------------- | ------------------------------------------------------------------------- |
| `01_sim_power_raw.csv`                 | Per-replicate p-values (AncCond `0→1`/`1→0`, phyloglm, PGLS).             |
| `01_sim_power_results.csv`             | Per-cell rejection rates (`rate`), significant counts, `metric` (FPR/power). |
| `01_sim_discard_log.csv`               | Rejected-draw counts per grid cell.                                       |
| `01_sim_asymm_power_raw.csv` / `_results.csv` | Raw and summarized asymmetric-rate results.                        |
| `02_sim_OU_adversarial_results.csv` / `.rds` | Adversarial false-positive summary (standard vs. pruned).          |
| `02_sim_OU_adversarial_raw.csv`        | Per-replicate adversarial p-values.                                       |
| `empirical_squamate_results.csv`       | Per-tree AncCond statistics and p-values for the squamate example.        |


## Working directories and paths

The scripts were written to run from within the folder that contains them, and they refer
to `AncCond.R` **by bare filename** (e.g. `source("AncCond.R")`) while writing data/figures
via relative paths such as `../data/` and `../results/`. Because the method now lives in
`R/`, do one of the following before running an analysis:

**Option A — copy the method next to the scripts:**

```bash
cp R/AncCond.R analyses/
cd analyses
Rscript 01_sim_power_analysis.R
```

**Option B — edit the `source()` line** in each script to point at `../R/AncCond.R`, then
run from `analyses/`.

Notes:

- Scripts `01`–`03` write their CSV outputs to the **current working directory**; move the
  files into `results/` (as in this repo) to keep things organized, or run the scripts from
  within `results/`.
- `04_empirical_squamates.R` already writes to `../results/` and `../figures/`, so run it
  from `analyses/`.
- Figure scripts read result CSVs **by bare filename**, so run them from `results/` (or copy
  the relevant CSVs into `figures/figure_scripts/` first).

## Citation and authors

If you use the AncCond test, please cite the accompanying manuscript and the original method as distributed in
the `evobiR` package.

**Authors:** Meghann McConnell, Nathan Anderson, Richard H. Adams, and Heath Blackmon.

Lab website: [coleoguy.github.io](http://coleoguy.github.io/)

## Contact

Questions about the method, the simulations, or the empirical example are best directed to
the corresponding authors listed in the manuscript. Issues and pull requests are welcome on the
[GitHub repository](https://github.com/coleoguy/anc.cond).
