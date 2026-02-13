# sequenceRCT

<img src="https://github.com/user-attachments/assets/24ce7601-a82b-4f59-9e15-d3bc143e62b9" width="50%" alt="sequenceRCTicon">

**Sequence Analysis Tools for Randomized Controlled Trials**

[![Version](https://img.shields.io/badge/version-0.2.0-blue.svg)](https://github.com/gav888/sequenceRCT/releases)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

## Overview

`sequenceRCT` is an R package for time-sensitive analysis of behavioral trajectories in randomized controlled trials. Instead of focusing only on endpoint differences, it helps you analyze how participants move across discrete states over time (e.g., relapse, stability, improvement).

The package supports an end-to-end workflow:

- state encoding from raw panel data
- group-wise transition matrix estimation
- sequence complexity profiling (entropy, turbulence, volatility)
- trajectory clustering
- inferential tests for between-group differences
- publication-ready visualizations

## Typical Use Case

Use `sequenceRCT` when your RCT data has:

- one row per participant
- repeated state measurements at ordered time points (e.g., `T1`, `T2`, `T3`, ...)
- a treatment/control grouping variable

This is common in behavioral public policy settings where timing, persistence, and relapse are central to intervention design.

## Main Functions

| Function | Purpose |
|---|---|
| `analyze_rct_sequences()` | Runs the full pipeline and returns sequences, transitions, complexity, tests, clusters, and plots |
| `encode_states()` | Converts raw state labels to coded sequence matrices |
| `calculate_empirical_transitions()` | Estimates transition probabilities between adjacent time points |
| `compute_complexity_measures()` | Computes entropy, turbulence, and volatility for each participant sequence |
| `run_sequence_clustering()` | Clusters trajectories using TraMineR distances + hierarchical clustering |
| `test_complexity_measures()` | Tests group differences in complexity metrics |
| `test_cluster_membership()` | Tests group-cluster association (chi-square) |
| `plot_*()` helpers | Visualizes state distributions, sequence indices, transitions, complexity, and clusters |

## Installation

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("gav888/sequenceRCT")
```

Or with `remotes`:

```r
# install.packages("remotes")
remotes::install_github("gav888/sequenceRCT")
```

## Quick Start

```r
library(sequenceRCT)

# Example RCT dataset with three repeated state measurements
df <- data.frame(
  group = rep(c("Treatment", "Control"), each = 10),
  T1 = sample(c("Lost", "Stable", "Gained"), 20, replace = TRUE),
  T2 = sample(c("Lost", "Stable", "Gained"), 20, replace = TRUE),
  T3 = sample(c("Lost", "Stable", "Gained"), 20, replace = TRUE)
)

# Run the full sequence analysis pipeline
results <- analyze_rct_sequences(
  data = df,
  group_col = "group",
  state_cols = c("T1", "T2", "T3"),
  state_labels = c("Lost", "Stable", "Gained")
)

# Core outputs
head(results$complexity)
results$complexity_tests
results$cluster_tests

# Default plots
print(results$plots$state_distribution)
print(results$plots$sequence_index)
print(results$plots$markov_transitions)
print(results$plots$complexity_measures)
print(results$plots$cluster_membership)
```

## Citation

If you use `sequenceRCT`, please cite:

Veltri GA (2026) *Time-sensitive RCTs in behavioral public policy: a pragmatic framework using sequence methods, personalization, and reinforcement learning*. Frontiers in Behavioral Economics, 5:1684887. [https://doi.org/10.3389/frbhe.2026.1684887](https://doi.org/10.3389/frbhe.2026.1684887)

BibTeX:

```bibtex
@article{veltri2026sequenceRCT,
  author = {Veltri, Giuseppe Alessandro},
  title = {Time-sensitive RCTs in behavioral public policy: a pragmatic framework using sequence methods, personalization, and reinforcement learning},
  journal = {Frontiers in Behavioral Economics},
  year = {2026},
  volume = {5},
  pages = {1684887},
  doi = {10.3389/frbhe.2026.1684887},
  url = {https://doi.org/10.3389/frbhe.2026.1684887}
}
```

## Contributing

Issues and pull requests are welcome at [https://github.com/gav888/sequenceRCT](https://github.com/gav888/sequenceRCT).

## License

This project is licensed under the [GPL-3 License](https://www.gnu.org/licenses/gpl-3.0.en.html).
