# sequenceRCT

**Sequence Analysis Tools for Randomized Controlled Trials**

[![Version](https://img.shields.io/badge/version-0.1.0-blue.svg)](https://github.com/gav888/sequenceRCT/releases)  
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

## Overview

`sequenceRCT` provides agnostic functions to analyze behavioral sequence data from randomized controlled trials. It includes utilities to:

- Encode raw state data (`encode_states`)
- Estimate empirical transition probabilities (`calculate_empirical_transitions`)
- Compute complexity measures: entropy, turbulence, volatility (`compute_complexity_measures`)
- Perform hierarchical clustering on sequences (`run_sequence_clustering`)
- Run statistical tests for group differences (Wilcoxon, Cliff’s delta, Chi-square) (`test_complexity_measures`, `test_cluster_membership`)
- Format and visualize results:
  - Sequence plots (`plot_state_distribution`, `plot_sequence_index`)
  - Alluvial transition diagrams (`plot_markov_transitions`)
  - Complexity boxplots (`plot_complexity_measures`)
  - Cluster membership bar charts (`plot_cluster_membership`)
- A single high-level pipeline function (`analyze_rct_sequences`) that ties everything together.

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

## Usage

```r
library(sequenceRCT)

# Prepare a sample RCT dataset
df <- data.frame(
  group = rep(c("Treatment", "Control"), each = 10),
  T1 = sample(c("Lost", "Stable", "Gained"), 20, replace = TRUE),
  T2 = sample(c("Lost", "Stable", "Gained"), 20, replace = TRUE),
  T3 = sample(c("Lost", "Stable", "Gained"), 20, replace = TRUE)
)

# Run the full analysis pipeline
results <- analyze_rct_sequences(
  data = df,
  group_col = "group",
  state_cols = c("T1", "T2", "T3"),
  state_labels = c("Lost", "Stable", "Gained")
)

# View complexity measures
head(results$complexity)

# Display default plots
print(results$plots$state_distribution)
print(results$plots$sequence_index)
print(results$plots$markov_transitions)
print(results$plots$complexity_measures)
print(results$plots$cluster_membership)
```

## Contributing

Please file issues or pull requests at [https://github.com/gav888/sequenceRCT](https://github.com/gav888/sequenceRCT).

## License

This project is licensed under the [GPL-3 License](https://www.gnu.org/licenses/gpl-3.0.en.html).
