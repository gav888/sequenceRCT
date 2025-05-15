
#' Simulate behavior sequences using a Markov model with optional dropout
#'
#' @param n_individuals Integer; number of subjects to simulate.
#' @param n_time_points Integer; number of time points per subject.
#' @param base_transitions Numeric matrix; base transition probabilities (states × states).
#' @param characteristics Data frame; one row per individual containing named columns `responsiveness` and `prior_engagement`.
#' @param group Character vector of length n_individuals; group assignment per individual (e.g., "Treatment" or "Control").
#' @param dropout_rate Numeric between 0 and 1; probability of dropout at each time point.
#' @param scaling Numeric; treatment effect scaling factor (only applied for "Treatment" group).
#' @return Integer matrix of dimension n_individuals × n_time_points with state codes.
#' @export
simulate_sequences <- function(n_individuals,
                               n_time_points,
                               base_transitions,
                               characteristics,
                               group,
                               dropout_rate = 0,
                               scaling = 1) {
  seq_mat <- matrix(NA_integer_, nrow = n_individuals, ncol = n_time_points)
  for (i in seq_len(n_individuals)) {
    init_probs <- c(0.6, 0.3, 0.1)
    if (characteristics$prior_engagement[i] > 0.7) {
      init_probs <- c(0.3, 0.4, 0.3)
    }
    seq_mat[i, 1] <- sample(seq_along(init_probs) - 1, 1, prob = init_probs)
    dropout_time <- if (runif(1) < dropout_rate) sample(2:n_time_points, 1) else n_time_points + 1
    for (t in 2:min(dropout_time, n_time_points)) {
      indiv <- list(
        responsiveness = characteristics$responsiveness[i],
        prior_engagement = characteristics$prior_engagement[i]
      )
      trans <- base_transitions
      if (group[i] == "Treatment") {
        te <- indiv$responsiveness * exp(-0.2 * t) * scaling
        trans[1,2] <- trans[1,2] + te * 0.1
        trans[2,3] <- trans[2,3] + te * 0.15
        trans <- t(apply(trans, 1, function(x) x / sum(x)))
      }
      eng_eff <- indiv$prior_engagement * 0.2
      trans[,3] <- trans[,3] + eng_eff
      trans <- t(apply(trans, 1, function(x) x / sum(x)))
      prev_state <- seq_mat[i, t - 1] + 1
      seq_mat[i, t] <- sample(seq_len(ncol(trans)) - 1, 1, prob = trans[prev_state, ])
    }
  }
  seq_mat
}

#' Compute empirical transition probability matrix from sequences
#'
#' @param sequences Integer matrix; rows = individuals, cols = time points, values = state codes starting at 0.
#' @return Numeric matrix of size n_states × n_states with estimated transition probabilities.
#' @export
calculate_empirical_transitions <- function(sequences) {
  n_states <- length(unique(as.vector(sequences)))
  tm <- matrix(0, nrow = n_states, ncol = n_states)
  for (i in seq_len(nrow(sequences))) {
    for (t in seq_len(ncol(sequences) - 1)) {
      from <- sequences[i, t] + 1
      to   <- sequences[i, t + 1] + 1
      if (!is.na(from) && !is.na(to)) tm[from, to] <- tm[from, to] + 1
    }
  }
  tm / rowSums(tm)
}

#' Plot state distribution by group
#'
#' @param seqdata TraMineR sequence object created via \code{seqdef}.
#' @param group Factor or character vector indicating group for each sequence.
#' @param palette Optional color vector for states.
#' @param ... Additional args passed to \code{seqdplot}.
#' @export
plot_state_distribution <- function(seqdata, group, palette = NULL, ...) {
  seqdplot(
    seqdata,
    group = group,
    with.legend = TRUE,
    cpal = palette,
    ...
  )
}

#' Plot sequence index by group
#'
#' @param seqdata TraMineR sequence object.
#' @param group Group labels.
#' @param palette Optional palette.
#' @param ... Additional args for \code{seqiplot}.
#' @export
plot_sequence_index <- function(seqdata, group, palette = NULL, ...) {
  seqiplot(
    seqdata,
    group = group,
    with.legend = TRUE,
    cpal = palette,
    ...
  )
}

#' Plot Markov transition flows as an alluvial diagram
#'
#' @param transitions Data frame with columns From, To, Probability, Group.
#' @param palette Optional color vector matching levels of From.
#' @return A ggplot object.
#' @export
plot_markov_transitions <- function(transitions, palette = NULL) {
  ggplot(transitions,
         aes(axis1 = From, axis2 = To, y = Probability, fill = From)) +
    geom_alluvium(aes(alpha = Probability)) +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    facet_wrap(~ Group) +
    scale_fill_manual(values = palette) +
    theme_minimal() +
    theme(legend.position = "none")
}

#' Compute complexity measures: entropy, turbulence, volatility
#'
#' @param seqdata TraMineR sequence object.
#' @return Data frame with Entropy, Turbulence, Volatility for each sequence.
#' @export
compute_complexity_measures <- function(seqdata) {
  mat <- as.matrix(seqdata)
  n_time <- ncol(mat)
  data.frame(
    Entropy    = TraMineR::seqient(seqdata),
    Turbulence = TraMineR::seqST(seqdata),
    Volatility  = apply(mat, 1, function(x) sum(diff(x) != 0, na.rm = TRUE) / (n_time - 1))
  )
}

#' Perform hierarchical clustering on sequences
#'
#' @param seqdata TraMineR sequence object.
#' @param k Integer; number of clusters.
#' @param method Character; clustering method (default "ward.D2").
#' @param distance Character; seqdist method (default "LCS").
#' @return Integer vector of cluster membership.
#' @export
run_sequence_clustering <- function(seqdata, k = 3, method = "ward.D2", distance = "LCS") {
  distm <- TraMineR::seqdist(seqdata, method = distance)
  clust <- stats::hclust(as.dist(distm), method = method)
  cutree(clust, k = k)
}
