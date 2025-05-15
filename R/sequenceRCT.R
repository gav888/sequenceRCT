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
  # Extract the underlying numeric code matrix from the sequence object
  mat <- unclass(seqdata)
  if (!is.matrix(mat)) mat <- as.matrix(mat)
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

#' Encode raw state columns into integer matrix for TraMineR
#'
#' @param data Data frame containing state columns.
#' @param state_cols Character vector; names of state columns to encode.
#' @param state_labels Optional character vector; explicit labels for states.
#' @return A list with elements:
#'   \item{seq_mat}{Integer matrix of encoded states.}
#'   \item{labels}{Character vector of state labels.}
#'   \item{alphabet}{Integer vector of unique codes.}
#' @export
encode_states <- function(data, state_cols, state_labels = NULL) {
  raw <- data[, state_cols, drop = FALSE]
  if (!is.null(state_labels)) {
    # Use provided labels
    seq_mat <- apply(raw, 2, function(x) as.integer(factor(x, levels = state_labels)) - 1)
    labels <- state_labels
    alphabet <- seq_along(state_labels) - 1
  } else {
    # Infer labels from data
    uniq <- sort(unique(as.vector(raw)))
    labels <- as.character(uniq)
    seq_mat <- apply(raw, 2, function(x) as.integer(factor(x, levels = uniq)) - 1)
    alphabet <- uniq
  }
  list(seq_mat = seq_mat, labels = labels, alphabet = alphabet)
}

#' High-level analysis pipeline for RCT sequence data
#'
#' @param data Data frame containing one row per individual, with group and state columns.
#' @param group_col String; name of the column containing group labels (e.g. "Treatment" or "Control").
#' @param state_cols Character vector; names of columns representing states at each time point.
#' @param state_labels Character vector; optional labels for unique states in the order of \code{state_cols}.
#' @param palette Optional color vector for states.
#' @return A list with elements:
#'   \item{seqdata}{A TraMineR sequence object.}
#'   \item{transitions}{Data frame of From, To, Probability, Group for all groups.}
#'   \item{complexity}{Data frame of Entropy, Turbulence, Volatility, and Group.}
#'   \item{clusters}{Integer vector of cluster membership for each individual.}
#' @export
analyze_rct_sequences <- function(data, group_col, state_cols, state_labels = NULL, palette = NULL) {
  # Extract group factor
  group <- as.factor(data[[group_col]])
  # Encode raw state columns robustly
  enc <- encode_states(data, state_cols, state_labels)
  seq_mat <- enc$seq_mat
  labels <- enc$labels
  alphabet <- enc$alphabet
  # Create TraMineR sequence object
  seqdata <- TraMineR::seqdef(seq_mat, alphabet = alphabet, states = labels, labels = labels)
  # Compute empirical transitions by group
  transitions_list <- lapply(levels(group), function(g) {
    mat <- calculate_empirical_transitions(seq_mat[group == g, , drop = FALSE])
    df <- as.data.frame(mat)
    df$From <- labels
    df_long <- reshape2::melt(df, id.vars = "From", variable.name = "To", value.name = "Probability")
    df_long$Group <- g
    df_long
  })
  transitions <- do.call(rbind, transitions_list)
  # Compute complexity measures and add group
  complexity <- compute_complexity_measures(seqdata)
  complexity$Group <- group
  # Perform clustering (default k = 3)
  clusters <- run_sequence_clustering(seqdata)

  # Run statistical tests
  complexity_tests <- test_complexity_measures(complexity)
  cluster_tests    <- test_cluster_membership(clusters, group)

  # Generate default plots
  plots <- list(
    state_distribution   = plot_state_distribution(seqdata, group, palette = palette),
    sequence_index       = plot_sequence_index(seqdata, group, palette = palette),
    markov_transitions   = plot_markov_transitions(transitions, palette = palette),
    complexity_measures  = plot_complexity_measures(complexity),
    cluster_membership   = plot_cluster_membership(clusters, group)
  )

  # Return a list of results
  list(
    seqdata = seqdata,
    transitions = transitions,
    complexity = complexity,
    clusters = clusters,
    complexity_tests = complexity_tests,
    cluster_tests    = cluster_tests,
    plots            = plots
  )
}

#' Format a transition probability matrix for visualization
#'
#' @param trans_mat Numeric matrix; transition probabilities (states × states).
#' @param labels Character vector; state labels corresponding to rows/columns of \code{trans_mat}.
#' @return Data frame with columns \code{From}, \code{To}, and \code{Probability}.
#' @export
format_transitions_for_viz <- function(trans_mat, labels) {
  df <- as.data.frame(trans_mat)
  colnames(df) <- labels
  df$From <- labels
  df_long <- reshape2::melt(
    df,
    id.vars = "From",
    variable.name = "To",
    value.name = "Probability"
  )
  df_long
}

#' Test between-group differences in complexity measures
#'
#' @param complexity_df Data frame; output of \code{compute_complexity_measures()}, must include a \code{Group} column.
#' @param measures Character vector; names of complexity metrics to test (e.g., c("Entropy","Turbulence","Volatility")).
#' @return Named list of lists, each with elements \code{wilcox} (Wilcoxon test object) and \code{effect_size} (Cliff's delta object).
#' @export
test_complexity_measures <- function(complexity_df, measures = c("Entropy", "Turbulence", "Volatility")) {
  results <- lapply(measures, function(m) {
    formula <- stats::as.formula(paste(m, "~ Group"))
    wt <- stats::wilcox.test(formula, data = complexity_df)
    es <- effsize::cliff.delta(formula, data = complexity_df)
    list(wilcox = wt, effect_size = es)
  })
  names(results) <- measures
  results
}

#' Test association between cluster membership and group
#'
#' @param clusters Integer vector; cluster membership for each individual.
#' @param group Factor or character vector; group labels for each individual.
#' @return A list with elements \code{contingency_table} (table) and \code{chi_square} (chisq.test object).
#' @export
test_cluster_membership <- function(clusters, group) {
  tab <- table(Cluster = clusters, Group = group)
  chisq <- stats::chisq.test(tab)
  list(contingency_table = tab, chi_square = chisq)
}

#' Plot complexity measures by group
#'
#' @param complexity_df Data frame returned by \code{compute_complexity_measures()}, must include \code{Group}.
#' @param measures Character vector; complexity metrics to plot (default c("Entropy","Turbulence","Volatility")).
#' @return ggplot object with boxplots and facets for each measure.
#' @export
plot_complexity_measures <- function(complexity_df, measures = c("Entropy", "Turbulence", "Volatility")) {
  long <- reshape2::melt(
    complexity_df,
    id.vars = "Group",
    measure.vars = measures,
    variable.name = "Measure",
    value.name = "Value"
  )
  ggplot(long, aes(x = Group, y = Value, fill = Group)) +
    geom_boxplot() +
    facet_wrap(~ Measure, scales = "free_y") +
    theme_minimal() +
    theme(legend.position = "none")
}

#' Plot cluster membership distribution by group
#'
#' @param clusters Integer vector; output of \code{run_sequence_clustering} or element from \code{analyze_rct_sequences}.
#' @param group Factor or character vector; group labels corresponding to clusters.
#' @return ggplot object of cluster frequencies by group.
#' @export
plot_cluster_membership <- function(clusters, group) {
  df <- data.frame(Cluster = factor(clusters), Group = as.factor(group))
  tab <- as.data.frame(table(df))
  ggplot(tab, aes(x = Cluster, y = Freq, fill = Group)) +
    geom_col(position = "dodge") +
    labs(y = "Count") +
    theme_minimal()
}
