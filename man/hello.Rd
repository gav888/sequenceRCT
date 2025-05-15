
\name{simulate_sequences}
\alias{simulate_sequences}
\title{Simulate Behavioral Sequences Using a Markov Model}
\usage{
simulate_sequences(n_individuals, n_time_points, base_transitions, characteristics, group, dropout_rate = 0, scaling = 1)
}
\arguments{
  \item{n_individuals}{Integer, number of subjects to simulate.}
  \item{n_time_points}{Integer, number of time points per subject.}
  \item{base_transitions}{Numeric matrix, base transition probabilities (states × states).}
  \item{characteristics}{Data frame, one row per individual containing columns \code{responsiveness} and \code{prior_engagement}.}
  \item{group}{Character vector, group assignment per individual (e.g., "Treatment" or "Control").}
  \item{dropout_rate}{Numeric between 0 and 1, probability of dropout at each time point.}
  \item{scaling}{Numeric, treatment effect scaling factor applied for "Treatment" group.}
}
\description{
Simulate behavior sequences for individuals based on a Markov model with optional dropout and treatment effects.
}
\details{
The function returns an integer matrix of size \code{n_individuals} by \code{n_time_points}, where each entry is a state code (0-based). For the "Treatment" group, transition probabilities are modified by individual responsiveness and an exponential time decay scaled by \code{scaling}. All individuals also experience an engagement effect based on \code{prior_engagement}. Dropout can occur at a random time point with probability \code{dropout_rate}.
}
\examples{
# Example usage
set.seed(42)
char_df <- data.frame(responsiveness = rlnorm(5, 0, 0.5), prior_engagement = runif(5))
base_mat <- matrix(c(0.7,0.2,0.1, 0.2,0.6,0.2, 0.1,0.2,0.7), nrow = 3, byrow = TRUE)
group <- c("Treatment", "Control", "Treatment", "Control", "Treatment")
seqs <- simulate_sequences(5, 10, base_mat, char_df, group, dropout_rate = 0.1, scaling = 1.2)
}
\keyword{package}
