% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sequenceRCT.R
\name{plot_markov_transitions}
\alias{plot_markov_transitions}
\title{Plot Markov transition flows as an alluvial diagram}
\usage{
plot_markov_transitions(transitions, palette = NULL)
}
\arguments{
\item{transitions}{Data frame with columns From, To, Probability, Group.}

\item{palette}{Optional color vector matching levels of From.}
}
\value{
A ggplot object.
}
\description{
Plot Markov transition flows as an alluvial diagram
}
