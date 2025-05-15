# Avoid R CMD check notes for unbound variables in ggplot/alluvial calls
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("Cluster", "Freq", "Group",
      "From", "To", "Probability",
      "Value", "stratum")
  )
}
