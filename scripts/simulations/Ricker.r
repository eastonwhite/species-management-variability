# Deterministic Ricker model for local growth
#
Ricker <- function(Nt, R, alpha){
  Nt * R * exp(-1 * alpha * Nt)
}
