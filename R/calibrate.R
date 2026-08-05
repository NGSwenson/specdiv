# Probability calibration ----------------------------------------------------

#' Turn classifier scores into a probability matrix
#'
#' Clamp negative scores to zero and row-normalise, so each row sums to 1. Rows
#' that are all non-positive are made uniform. Useful for turning raw PLS-DA
#' predicted-response scores into probabilities.
#'
#' @param scores Numeric matrix, rows = individuals, columns = species.
#' @return A probability matrix the same shape as `scores`, dimnames kept.
#' @examples
#' s <- matrix(c(0.2, -0.1, 0.9, 0.3), 2, 2)
#' scores_to_prob(s)
#' @export
scores_to_prob <- function(scores) {
  scores <- as.matrix(scores)
  P <- pmax(scores, 0)
  rs <- rowSums(P)
  P[rs == 0, ] <- 1
  out <- P / rowSums(P)
  dimnames(out) <- dimnames(scores)
  out
}

#' Temperature (power) scaling of a probability matrix
#'
#' Sharpen or soften a probability matrix by raising each entry to the power
#' `gamma` and renormalising each row. `gamma > 1` sharpens (more confident),
#' `gamma < 1` softens, `gamma = 1` is a no-op.
#'
#' @param P Probability matrix (rows = individuals, columns = species).
#' @param gamma Positive scalar temperature exponent.
#' @return A probability matrix the same shape as `P`, dimnames kept.
#' @seealso [fit_temperature()] to choose `gamma`.
#' @examples
#' P <- matrix(c(0.6, 0.4, 0.55, 0.45), 2, 2, byrow = TRUE)
#' temperature_scale(P, 4)
#' @export
temperature_scale <- function(P, gamma) {
  P <- as.matrix(P)
  stopifnot(length(gamma) == 1, gamma > 0)
  Pg <- P^gamma
  out <- Pg / rowSums(Pg)
  dimnames(out) <- dimnames(P)
  out
}

#' Fit a temperature so mean confidence matches accuracy
#'
#' Choose the power-scaling exponent `gamma` (see [temperature_scale()]) so that
#' the mean maximum probability (mean confidence) of held-out predictions equals
#' their accuracy. This is one-parameter calibration: honest probabilities are
#' required because [propagate_diversity()] samples from them.
#'
#' Supply cross-validated probabilities and the true labels for the same rows.
#' If the held-out set does not resemble deployment (e.g. reference
#' cross-validation accuracy far exceeds deployment accuracy), the result will
#' be mildly overconfident and downstream intervals too tight; calibrate on a
#' deployment-like set where possible.
#'
#' @param cv_probs Cross-validated probability matrix (rows = held-out
#'   individuals, columns = species; column names = species labels).
#' @param cv_labels Character/factor vector of true labels for those rows.
#' @param grid Numeric vector of candidate `gamma` values to search
#'   (default `seq(1, 12, 0.5)`).
#' @return A list with `gamma` (chosen exponent), `cv_accuracy`, and
#'   `mean_confidence` (achieved mean max-probability at `gamma`).
#' @examples
#' set.seed(1)
#' P <- matrix(runif(30), 10, 3); P <- P / rowSums(P)
#' colnames(P) <- c("a", "b", "c")
#' y <- colnames(P)[max.col(P)]
#' fit_temperature(P, y)
#' @export
fit_temperature <- function(cv_probs, cv_labels, grid = seq(1, 12, 0.5)) {
  cv_probs <- as.matrix(cv_probs)
  sp <- colnames(cv_probs)
  if (is.null(sp)) stop("`cv_probs` must have species column names.")
  cv_acc <- mean(sp[max.col(cv_probs, ties.method = "first")] ==
                   as.character(cv_labels))
  conf <- vapply(grid, function(g)
    mean(apply(temperature_scale(cv_probs, g), 1, max)), numeric(1))
  gamma <- grid[which.min(abs(conf - cv_acc))]
  list(gamma = gamma, cv_accuracy = cv_acc,
       mean_confidence = conf[which.min(abs(conf - cv_acc))])
}
