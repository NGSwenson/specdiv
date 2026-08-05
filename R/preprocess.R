# Preprocessing for reflectance spectra (spectral front-end) -----------------

#' Standard normal variate (SNV) transform
#'
#' Centre and scale each spectrum (row) by its own mean and standard deviation,
#' removing per-spectrum multiplicative/additive scatter so that classes differ
#' by spectral *shape* rather than overall brightness.
#'
#' @param X Numeric matrix, rows = samples, columns = wavelengths.
#' @return A matrix the same shape as `X`, with wavelength column names kept.
#' @examples
#' X <- matrix(rnorm(20), 4, 5)
#' snv(X)
#' @export
snv <- function(X) {
  X <- as.matrix(X)
  out <- t(scale(t(X)))
  dimnames(out) <- dimnames(X)
  out
}

#' Savitzky-Golay smoothing derivative
#'
#' Fit a local polynomial in a sliding window and report a derivative of it,
#' smoothing noise and removing baseline shifts. Implemented from the normal
#' equations, so there is no package dependency.
#'
#' @param X Numeric matrix, rows = samples, columns = wavelengths (ordered).
#' @param window Odd integer window width (default 11).
#' @param poly Polynomial order (default 2).
#' @param deriv Derivative order to return (default 1; use 0 for smoothing only).
#' @return A matrix the same shape as `X`, with dimnames kept.
#' @examples
#' X <- matrix(cumsum(rnorm(100)), 2, 50)
#' d <- sg_derivative(X)
#' @export
sg_derivative <- function(X, window = 11, poly = 2, deriv = 1) {
  X <- as.matrix(X)
  stopifnot(window %% 2 == 1, poly >= deriv, ncol(X) >= window)
  hw <- (window - 1) / 2
  z  <- -hw:hw
  A  <- outer(z, 0:poly, `^`)                       # Vandermonde design
  co <- (solve(crossprod(A), t(A)) * factorial(deriv))[deriv + 1, ]
  out <- apply(X, 1, function(y) {
    yp <- c(rep(y[1], hw), y, rep(y[length(y)], hw)) # pad ends
    v  <- as.numeric(stats::filter(yp, rev(co), sides = 2))
    v[(hw + 1):(hw + length(y))]
  })
  out <- t(out)
  dimnames(out) <- dimnames(X)
  out
}

#' Preprocess a spectra matrix
#'
#' Apply a named preprocessing pipeline to a spectra matrix. Reference and
#' deployment spectra must be processed identically. The default `"deriv"`
#' (SNV followed by a first derivative) is what the leaf-clip case study uses.
#'
#' @param X Numeric matrix, rows = samples, columns = wavelengths.
#' @param method One of `"deriv"` (SNV + 1st derivative, default), `"snv"`, or
#'   `"raw"` (no transform).
#' @param ... Passed to [sg_derivative()] (e.g. `window`, `poly`).
#' @return A matrix the same shape as `X`, with dimnames kept.
#' @examples
#' X <- matrix(cumsum(rnorm(200)), 4, 50)
#' preprocess_spectra(X, "snv")
#' @export
preprocess_spectra <- function(X, method = c("deriv", "snv", "raw"), ...) {
  method <- match.arg(method)
  X <- as.matrix(X)
  res <- switch(method,
                raw   = X,
                snv   = snv(X),
                deriv = sg_derivative(snv(X), ...))
  dimnames(res) <- dimnames(X)
  res
}

#' Cross-campaign batch alignment
#'
#' Shift a batch of new spectra so that their per-wavelength means match the
#' reference (training) means, removing a whole-batch offset such as a different
#' measurement season or campaign. Unsupervised: uses only the new spectra, no
#' labels. Apply *after* [preprocess_spectra()].
#'
#' Use this when the new data are from a **different** campaign than the library
#' (e.g. lab 2025 -> field 2026). Do **not** use it within a single campaign:
#' there is no offset to remove and it distorts an imbalanced batch, hurting the
#' dominant class. Where possible, re-anchoring the library with a small
#' same-campaign sample is more robust than any global correction.
#'
#' @param new_x Numeric matrix of new/test spectra (rows = samples).
#' @param ref_x Numeric matrix of reference/training spectra (same columns).
#' @return A matrix the same shape as `new_x`.
#' @examples
#' ref <- matrix(rnorm(50), 10, 5)
#' new <- matrix(rnorm(15) + 2, 3, 5)   # offset batch
#' batch_align(new, ref)
#' @export
batch_align <- function(new_x, ref_x) {
  new_x <- as.matrix(new_x); ref_x <- as.matrix(ref_x)
  stopifnot(ncol(new_x) == ncol(ref_x))
  out <- sweep(sweep(new_x, 2, colMeans(new_x)), 2, colMeans(ref_x), "+")
  dimnames(out) <- dimnames(new_x)
  out
}
