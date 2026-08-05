# Spectral front-end: PLS-DA classifier (optional, needs mixOmics) -----------

#' Fit a PLS-DA species classifier on spectra
#'
#' Convenience wrapper around `mixOmics::plsda()` for the leaf-spectra case.
#' Trains a partial least squares discriminant analysis mapping preprocessed
#' spectra to species. Any probabilistic classifier can drive the rest of the
#' framework; this is the one used in the case study.
#'
#' Requires the \pkg{mixOmics} package (Bioconductor).
#'
#' @param x Numeric matrix of preprocessed training spectra (rows = individuals,
#'   columns = wavelengths). See [preprocess_spectra()].
#' @param y Species labels (character/factor), length `nrow(x)`.
#' @param ncomp Number of PLS components (default 24).
#' @return A `mixOmics` `plsda` model object, with attribute `"species"` giving
#'   the class levels in column order.
#' @seealso [predict_probs()], [preprocess_spectra()].
#' @examples
#' \dontrun{
#' Xtr <- preprocess_spectra(train_spectra)
#' m <- fit_spectral_classifier(Xtr, train_species, ncomp = 24)
#' }
#' @export
fit_spectral_classifier <- function(x, y, ncomp = 24) {
  if (!requireNamespace("mixOmics", quietly = TRUE))
    stop("fit_spectral_classifier() needs the 'mixOmics' package ",
         "(Bioconductor): BiocManager::install('mixOmics').")
  y <- factor(y)
  model <- mixOmics::plsda(as.matrix(x), y, ncomp = ncomp)
  attr(model, "species") <- levels(y)
  model
}

#' Predict calibrated-ready species probabilities from a spectral classifier
#'
#' Predict a probability matrix (rows = individuals, columns = species) from a
#' model fitted by [fit_spectral_classifier()]. Raw PLS-DA response scores are
#' turned into probabilities with [scores_to_prob()]; sharpen them with
#' [fit_temperature()] / [temperature_scale()] before propagation.
#'
#' @param model A model from [fit_spectral_classifier()].
#' @param newx Numeric matrix of preprocessed spectra to predict (same
#'   wavelengths as training).
#' @param ncomp Component at which to predict (default: the model's `ncomp`).
#' @return A probability matrix with species column names.
#' @seealso [fit_spectral_classifier()], [scores_to_prob()].
#' @export
predict_probs <- function(model, newx, ncomp = NULL) {
  if (!requireNamespace("mixOmics", quietly = TRUE))
    stop("predict_probs() needs the 'mixOmics' package.")
  species <- attr(model, "species")
  nc <- if (is.null(ncomp)) model$ncomp else ncomp
  Yh <- stats::predict(model, as.matrix(newx))$predict[, , nc]
  Yh <- as.matrix(Yh)
  if (is.null(colnames(Yh))) colnames(Yh) <- species
  scores_to_prob(Yh)
}
