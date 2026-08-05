#' specdiv: uncertainty-aware diversity from predicted species
#'
#' `specdiv` turns per-individual *probabilistic* species predictions into
#' per-site taxonomic, phylogenetic and functional diversity that carries its
#' own uncertainty. The core is a Monte Carlo that draws an identity for each
#' individual from its calibrated probability vector, reconstitutes each site,
#' and recomputes every diversity measure, so each estimate arrives with a
#' confidence interval. All three axes are reported as Hill numbers of matched
#' order q (a common "effective number of species" currency), alongside
#' standardized effect sizes of assemblage structure.
#'
#' The engine ([propagate_diversity()], [ses_structure()]) is agnostic to how
#' the probabilities are produced. A spectral front-end
#' ([fit_spectral_classifier()], [predict_probs()], [preprocess_spectra()]) is
#' provided for the leaf-clip reflectance case the framework was built for, with
#' [batch_align()] for cross-campaign transfer and [apply_abundance_prior()] for
#' the open-species-pool setting.
#'
#' @section Typical workflow:
#' \enumerate{
#'   \item Get a calibrated probability matrix (rows = individuals, columns =
#'     species). Either bring your own, or use [fit_spectral_classifier()] +
#'     [predict_probs()], then [fit_temperature()] + [temperature_scale()].
#'   \item Build a phylogeny (tip labels = species) and a functional distance
#'     matrix ([trait_distance()]).
#'   \item Call [propagate_diversity()] for Hill numbers with intervals, and
#'     [ses_structure()] for MPD/MNTD effect sizes.
#' }
#'
#' @keywords internal
"_PACKAGE"
