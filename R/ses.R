# Assemblage structure: standardized effect sizes (MPD / MNTD) --------------

# Internal: mean, lo, hi over finite values; degenerate -> 0.
.qci0 <- function(x, ci) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(0, 0, 0))
  c(mean(x), stats::quantile(x, ci[1]), stats::quantile(x, ci[2]))
}
.z0 <- function(x) { x[!is.finite(x)] <- 0; x }

#' Standardized effect sizes of assemblage structure, with uncertainty
#'
#' Standardized effect sizes (SES) of mean pairwise distance (MPD) and mean
#' nearest-taxon distance (MNTD) on the phylogenetic and/or functional axes,
#' computed with \pkg{picante} and propagated through the classification Monte
#' Carlo. Negative SES indicates clustering (co-occurring individuals more
#' related / more similar than chance), positive indicates overdispersion
#' (`NRI = -SES.MPD`, `NTI = -SES.MNTD`). MPD captures whole-assemblage (basal)
#' structure and MNTD terminal (nearest-relative) structure; both stay
#' informative for two-species sites, unlike the functional Hill number at
#' `q = 1` (see [hill_functional()]).
#'
#' Requires the \pkg{picante} package.
#'
#' @param probs Numeric probability matrix (rows = individuals, columns =
#'   species; column names = species labels).
#' @param site Vector of length `nrow(probs)` giving each individual's site.
#' @param tree Optional ape `phylo` (tip labels = species) for the phylogenetic
#'   axis (`pMPD`, `pMNTD`).
#' @param trait_dist Optional functional distance matrix (dimnames = species)
#'   for the functional axis (`fMPD`, `fMNTD`).
#' @param truth Optional vector of known labels for validation (`true_*`
#'   columns).
#' @param n_draws Monte Carlo draws (default 200).
#' @param runs \pkg{picante} null-model permutations per draw (default 499).
#' @param null_model \pkg{picante} null model (default `"taxa.labels"`).
#' @param abundance_weighted Abundance-weight the SES (default `TRUE`).
#' @param ci Length-2 lower/upper interval probabilities (default
#'   `c(0.025, 0.975)`).
#' @param min_k Keep sites with at least this many species (default 2; SES is
#'   undefined below 2). Richness is taken from `truth` when supplied, else from
#'   the most-probable label per individual.
#' @param seed Optional integer seed.
#' @return A data frame, one row per retained site: `site`, `k_species`,
#'   `n_seedlings`, and for each of `pMPD`, `pMNTD`, `fMPD`, `fMNTD` (whichever
#'   axes were requested) a `_ses` mean plus `_lo`/`_hi` interval columns, and
#'   `true_*` columns when `truth` is given.
#' @seealso [propagate_diversity()].
#' @examples
#' \dontrun{
#' ses_structure(P, site, tree = tr, trait_dist = D, truth = y)
#' }
#' @export
ses_structure <- function(probs, site, tree = NULL, trait_dist = NULL,
                          truth = NULL, n_draws = 200, runs = 499,
                          null_model = "taxa.labels",
                          abundance_weighted = TRUE, ci = c(0.025, 0.975),
                          min_k = 2, seed = NULL) {
  if (!requireNamespace("picante", quietly = TRUE))
    stop("ses_structure() needs the 'picante' package. Install it with ",
         "install.packages('picante').")
  if (is.null(tree) && is.null(trait_dist))
    stop("Supply at least one of `tree` or `trait_dist`.")
  probs <- as.matrix(probs)
  species <- colnames(probs)
  if (is.null(species)) stop("`probs` must have species column names.")
  if (!is.null(seed)) set.seed(seed)

  Dphy <- if (!is.null(tree)) stats::cophenetic(tree)[species, species] else NULL
  Dfun <- if (!is.null(trait_dist)) trait_dist[species, species] else NULL
  sites <- sort(unique(site))

  build_comm <- function(lab) {
    m <- matrix(0L, length(sites), length(species),
                dimnames = list(as.character(sites), species))
    for (s in sites)
      m[as.character(s), ] <-
        as.integer(table(factor(lab[site == s], levels = species)))
    m
  }
  one <- function(comm, D)
    picante::ses.mpd(comm, D, null.model = null_model,
                     abundance.weighted = abundance_weighted, runs = runs)$mpd.obs.z
  oneN <- function(comm, D)
    picante::ses.mntd(comm, D, null.model = null_model,
                      abundance.weighted = abundance_weighted, runs = runs)$mntd.obs.z
  metrics <- character(0)
  if (!is.null(Dphy)) metrics <- c(metrics, "pMPD", "pMNTD")
  if (!is.null(Dfun)) metrics <- c(metrics, "fMPD", "fMNTD")

  ses_all <- function(lab) {
    comm <- build_comm(lab); out <- list()
    if (!is.null(Dphy)) { out$pMPD <- .z0(one(comm, Dphy)); out$pMNTD <- .z0(oneN(comm, Dphy)) }
    if (!is.null(Dfun)) { out$fMPD <- .z0(one(comm, Dfun)); out$fMNTD <- .z0(oneN(comm, Dfun)) }
    out
  }

  # Monte Carlo over classification uncertainty
  MC <- lapply(metrics, function(.) matrix(NA_real_, n_draws, length(sites)))
  names(MC) <- metrics
  for (i in seq_len(n_draws)) {
    drawn <- species[apply(probs, 1,
                           function(pr) sample.int(length(species), 1L, prob = pr))]
    s <- ses_all(drawn)
    for (m in metrics) MC[[m]][i, ] <- s[[m]]
  }

  # richness for the min_k filter, and truth SES if requested
  if (!is.null(truth)) {
    kcount <- vapply(sites, function(s)
      length(unique(truth[site == s])), integer(1))
    truthS <- ses_all(as.character(truth))
  } else {
    hard <- species[max.col(probs, ties.method = "first")]
    kcount <- vapply(sites, function(s)
      length(unique(hard[site == s])), integer(1))
    truthS <- NULL
  }
  nseed <- vapply(sites, function(s) sum(site == s), integer(1))

  res <- data.frame(site = sites, k_species = kcount, n_seedlings = nseed,
                    stringsAsFactors = FALSE, check.names = FALSE)
  for (m in metrics) {
    e <- t(apply(MC[[m]], 2, .qci0, ci = ci))
    res[[paste0(m, "_ses")]] <- e[, 1]
    res[[paste0(m, "_lo")]]  <- e[, 2]
    res[[paste0(m, "_hi")]]  <- e[, 3]
    if (!is.null(truthS)) res[[paste0("true_", m)]] <- truthS[[m]]
  }
  res <- res[res$k_species >= min_k, , drop = FALSE]
  rownames(res) <- NULL
  res
}
