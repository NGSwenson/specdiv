# specdiv

<!-- badges: start -->
<!-- badges: end -->

**Uncertainty-aware taxonomic, phylogenetic and functional diversity from
predicted species.**

`specdiv` turns per-individual *probabilistic* species predictions into
per-site diversity that carries its own uncertainty. A calibrated classifier
assigns each individual a probability vector over candidate species; a Monte
Carlo repeatedly draws identities, reconstitutes each site, and recomputes
diversity, so every estimate arrives with a confidence interval. All three axes
are returned as **Hill numbers** of matched order — a common "effective number
of species" currency — alongside standardized effect sizes (MPD/MNTD) of
assemblage structure.

The engine is agnostic to how the probabilities arise. It was developed for
nondestructive leaf-clip reflectance spectra of forest seedlings, and a spectral
front-end (PLS-DA via **mixOmics**, preprocessing, cross-campaign batch
alignment, and an independent abundance prior) is included for that case.

## Installation

```r
# install.packages("remotes")
remotes::install_github("NGSwenson/specdiv")
```

The core engine depends only on `ape`. Two optional features pull heavier
packages, loaded only when used:

- `ses_structure()` needs [`picante`](https://cran.r-project.org/package=picante)
  (`install.packages("picante")`).
- the spectral classifier (`fit_spectral_classifier()`, `predict_probs()`)
  needs [`mixOmics`](https://bioconductor.org/packages/mixOmics)
  (`BiocManager::install("mixOmics")`).

## Quick start

```r
library(specdiv)

# P: individuals x species probability matrix (colnames = species)
# tree: ape phylo with those species as tips
# Dfun: functional distance matrix (trait_distance(traits))

div <- propagate_diversity(
  P, site,
  tree = tree, trait_dist = Dfun, truth = truth,
  q = c(0, 1, 2), n_draws = 500
)
# -> per-site TD/PD/FD Hill numbers with _lo/_hi intervals (+ true_* if truth given)

ses <- ses_structure(P, site, tree = tree, trait_dist = Dfun, truth = truth)
# -> per-site pMPD/pMNTD/fMPD/fMNTD effect sizes with intervals
```

See `vignette("specdiv")` for the full walkthrough, including calibration and
the spectral front-end.

## What each piece does

| Function | Role |
|---|---|
| `propagate_diversity()` | Monte Carlo Hill numbers (3 axes, orders q) with CIs |
| `ses_structure()` | MPD/MNTD standardized effect sizes with CIs (needs `picante`) |
| `diversity_from_labels()` | diversity from known identities (the validation target) |
| `hill_taxonomic()`, `hill_phylo()`, `hill_functional()` | the underlying Hill-number formulas |
| `fit_temperature()`, `temperature_scale()`, `scores_to_prob()` | probability calibration |
| `trait_distance()` | functional distances from traits |
| `apply_abundance_prior()` | Bayesian prior for the open-species-pool setting |
| `preprocess_spectra()`, `snv()`, `sg_derivative()` | spectral preprocessing |
| `batch_align()` | cross-campaign batch correction |
| `fit_spectral_classifier()`, `predict_probs()` | PLS-DA front-end (needs `mixOmics`) |

## Citation

If you use `specdiv`, please cite the accompanying paper (Swenson et al., in
prep) and this package.

## License

MIT © Nathan G. Swenson
