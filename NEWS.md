# CBFM (development version)

## New features

* `Sigma_control` gains three optional arguments — `prior_precision_space`,
  `prior_precision_time`, and `prior_precision_spacetime` — each a square
  precision matrix that, if supplied, is added directly to the precision of
  the corresponding basis-coefficient random effects. The augmentation feeds
  both the REML lambda-estimation objective in `.update_Sigma_fn` and the
  `invcov` returned by `.update_LoadingSigma_fn`, so the effect propagates
  through downstream covariance / standard error / effective degrees-of-freedom
  calculations without further intervention. Validation (dimension, symmetry,
  positive semi-definiteness) is performed up front; rank-deficient priors are
  accepted. Existing users see no change when the new arguments are omitted.

  The feature has two equivalent readings. Under a penalized-likelihood
  framing, `prior_precision_*` is a user-specified fixed quadratic penalty
  analogous to the `H` argument in `mgcv::gam()`. Under a Bayesian framing,
  it is a zero-mean Gaussian prior on the basis coefficients with the supplied
  precision matrix. A common ridge-style choice is `tau * diag(p)` for a
  positive scalar `tau` and `p = ncol(B_*)`.

## Migration from the internal prototype

Users who had been patching against the earlier two-argument prototype
(`fixed_space` + `fixed_space_tau`, and parallel entries for `time` and
`spacetime`) should switch to the single precision-matrix form:

```r
# Old prototype
Sigma_control <- list(fixed_space = K, fixed_space_tau = tau)

# Equivalent new API
Sigma_control <- list(prior_precision_space = MASS::ginv(K) / tau)
# or, if K is full rank:
Sigma_control <- list(prior_precision_space = solve(K) / tau)
```

The new API takes the precision matrix directly, so users who already have a
precision matrix in hand can supply it as-is rather than inverting to a
covariance first.
