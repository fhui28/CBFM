# DESIGN: Prior precision matrix for basis-coefficient random effects

**Branch:** `feature/precision-prior` (forked from `allowingweights`)
**Author:** Christopher Haak (chrishaak)
**Status:** Design — not yet implemented

---

## 1. Feature statement

Allow users of CBFM to supply a fixed, user-specified precision matrix that is
added to the precision of the basis-coefficient random effects during
estimation. The fixed precision is *not* subject to lambda / G estimation —
it enters the model as supplied.

This extends `Sigma_control` with three optional new entries:

- `prior_precision_space`
- `prior_precision_time`
- `prior_precision_spacetime`

Each, if supplied, is a square matrix of the same dimension as the
corresponding basis (`ncol(B_space)`, etc.).

## 2. Motivation and framing

The feature has two natural interpretations, both valid:

**Frequentist / smoothing-penalty framing.** Analogous to mgcv's `H`
argument — a user-specified supplementary fixed penalty added to the
data-driven smoothing penalty. Unlike mgcv's `S` matrices (which have
associated smoothing parameters estimated from the data), `H` enters at a
fixed scale chosen by the user.

**Bayesian framing.** A prior on the basis-coefficient random effects.
If the existing loading-nugget precision is `P_data` and the user supplies
`prior_precision_space = M`, the joint prior on the loadings becomes
multivariate Gaussian with precision `P_data + M`. Users who think in
Bayesian terms can specify `M` directly as the prior precision.

Both framings produce identical behavior. Documentation should
acknowledge both.

## 3. API

### 3.1 Chosen API (after review with Francis Hui)

```r
Sigma_control <- list(
    prior_precision_space = M_space,         # p_s x p_s matrix, optional
    prior_precision_time = M_time,           # p_t x p_t matrix, optional
    prior_precision_spacetime = M_spacetime  # p_st x p_st matrix, optional
)
```

Each argument, if supplied, is a **precision matrix supplied directly** —
no transformation, no separate scale parameter. The user is responsible for
constructing whatever precision matrix they want.

### 3.2 Rejected alternative APIs and why

**Two-argument API (original prototype):**

```r
Sigma_control <- list(
    fixed_space = K_space,         # covariance-shaped
    fixed_space_tau = tau_space    # scale
)
# Internally: precision = .pinv(K_space) / tau_space
```

Rejected because:

- Imposes a specific decomposition (`pinv(K) / tau`) that users may not want.
- Adds complexity (two arguments per block, six total).
- Users wanting a generic precision matrix have to invert their own precision,
  pass it as `K`, and set `tau = 1` — a pointless round trip.
- "Fixed" is a loaded term that overlaps confusingly with "fixed effect".

**Alternative naming `H_{space,time,spacetime}`** (mgcv mimicry):
Considered. Rejected in favor of `prior_precision_*` for clarity; `H_*` is
only meaningful to users who know mgcv internals. Open to revisiting if
Francis prefers.

### 3.3 Migration from the prototype

Users who already have the old-API inputs `K` (covariance-shaped) and
`tau` (scale):

```r
# Old prototype
Sigma_control <- list(fixed_space = K, fixed_space_tau = tau)

# Equivalent new API
Sigma_control <- list(prior_precision_space = MASS::ginv(K) / tau)
# or, if K is full-rank:
Sigma_control <- list(prior_precision_space = solve(K) / tau)
```

This migration snippet goes in `NEWS.md`.

## 4. Validation

Implemented in `R/checkfillfns.R`. For each of `prior_precision_{space,
time, spacetime}`, if supplied:

1. Corresponding `B_` basis must be in use (error if not).
2. Must be a single matrix (not a list; errors otherwise).
3. Must be square with dimensions matching `ncol(B_{space,time,spacetime})`.
4. Must be symmetric (up to `sqrt(.Machine$double.eps)`).
5. Must be positive semi-definite (eigenvalues >= `-sqrt(.Machine$double.eps)`).

Symmetry and PSD checks are additions over the prototype — the prototype
had dimension checks only. Both checks are cheap relative to model fitting
and catch common user errors.

## 5. Implementation

The feature touches three files. At a high level, the precision
contribution is added at two distinct points in the estimator:

- **Inside the REML objective `fn_lambda`** in `updateGSigma.R`, so that
  the data-driven `lambda` parameters are estimated *accounting for* the
  prior.
- **After the data-driven update** of the loading-nugget `invcov` in both
  `updateGSigma.R` (`.update_LoadingSigma_fn`) and `CBFM.R` (the outer
  iteration loop), so that downstream operations see the augmented precision.

### 5.1 `R/CBFM.R` — post-update augmentation

After each of `new_LoadingnuggetSigma_{space,time,spacetime}` is updated
in the outer iteration, augment the precision with the fixed prior:

```r
if (!is.null(Sigma_control[["prior_precision_space"]])) {
    new_LoadingnuggetSigma_space$invcov <-
        new_LoadingnuggetSigma_space$invcov +
        Sigma_control[["prior_precision_space"]]
    new_LoadingnuggetSigma_space$cov <- .pinv(new_LoadingnuggetSigma_space$invcov)
}
```

(Parallel blocks for time and spacetime.)

### 5.2 `R/updateGSigma.R` — augmentation inside REML objective

In both the "single matrix" and "list of matrices" code paths in
`.update_Sigma_fn`, the `fn_lambda` closure is modified to use a
prior-augmented precision matrix `Sigmainv_total`:

- Pre-compute `prior_precision_*` once outside the closure (it's constant
  across optim iterations).
- Inside the closure:
  - `Sigmainv_total <- data_driven_precision + prior_precision_*`
  - Log-determinant computed from eigenvalues of `Sigmainv_total` with
    a floor at `.Machine$double.eps`, replacing the original closed-form
    `num_basisfns * log(expx)`.
  - `trace_quantity` and `kronecker(Ginv, ...)` use `Sigmainv_total`.

**Conditional fast path.** When no prior is supplied, the code should
take the original (cheaper) code path using `num_basisfns * log(expx)`
directly. The eigenvalue-based log-det is only used when the prior is
non-NULL. This ensures zero performance regression for existing users.

```r
if (is.null(prior_precision_space)) {
    # Original fast path — unchanged
    out <- 0.5*num_spp*num_basisfns*log(expx) - 0.5*expx*trace_quantity
    out <- out - 0.5*determinant(BtKB + expx*KronGinvSigmainv)$mod
} else {
    # Augmented path
    Sigmainv_total <- expx * Sigmainv_space + prior_precision_space
    e <- eigen(Sigmainv_total, only.values = TRUE, symmetric = TRUE)$values
    trace_quantity <- sum(diag(Sigmainv_total %*% AT_Ginv_A))
    out <- 0.5*num_spp*sum(log(e[e > .Machine$double.eps])) - 0.5*trace_quantity
    out <- out - 0.5*determinant(BtKB + kronecker(Ginv, Sigmainv_total))$mod
}
```

### 5.3 `R/updateGSigma.R` — `.update_LoadingSigma_fn` augmentation

After the data-driven `invcov` is assembled from `custom_*` and the
estimated `lambdas`, augment with the prior and recompute `cov`:

```r
if (!is.null(Sigma_control[["prior_precision_space"]])) {
    out$invcov <- out$invcov + Sigma_control[["prior_precision_space"]]
    out$cov <- .pinv(out$invcov)
}
```

Parallel for time and spacetime. In the prototype this was factored into
a `.add_fixed_penalty` helper; worth keeping that factoring since it's
called at three parallel sites. Rename to `.add_prior_precision` for
consistency with the new naming.

### 5.4 Factoring

The prototype has ~6 near-identical code blocks across the three
basis types in `fn_lambda`, plus three in `.update_LoadingSigma_fn`.
Some factoring is worth doing (especially `.update_LoadingSigma_fn`).
`fn_lambda` factoring is more subtle because of the closure variables;
judge by whether the resulting helper is clearer than the duplication.

## 6. Open technical questions

### 6.1 EDF bookkeeping

CBFM tracks effective degrees of freedom (EDF) via precision matrices.
The augmented `invcov` affects any downstream EDF calculation. Need to:

1. Grep the codebase for EDF / edf / effective.df references.
2. Determine which precision matrix (augmented or un-augmented) each
   reference should use.
3. Verify the prototype's behavior; fix if incorrect.

This is the single most important technical question to resolve during
implementation. It was not addressed in the prototype.

### 6.2 Performance for large bases

`eigen()` in the augmented `fn_lambda` path is O(p^3) per optim iteration.
For `p = 30` (surfclam spatial basis) this is trivial. For larger bases
it could matter. Cholesky-based log-determinant would be faster but less
robust to near-singular matrices.

Resolution: accept `eigen()` for now. Document performance consideration
in the roxygen. Optimize only if a user reports a bottleneck.

### 6.3 Interaction with `custom_*` penalties

Users can already supply `Sigma_control$custom_space` (a single matrix or
list). The prior precision is *additive* to whatever the custom / estimated
structure produces. Need to verify this is the intended semantics and that
edge cases (e.g., `custom_space` matrix that already has `prior_precision`
baked in) don't produce double-counting. Most likely fine as long as
documentation is clear.

## 7. Documentation

Roxygen additions to whatever function(s) document `Sigma_control`:

- Describe each of `prior_precision_{space,time,spacetime}` as a fixed
  precision matrix added to the basis-coefficient random-effects precision.
- Include both framings (supplementary penalty; Bayesian prior).
- Document dimension, symmetry, and PSD requirements.
- Cross-reference mgcv's `H` argument for readers familiar with it.
- Include a brief example construction (e.g., `tau * diag(p)` as a
  ridge-style prior).

`NEWS.md`:

- Entry describing the new feature.
- Include the migration snippet from section 3.3 for anyone who had been
  patching against a prototype API.

## 8. Testing

`tests/testthat/test-prior-precision.R` (new file):

**Test 1 — no-op when NULL.** Fit a small CBFM model with and without
`Sigma_control$prior_precision_space = NULL`. Results must be bit-identical
(confirms the NULL branch is truly a no-op).

**Test 2 — strong prior shrinks loadings.** Fit with
`prior_precision_space = 1e6 * diag(p)` (very strong ridge-style prior on
space loadings). Verify loadings are near zero vs. the unconstrained fit.

**Test 3 — input validation.** Confirm errors are raised for:

- Non-matrix input
- Wrong-dimension matrix
- Non-symmetric matrix
- Non-PSD matrix
- Supplied without corresponding `B_` basis

**Test 4 (optional) — equivalence to custom penalty.** Compare a fit
using `prior_precision_space = M` against a fit using an equivalent
structure through `Sigma_control$custom_space`, where equivalence can be
constructed. This checks that the two mechanisms compose sensibly.

Tests should use a tiny simulated dataset so they run in seconds.

## 9. Implementation plan (commit sequence)

Proposed commits on `feature/precision-prior`:

1. **`docs: add DESIGN.md for prior precision feature`** — this document.
2. **`feat: add input validation for prior_precision_* in checkfillfns.R`** —
   validation block only, no functional change.
3. **`feat: augment LoadingnuggetSigma invcov with prior_precision in CBFM.R`** —
   the three 6-line blocks in the outer loop.
4. **`feat: augment Sigmainv_total in fn_lambda with prior_precision`** —
   the REML objective modifications, with the conditional fast path.
5. **`feat: augment invcov in .update_LoadingSigma_fn with .add_prior_precision helper`** —
   post-update augmentation + the helper.
6. **`fix: adjust EDF bookkeeping for prior_precision`** (if needed, based
   on investigation in 6.1).
7. **`docs: roxygen documentation for prior_precision_*`**.
8. **`test: add prior_precision feature tests`**.
9. **`docs: NEWS.md entry for prior_precision feature`**.

Each commit should leave the package in a working state (tests pass, or
at least don't introduce new failures).

## 10. Things deferred to other PRs

Discovered in the prototype but *not* part of this PR:

- **Adaptive G damping** (`G_dampen`, `adaptive_G_dampen`): a separate
  numerical-stability feature touching the same files. Will be a future
  PR if kept at all. Not merged into this branch.
- **Removal of debug `cat()` statements**: the prototype has ~20 such
  statements. None are carried forward into this branch.

## 11. Open API questions for Francis (to resolve before merging)

- Naming: `prior_precision_*` vs `H_*` — defaulting to the former;
  happy to rename if preferred.
- Placement within `Sigma_control` vs a top-level argument — the
  prototype uses `Sigma_control`; this seems right (the feature
  concerns the Sigma structure) but open to alternatives.
