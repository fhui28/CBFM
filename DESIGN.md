# DESIGN: Prior precision matrix for basis-coefficient random effects

**Branch:** `feature/precision-prior` (forked from `allowingweights`)
**Author:** Christopher Haak (chrishaak)
**Status:** Design — implementation complete, pending review.

---

## 1. Feature statement

Extend `Sigma_control` with three optional entries:

- `prior_precision_space`
- `prior_precision_time`
- `prior_precision_spacetime`

Each, if supplied, is a square symmetric positive semi-definite matrix of
dimension matching the corresponding basis (`ncol(B_space)`, etc.). CBFM
adds it directly to the precision of the basis-coefficient random effects,
bypassing the lambda-estimation and G-estimation machinery entirely.

The motivating use case is resolving null-space identifiability problems in
CBFM's two-block estimator. The feature is mathematically general and admits
multiple equivalent interpretations (see Section 2.4).

## 2. Motivation

### 2.1 The null-space problem in two-block estimation

CBFM estimates species-specific spatial, temporal, and spatio-temporal
effects through basis function regression, with a community-level covariance
structure `Sigma` (and optional between-species correlation `G`) governing
how coefficients are shared across species. When a user supplies a
structural penalty through `Sigma_control$custom_*` — e.g., from `mgcv` — that
penalty defines which directions of the coefficient space are constrained
by smoothness and which are left free.

Most structural penalties are **rank-deficient**: they do not constrain every
direction. The unconstrained directions form the penalty's **null space**.
The most common instance is the constant direction: a uniform shift across
all basis functions that produces a flat component of the fitted surface.

The intercept in the GAM (`X`) block is already present to capture each
species' baseline abundance. A constant component of the basis block does
the same job. The two components are mathematically equivalent for
prediction purposes — but they are not jointly identified. During PQL
estimation, the two blocks can oscillate between equivalent representations,
inflating standard errors and degrading convergence.

This problem is not specific to spatial smooths. It appears across many
basis types used in CBFM:

- Spatial smooths (Duchon, TPS, GP): null space = constant direction.
- Random walk on time (RW1): null space = constant direction.
- Tensor product smooths (`t2`): null space = intercept-like direction
  arising from the tensor construction.
- Nested factor random effects: null space = constant direction, orthogonal
  to both between-group and within-group contrasts.

### 2.2 Why the null-space penalty must be separated from the estimated system

A natural first attempt is to treat the null-space penalty as just another
structural penalty: add S_null alongside the real penalties as one element
of `Sigma_control$custom_*`, let REML estimate its lambda (possibly with an
upper bound to prevent it from collapsing to zero penalty).

This works when `G_control$structure = "identity"` — each penalty gets its
own scalar variance, so REML can drive S_null's lambda against a cap
without interfering with the structural penalties' scales.

It fails when `G_control$structure` is `"homogeneous"` or `"unstructured"`.
Those settings make G estimation attempt to fit a shared correlation
structure across **all** penalties in `custom_*`, including S_null. But
S_null's role is identifiability, not smoothness. It is conceptually
different from the structural penalties: it represents a fixed prior belief
that a specific direction should be penalized, not a data-driven shrinkage
target whose variance should be estimated jointly with the others.
Including it in G estimation produces systematic failures — the optimizer
tries to model correlation on a dimension that shouldn't be modeled at all.

A related failure mode: even where REML technically works (identity G),
the optimizer spends effort pushing S_null's lambda toward its cap on
every iteration. The penalty's "correct" strength is not a quantity to
estimate; it is a user-specified constant. Routing it through REML wastes
optimization and creates spurious coupling between its scale and the
structural scales.

The fixed-precision approach resolves both issues by taking the
null-space penalty out of the estimated system entirely. It is added to
the basis-coefficient precision at a fixed scale the user chooses, where
it influences estimation but is invisible to both the REML lambda
optimization (see Section 5.2 on why lambda estimation remains consistent)
and the G optimization.

### 2.3 The interface

    Sigma_control <- list(
        custom_space = list(S1, S2),               # structural, lambda-estimated
        prior_precision_space = prior_precision    # fixed, bypasses lambda and G
    )

The structural penalties in `custom_*` are estimated as before. The prior
precision in `prior_precision_*` is a fixed quantity controlled entirely
by the user.

### 2.4 Interpretations

The feature admits three equivalent mathematical framings:

**Null-space identifiability (primary use case).** As described in 2.1–2.3.
The prior is chosen to target directions the structural penalties leave
unconstrained, resolving block-level ambiguity.

**mgcv `H`-matrix analogue.** mgcv's `gam()` accepts an `H` argument — a
user-specified supplementary fixed penalty added to the data-driven
smoothing penalty. `prior_precision_*` is the CBFM equivalent for
basis-coefficient random effects.

**Bayesian prior on basis coefficients.** If the data-driven precision on
the basis coefficients is `P_data`, supplying `prior_precision_* = M`
produces a posterior precision `P_data + M`. Equivalently, this is a
zero-mean Gaussian prior on the coefficients with precision `M`.

Users are free to choose whichever framing they find most natural. The
implementation is identical across framings.

### 2.5 Orthogonality: why this does not distort structural estimation

Null-space penalties constructed for the use case in 2.1 are orthogonal to
the structural penalties by construction. For an RW1, S_null projects
onto the constant vector, which is in the null space of the difference
operator S_rw1. For a tensor product smooth, S_null derived from
column means targets the intercept-like direction the structural penalties
do not constrain. For a nested factor RE, S_null projects onto the
constant, orthogonal to the group contrasts.

Because the prior adds precision only on directions where the structural
penalties contribute zero, the REML gradient with respect to the lambdas
is unchanged — lambda estimates are, in the idealized limit, unaffected
by the presence of a well-constructed null-space prior. G estimation is
likewise unaffected: G operates on the structurally-penalized dimensions
and the prior lies outside its scope.

This property is the reason the fixed-precision approach works: it provides
identifiability without contaminating the estimation of the parts of the
model the user actually wants data-driven.

## 3. API

### 3.1 Chosen API

    Sigma_control <- list(
        prior_precision_space = M_space,         # p_s x p_s matrix, optional
        prior_precision_time = M_time,           # p_t x p_t matrix, optional
        prior_precision_spacetime = M_spacetime  # p_st x p_st matrix, optional
    )

Each argument, if supplied, is a **precision matrix supplied directly** —
no inversion, no transformation, no separate scale parameter. Users who
conceptually separate the shape of the penalty from its strength can
express this at the construction site:

    prior_precision_space <- strength * base_penalty_matrix

The API is deliberately thin: CBFM accepts a precision matrix; how the
matrix is constructed is the user's responsibility. For the null-space
use case, construction helpers live naturally in downstream tooling
(e.g., alongside `construct_smooth_basis`, `construct_random_walk`) rather
than inside CBFM itself.

### 3.2 Rejected alternatives

**Two-argument API `(fixed_*, fixed_*_tau)`.** Original prototype form.
Rejected on API review with Francis Hui: it imposes a specific
decomposition `precision = pinv(K) / tau` that users may not want;
requires six arguments rather than three; and forces a pointless
round-trip (invert, pass, un-invert) for users who already have a
precision matrix in hand.

**Naming as `H_*` to mimic mgcv.** Considered. Rejected in favor of
`prior_precision_*` for self-documentation; `H_*` is only transparent
to readers familiar with mgcv internals.

**Placement as top-level arguments rather than inside `Sigma_control`.**
Considered. Rejected because the feature concerns the Sigma structure
and naturally belongs there, alongside `custom_*`.

### 3.3 Migration from the internal prototype

Users who had been patching against the original two-argument prototype
can migrate:

    # Old prototype
    Sigma_control <- list(fixed_space = K, fixed_space_tau = tau)

    # Equivalent new API
    Sigma_control <- list(prior_precision_space = MASS::ginv(K) / tau)
    # or, if K is full-rank:
    Sigma_control <- list(prior_precision_space = solve(K) / tau)

This migration note is captured in `NEWS.md`.

## 4. Validation

Implemented in `R/checkfillfns.R`. For each of `prior_precision_{space,
time, spacetime}`, if supplied:

1. Corresponding `B_` basis must be in use (error if not).
2. Must be a single matrix (not a list; errors otherwise).
3. Must be square with dimensions matching `ncol(B_{space,time,spacetime})`.
4. Must be symmetric (up to `sqrt(.Machine$double.eps)`).
5. Must be positive semi-definite — eigenvalues >= `-sqrt(.Machine$double.eps)`.

**Rank-deficient priors are explicitly accepted.** The canonical
null-space use case produces rank-1 matrices (e.g., `tcrossprod(col_means)`);
rejecting rank-deficient inputs would render the feature useless for its
primary purpose. Validation uses `>= 0` eigenvalue tolerance, and the
implementation handles pseudo-inverse operations via `.pinv()` wherever
`solve()` would fail on a singular input.

## 5. Implementation

The feature touches three files at three code sites. The architecture
mirrors CBFM's existing data flow: the prior must be visible to the REML
objective during lambda estimation, present in the returned `invcov` from
`.update_LoadingSigma_fn` so downstream computations inherit it, and
applied at initialization for the control path where the outer Sigma
update is skipped.

### 5.1 `R/CBFM.R` — initialization-time augmentation

For the code path where `Sigma_control$custom_*` is supplied as a single
matrix under an unstructured G (case 3 in §5.4), the outer Sigma update
is skipped entirely. The `invcov` assembled at initialization from
`.pinv(custom_*)` is the final value used throughout estimation. Without
an initialization-time hook, the prior would never apply for this case.

Three blocks in `CBFM.R` (one each for space, time, spacetime), following
the initial `.pinv(custom_*)` assignment, augment `invcov`:

    if (!is.null(Sigma_control[["prior_precision_space"]])) {
        new_LoadingnuggetSigma_space$invcov <-
            new_LoadingnuggetSigma_space$invcov +
            Sigma_control[["prior_precision_space"]]
        new_LoadingnuggetSigma_space$cov <-
            .pinv(new_LoadingnuggetSigma_space$invcov)
    }

(Analogous blocks for time and spacetime.)

### 5.2 `R/updateGSigma.R` — `fn_lambda` augmentation inside `.update_Sigma_fn`

In both the "single matrix" and "list of matrices" code paths of the
REML objective function `fn_lambda`, the prior is included in the
objective's precision term:

    Sigmainv_total <- expx * Sigmainv_space + prior_precision_space

Consequences:

- `trace_quantity` is computed against the augmented precision.
- Log-determinant uses `Sigmainv_total` via an eigenvalue-based
  computation (`sum(log(e[e > .Machine$double.eps]))`), replacing the
  original `num_basisfns * log(expx)` closed-form. The eigenvalue
  approach is robust to rank-deficient priors.
- The `kronecker(Ginv, ...)` term uses `Sigmainv_total`.

**On the orthogonality argument.** Section 2.5 noted that a well-constructed
null-space prior leaves the REML gradient with respect to lambdas unchanged,
and thus would not in theory affect lambda estimates. The prior is
nevertheless included in `fn_lambda` for two reasons:

1. **Numerical robustness.** Orthogonality is only exact in idealized
   settings; user-supplied priors may deviate slightly. Including the
   prior in the objective ensures consistent accounting regardless of
   exact orthogonality.
2. **Consistency.** The prior is part of the precision; the REML
   objective evaluates the likelihood implied by the precision. Omitting
   the prior would evaluate an objective inconsistent with the prior
   that actually applies to the coefficients.

### 5.3 `R/updateGSigma.R` — `.add_prior_precision` helper in `.update_LoadingSigma_fn`

For the main estimation paths (cases 1, 2, and 4 in §5.4), the data-driven
`invcov` assembled inside `.update_LoadingSigma_fn` is augmented with the
prior via a helper:

    .add_prior_precision <- function(out, prior_precision) {
        if (!is.null(prior_precision)) {
            out$invcov <- out$invcov + prior_precision
            out$cov <- .pinv(out$invcov)
        }
        out
    }

Called in both branches (`estimate_lambda_not_Sigma = 0` and `= 1`), after
`invcov` is assembled from the structural penalties and estimated lambdas
but before the function returns. This ensures the augmented `invcov`
propagates to every downstream consumer — standard error computation,
effective degrees of freedom, `make_tidibits_data()` — without further
intervention elsewhere.

### 5.4 Relationship between the three sites across estimation cases

CBFM's estimator follows different code paths depending on whether the
user supplies a `custom_*` Sigma and whether the block's G structure is
identity/homogeneous or unstructured. This gives four cases:

| Case | `custom_*`? | G structure | Outer Sigma update | Augmentation applied via |
|---|---|---|---|---|
| 1 | no | unstructured | `.update_LoadingSigma_fn(estimate_lambda_not_Sigma = 0)` | §5.3 (branch 1) |
| 2 | no | identity/homogeneous | `.update_LoadingSigma_fn(estimate_lambda_not_Sigma = 0)` | §5.3 (branch 1) |
| 3 | yes, as matrix | unstructured | **skipped** — `invcov` is fixed at initialization | **§5.1** (init-time) |
| 4 | yes | identity/homogeneous | `.update_LoadingSigma_fn(estimate_lambda_not_Sigma = 1)` | §5.3 (branch 2) |

§5.2 (`fn_lambda` augmentation) applies whenever the outer Sigma update
is called — i.e., cases 1, 2, and 4.

Cases 1 and 2 share the same augmentation path; the G-structure difference
affects other parts of the estimator but not where the prior enters.
Case 3 is the only case where the outer update is skipped entirely,
which is why the §5.1 initialization-time hook exists.

Each augmentation site fires exactly once per outer iteration for its
applicable case. No double-counting; no case left uncovered.

## 6. EDF bookkeeping

CBFM tracks effective degrees of freedom and standard errors via the
stored `invcov` returned by `.update_LoadingSigma_fn` and subsequently
read by `make_tidibits_data()` and the various downstream summary
functions. Because augmentation happens inside `.update_LoadingSigma_fn`
(§5.3) and at initialization for the case-3 path (§5.1), the returned
`invcov` already incorporates the prior. Downstream EDF, SE, and variance
partitioning calculations consume this augmented `invcov` and thus
reflect the prior's contribution automatically — no separate EDF-related
modifications were required.

Using the augmented precision for uncertainty quantification is the
statistically coherent choice under both the mgcv-H and Bayesian framings:
the prior contributes information that sharpens the effective parameter
count and narrows the implied uncertainty. For the null-space use case
specifically, this is analogous to how mgcv's post-constraint SEs are
reported (the sum-to-zero constraint informs the SEs rather than being
ignored).

## 7. Documentation

Roxygen additions in the commits below document the three new
`Sigma_control` entries: their type, dimensions, validation requirements,
both the mgcv-`H` and Bayesian framings, the null-space identifiability
use case, and a minimal construction example.

The `NEWS.md` entry describes the feature, the equivalence of framings,
and the migration from the internal prototype. A version heading will be
set by the maintainer at release time.

## 8. Testing

`tests/testthat/test-prior-precision.R` exercises:

- **No-op when NULL.** Fit identical with `prior_precision_*` unset versus
  explicitly `NULL`.
- **Strong ridge prior shrinks loadings.** A large `strength * diag(p)`
  prior drives loadings toward zero as expected.
- **Input validation.** Non-matrix input, wrong-dimension matrix,
  non-symmetric matrix, non-PSD matrix, and supplying the argument
  without a corresponding `B_` basis all raise informative errors.

The suite (15 tests total) runs in under 30 seconds on a small simulated
dataset.

## 9. Commit sequence (as landed)

    1c067be  docs: add DESIGN.md for prior precision feature
    7955375  feat: add input validation for prior_precision_* in checkfillfns.R
    7709b79  feat: augment LoadingnuggetSigma invcov with prior_precision in CBFM.R
    6192c2c  feat: augment Sigmainv_total in fn_lambda with prior_precision
    2966d42  feat: augment invcov in .update_LoadingSigma_fn via .add_prior_precision
    54bc496  feat: relocate prior_precision augmentation sites in CBFM.R
    0196386  test: add prior_precision feature tests
    85db3d5  docs: roxygen documentation for prior_precision_*
    278a48d  docs: NEWS.md entry for prior_precision feature
    85f51db  chore: add .claude to .gitignore
    955687a  build: add DESIGN.md and .claude to .Rbuildignore

Commits `7709b79` and `54bc496` form a refactor: the former initially
placed augmentation at the outer-loop sites; the latter relocated it to
the correct initialization-time locations once the case analysis was
worked out (§5.4), to avoid double-counting in case 3.

## 10. Open API questions for review

- Naming: `prior_precision_*` reflects Francis's feedback and the project's
  naming preference; `H_*` remains an alternative if closer alignment with
  mgcv is preferred.
- Placement: within `Sigma_control` rather than top-level. The feature
  concerns the Sigma structure and naturally belongs there alongside
  `custom_*`.
