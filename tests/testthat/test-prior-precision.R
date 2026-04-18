# Tests for the prior_precision_{space,time,spacetime} feature
# (Sigma_control$prior_precision_*). See DESIGN.md for the design.

# ---------------------------------------------------------------------------
# Validation — cheap, no model fit required
# ---------------------------------------------------------------------------

test_that(".check_prior_precision accepts NULL and valid PSD matrices", {
    expect_silent(
        CBFM:::.check_prior_precision(
            NULL, in_use = TRUE, dim_B = 5,
            which_name = "prior_precision_space", basis_name = "B_space"
        )
    )
    # rank-deficient PSD (eigenvalues 2 and 0)
    M <- matrix(c(1, 1, 1, 1), 2, 2)
    expect_silent(
        CBFM:::.check_prior_precision(
            M, in_use = TRUE, dim_B = 2,
            which_name = "prior_precision_space", basis_name = "B_space"
        )
    )
    # full-rank PSD
    expect_silent(
        CBFM:::.check_prior_precision(
            diag(5), in_use = TRUE, dim_B = 5,
            which_name = "prior_precision_space", basis_name = "B_space"
        )
    )
})

test_that(".check_prior_precision rejects all documented bad inputs", {
    # prior supplied without corresponding basis
    expect_error(
        CBFM:::.check_prior_precision(
            diag(5), in_use = FALSE, dim_B = 5,
            which_name = "prior_precision_space", basis_name = "B_space"
        ),
        "if B_space is also not supplied"
    )
    # list instead of a matrix
    expect_error(
        CBFM:::.check_prior_precision(
            list(diag(5)), in_use = TRUE, dim_B = 5,
            which_name = "prior_precision_space", basis_name = "B_space"
        ),
        "single matrix, not a list"
    )
    # non-matrix (vector)
    expect_error(
        CBFM:::.check_prior_precision(
            1:25, in_use = TRUE, dim_B = 5,
            which_name = "prior_precision_space", basis_name = "B_space"
        ),
        "should be a matrix"
    )
    # wrong dimension
    expect_error(
        CBFM:::.check_prior_precision(
            diag(4), in_use = TRUE, dim_B = 5,
            which_name = "prior_precision_space", basis_name = "B_space"
        ),
        "same dimensions as ncol"
    )
    # non-symmetric
    expect_error(
        CBFM:::.check_prior_precision(
            matrix(c(1, 0.5, 0, 1), 2, 2), in_use = TRUE, dim_B = 2,
            which_name = "prior_precision_space", basis_name = "B_space"
        ),
        "symmetric"
    )
    # non-PSD (eigenvalues 3, -1)
    expect_error(
        CBFM:::.check_prior_precision(
            matrix(c(1, 2, 2, 1), 2, 2), in_use = TRUE, dim_B = 2,
            which_name = "prior_precision_space", basis_name = "B_space"
        ),
        "positive semi-definite"
    )
})


# ---------------------------------------------------------------------------
# End-to-end fits — skipped on CRAN and when autoFRK is unavailable
# ---------------------------------------------------------------------------

.make_tiny_cbfm_dataset <- function() {
    set.seed(1)
    N <- 60
    S <- 6
    xy <- data.frame(x = stats::runif(N, 0, 5), y = stats::runif(N, 0, 5))
    X <- matrix(stats::rnorm(N * 2), N)
    colnames(X) <- c("a", "b")
    dat <- cbind(xy, as.data.frame(X))
    Bsp <- as.matrix(autoFRK::mrts(dat[, c("x", "y")], 6))[, -1, drop = FALSE]
    y <- matrix(
        stats::rbinom(N * S, 1, stats::plogis(stats::rnorm(N * S, sd = 0.5))),
        N, S
    )
    colnames(y) <- paste0("spp", seq_len(S))
    list(y = y, data = dat, Bsp = Bsp)
}

test_that("prior_precision_space = NULL is identical to the argument being omitted", {
    skip_on_cran()
    skip_if_not_installed("autoFRK")
    d <- .make_tiny_cbfm_dataset()

    fit_missing <- CBFM(
        y = d$y, formula = ~ a + b, data = d$data, B_space = d$Bsp,
        family = stats::binomial(),
        control = list(trace = 0, maxit = 5),
        stderrors = FALSE, ncores = 1,
        Sigma_control = list(rank = 2),
        G_control = list(rank = 2)
    )
    fit_null <- CBFM(
        y = d$y, formula = ~ a + b, data = d$data, B_space = d$Bsp,
        family = stats::binomial(),
        control = list(trace = 0, maxit = 5),
        stderrors = FALSE, ncores = 1,
        Sigma_control = list(rank = 2, prior_precision_space = NULL),
        G_control = list(rank = 2)
    )

    expect_identical(fit_missing$basis_effects_mat, fit_null$basis_effects_mat)
    expect_identical(fit_missing$betas, fit_null$betas)
    expect_identical(fit_missing$Sigma_space, fit_null$Sigma_space)
    expect_identical(fit_missing$logLik, fit_null$logLik)
})

test_that("a strong prior_precision_space shrinks loadings and Sigma_space", {
    skip_on_cran()
    skip_if_not_installed("autoFRK")
    d <- .make_tiny_cbfm_dataset()
    p <- ncol(d$Bsp)

    fit_null <- CBFM(
        y = d$y, formula = ~ a + b, data = d$data, B_space = d$Bsp,
        family = stats::binomial(),
        control = list(trace = 0, maxit = 5),
        stderrors = FALSE, ncores = 1,
        Sigma_control = list(rank = 2),
        G_control = list(rank = 2)
    )
    fit_strong <- CBFM(
        y = d$y, formula = ~ a + b, data = d$data, B_space = d$Bsp,
        family = stats::binomial(),
        control = list(trace = 0, maxit = 5),
        stderrors = FALSE, ncores = 1,
        Sigma_control = list(rank = 2, prior_precision_space = 1e6 * diag(p)),
        G_control = list(rank = 2)
    )

    # The augmentation invcov ≈ data-driven + 1e6 * I, so cov ≈ 1e-6 * I
    expect_lt(max(abs(fit_strong$Sigma_space)), 1e-4)
    # Basis loadings must contract vs the unconstrained fit
    expect_lt(
        max(abs(fit_strong$basis_effects_mat)),
        0.9 * max(abs(fit_null$basis_effects_mat))
    )
})
