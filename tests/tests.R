# tests/tests.R
# Run automatically by R CMD check. Any error (stopifnot failure) fails the
# check. Structure: one block per procedure, each validated against
# (a) an independent implementation in base R (stats::power.t.test),
# (b) closed-form values computed from the published formulas, and
# (c) internal consistency: power(n) >= target and power(n - 1) < target.

library(clinss)

tol <- 1e-8

# ---------------------------------------------------------------------------
# 1. Two means, superiority, t-test: cross-validate against power.t.test
#    Julious (2004): standardized effect 0.5, 90% power, two-sided 0.05
#    gives 86 per group for the t-test.
# ---------------------------------------------------------------------------
res <- ss_two_means(effect = 5, sd = 10, power = 0.90, alpha = 0.05,
                    hypothesis = "superiority", sides = 2, test = "t")
stopifnot(res$n[["group 1"]] == 86L, res$n[["group 2"]] == 86L,
          res$n_total == 172L)

ref <- stats::power.t.test(n = 86, delta = 5, sd = 10, sig.level = 0.05,
                           type = "two.sample", alternative = "two.sided", strict = TRUE)
stopifnot(abs(res$power_achieved - ref$power) < 1e-6)

# smallest n property
p_at <- function(n) power_two_means(n, n, effect = 5, sd = 10, alpha = 0.05,
                                    hypothesis = "superiority", sides = 2,
                                    test = "t")
stopifnot(p_at(86) >= 0.90, p_at(85) < 0.90)

# power engine agrees with power.t.test across a grid of n
for (n in c(10, 25, 50, 100, 400)) {
  ours   <- power_two_means(n, n, effect = 3, sd = 8, alpha = 0.05,
                            hypothesis = "superiority", sides = 2, test = "t")
  theirs <- stats::power.t.test(n = n, delta = 3, sd = 8,
                                sig.level = 0.05, strict = TRUE)$power
  stopifnot(abs(ours - theirs) < 1e-6)
}

# ---------------------------------------------------------------------------
# 2. Two means, superiority, z-test: closed form
#    n1 = (1 + 1/r) sigma^2 (z_{1-a/2} + z_{1-b})^2 / delta^2
#    Chow, Shao & Wang (2008) eq. 3.2.3 form; Julious (2004) eq. (5).
# ---------------------------------------------------------------------------
res_z <- ss_two_means(effect = 5, sd = 10, power = 0.90, alpha = 0.05,
                      sides = 2, test = "z")
n_closed <- 2 * (10 / 5)^2 * (qnorm(0.975) + qnorm(0.90))^2
stopifnot(res_z$n[["group 1"]] == ceiling(n_closed))   # 85

# ---------------------------------------------------------------------------
# 3. Unequal allocation: ratio = 2 (n2 = 2 n1)
#    n1 = (1 + 1/2) sigma^2 (za + zb)^2 / delta^2
# ---------------------------------------------------------------------------
res_r <- ss_two_means(effect = 5, sd = 10, power = 0.90, alpha = 0.05,
                      sides = 2, test = "z", ratio = 2)
n1_closed <- (1 + 1 / 2) * (10 / 5)^2 * (qnorm(0.975) + qnorm(0.90))^2
stopifnot(res_r$n[["group 1"]] == ceiling(n1_closed),
          res_r$n[["group 2"]] == 2 * res_r$n[["group 1"]])

# ---------------------------------------------------------------------------
# 4. Noninferiority, z-test: closed form with shifted effect
#    Julious (2004), Section 4: n1 = (1+1/r) sigma^2 (za+zb)^2/(d+m)^2, 1-sided
# ---------------------------------------------------------------------------
res_ni <- ss_two_means(effect = 0, sd = 20, margin = 10, alpha = 0.025,
                       power = 0.90, hypothesis = "noninferiority", test = "z")
n_ni <- 2 * (20 / 10)^2 * (qnorm(0.975) + qnorm(0.90))^2
stopifnot(res_ni$n[["group 1"]] == ceiling(n_ni), res_ni$sides == 1L)

# t-test version: internal consistency and dominance over z
res_ni_t <- ss_two_means(effect = 0, sd = 20, margin = 10, alpha = 0.025,
                         power = 0.90, hypothesis = "noninferiority",
                         test = "t")
stopifnot(res_ni_t$n[["group 1"]] >= res_ni$n[["group 1"]])
p_ni <- function(n) power_two_means(n, n, effect = 0, sd = 20, alpha = 0.025,
                                    hypothesis = "noninferiority",
                                    margin = 10, test = "t")
stopifnot(p_ni(res_ni_t$n[["group 1"]]) >= 0.90,
          p_ni(res_ni_t$n[["group 1"]] - 1L) < 0.90)

# ---------------------------------------------------------------------------
# 5. Superiority by a margin equals superiority with effect - margin (z-test)
# ---------------------------------------------------------------------------
a <- ss_two_means(effect = 8, sd = 10, margin = 3, alpha = 0.025,
                  power = 0.80, hypothesis = "superiority_margin", test = "z")
b <- ss_two_means(effect = 5, sd = 10, alpha = 0.025, power = 0.80,
                  sides = 1, hypothesis = "superiority", test = "z")
stopifnot(a$n_total == b$n_total)

# ---------------------------------------------------------------------------
# 6. Equivalence (TOST)
#    Zero true difference, z-test closed form (Julious 2004, Section 6):
#    n1 = (1+1/r) sigma^2 (z_{1-a} + z_{1-b/2})^2 / m^2
# ---------------------------------------------------------------------------
res_eq <- ss_two_means(effect = 0, sd = 10, margin = 5, alpha = 0.05,
                       power = 0.80, hypothesis = "equivalence", test = "z")
n_eq <- 2 * (10 / 5)^2 * (qnorm(0.95) + qnorm(0.90))^2
stopifnot(res_eq$n[["group 1"]] == ceiling(n_eq))

# TOST power is symmetric in the sign of the true effect
p1 <- power_two_means(60, 60, effect =  1, sd = 10, alpha = 0.05,
                      hypothesis = "equivalence", margin = 5, test = "t")
p2 <- power_two_means(60, 60, effect = -1, sd = 10, alpha = 0.05,
                      hypothesis = "equivalence", margin = 5, test = "t")
stopifnot(abs(p1 - p2) < tol)

# TOST power at the margin boundary is at most alpha
pb <- power_two_means(200, 200, effect = 5, sd = 10, alpha = 0.05,
                      hypothesis = "equivalence", margin = 5, test = "z")
stopifnot(pb <= 0.05 + 1e-9)

# ---------------------------------------------------------------------------
# 7. One mean and paired means: cross-validate against power.t.test
# ---------------------------------------------------------------------------
res1 <- ss_one_mean(effect = 5, sd = 10, power = 0.90, alpha = 0.05, sides = 2)
ref1 <- stats::power.t.test(delta = 5, sd = 10, sig.level = 0.05,
                            power = 0.90, type = "one.sample")
stopifnot(res1$n[["subjects"]] == ceiling(ref1$n))

resp <- ss_paired_means(effect = 2, sd_diff = 6, power = 0.80, alpha = 0.05,
                        sides = 2)
refp <- stats::power.t.test(delta = 2, sd = 6, sig.level = 0.05,
                            power = 0.80, type = "paired")
stopifnot(resp$n[["pairs"]] == ceiling(refp$n))

# ---------------------------------------------------------------------------
# 8. Dropout inflation: n_enrolled = ceiling(n_total / (1 - dropout))
# ---------------------------------------------------------------------------
res_d <- ss_two_means(effect = 5, sd = 10, power = 0.90, alpha = 0.05,
                      sides = 2, dropout = 0.15)
stopifnot(res_d$n_enrolled == ceiling(res_d$n_total / 0.85))

# ---------------------------------------------------------------------------
# 9. Result object and methods
# ---------------------------------------------------------------------------
stopifnot(inherits(res, "clinss_result"))
out <- capture.output(print(res))
stopifnot(length(out) > 3)
df <- as.data.frame(res)
stopifnot(is.data.frame(df), df$n_total == 172L)
txt <- report(res)
stopifnot(is.character(txt), grepl("172", txt), grepl("90%", txt))
txt_d <- report(res_d)
stopifnot(grepl("dropout", txt_d))

# ---------------------------------------------------------------------------
# 10. Input validation errors
# ---------------------------------------------------------------------------
expect_error <- function(expr) {
  ok <- tryCatch({ expr; FALSE }, error = function(e) TRUE)
  stopifnot(ok)
}
expect_error(ss_two_means(effect = 0, sd = 10))                    # zero effect
expect_error(ss_two_means(effect = 5, sd = -1))                    # bad sd
expect_error(ss_two_means(effect = 5, sd = 10, power = 1.2))       # bad power
expect_error(ss_two_means(effect = -3, sd = 10, margin = 2,
                          hypothesis = "noninferiority"))          # wrong side
expect_error(ss_two_means(effect = 6, sd = 10, margin = 5,
                          hypothesis = "equivalence"))             # |eff|>=m

cat("All clinss tests passed.\n")

# ===========================================================================
# v0.2.0: Binary endpoints
# ===========================================================================

# ---------------------------------------------------------------------------
# 11. Two proportions, superiority: cross-validate against power.prop.test
#     (base R independent implementation; same pooled-variance formula)
# ---------------------------------------------------------------------------
for (cfg in list(c(0.60, 0.40, 50), c(0.65, 0.45, 80), c(0.30, 0.15, 120))) {
  p1 <- cfg[1]; p2 <- cfg[2]; n <- cfg[3]
  ours   <- power_two_proportions(n1 = n, p1 = p1, p2 = p2, alpha = 0.05,
                                  hypothesis = "superiority", sides = 2)
  theirs <- stats::power.prop.test(n = n, p1 = p1, p2 = p2,
                                   sig.level = 0.05, strict = TRUE)$power
  # strict = TRUE makes power.prop.test include the wrong-tail rejection
  # probability, matching our exact two-sided power.
  stopifnot(abs(ours - theirs) < 1e-6)
}

res_2p <- ss_two_proportions(p1 = 0.60, p2 = 0.40, power = 0.90,
                             alpha = 0.05, sides = 2)
refn <- stats::power.prop.test(p1 = 0.60, p2 = 0.40, power = 0.90,
                               sig.level = 0.05, strict = TRUE)$n
stopifnot(res_2p$n[["group 1"]] == ceiling(refn))

# ---------------------------------------------------------------------------
# 12. Farrington-Manning restricted MLE: closed-form cubic vs independent
#     numerical maximization of the restricted binomial likelihood
# ---------------------------------------------------------------------------
fm_mle_num <- function(p1, p2, n1, n2, d0) {
  negll <- function(p2t) {
    p1t <- p2t + d0
    -(n1 * (p1 * log(p1t) + (1 - p1) * log(1 - p1t)) +
      n2 * (p2 * log(p2t) + (1 - p2) * log(1 - p2t)))
  }
  lo <- max(1e-10, -d0 + 1e-10); hi <- min(1 - 1e-10, 1 - d0 - 1e-10)
  opt <- stats::optimize(negll, c(lo, hi), tol = 1e-12)
  c(p1t = opt$minimum + d0, p2t = opt$minimum)
}

fm_cf <- clinss:::fm_restricted_mle
set.seed(1)
for (i in 1:50) {
  p1 <- runif(1, 0.10, 0.90); p2 <- runif(1, 0.10, 0.90)
  n1 <- sample(20:400, 1);    n2 <- sample(20:400, 1)
  d0 <- runif(1, max(-0.3, -p2 + 0.05), min(0.3, 1 - p2 - 0.05, p1 - 0.05))
  a <- fm_cf(p1, p2, n1, n2, d0)
  b <- fm_mle_num(p1, p2, n1, n2, d0)
  stopifnot(max(abs(a - b)) < 1e-6)
  stopifnot(abs(a["p1t"] - a["p2t"] - d0) < 1e-10)
  stopifnot(all(a > 0), all(a < 1))
}

# d0 = 0 must reduce to the pooled proportion
a0 <- fm_cf(0.6, 0.4, 100, 200, 0)
stopifnot(abs(a0[["p1t"]] - (100 * 0.6 + 200 * 0.4) / 300) < 1e-12)

# ---------------------------------------------------------------------------
# 13. Two proportions, non-inferiority (Farrington-Manning)
# ---------------------------------------------------------------------------
res_ni2 <- ss_two_proportions(p1 = 0.85, p2 = 0.85, margin = 0.10,
                              alpha = 0.025, power = 0.90,
                              hypothesis = "noninferiority")
stopifnot(res_ni2$sides == 1L)
p_at <- function(n) power_two_proportions(n, n, p1 = 0.85, p2 = 0.85,
                                          margin = 0.10, alpha = 0.025,
                                          hypothesis = "noninferiority")
stopifnot(p_at(res_ni2$n[["group 1"]]) >= 0.90,
          p_at(res_ni2$n[["group 1"]] - 1L) < 0.90)

# ---------------------------------------------------------------------------
# 14. Two proportions, equivalence: symmetry and boundary behaviour
# ---------------------------------------------------------------------------
pe1 <- power_two_proportions(300, 300, p1 = 0.72, p2 = 0.70, margin = 0.10,
                             alpha = 0.05, hypothesis = "equivalence")
pe2 <- power_two_proportions(300, 300, p1 = 0.70, p2 = 0.72, margin = 0.10,
                             alpha = 0.05, hypothesis = "equivalence")
stopifnot(abs(pe1 - pe2) < 1e-9)

pb2 <- power_two_proportions(2000, 2000, p1 = 0.80, p2 = 0.70, margin = 0.10,
                             alpha = 0.05, hypothesis = "equivalence")
stopifnot(pb2 <= 0.05 + 1e-9)

# ---------------------------------------------------------------------------
# 15. One proportion, exact test: rejection threshold agrees with binom.test
# ---------------------------------------------------------------------------
res_1p <- ss_one_proportion(p = 0.35, p0 = 0.20, power = 0.80, alpha = 0.05,
                            sides = 1)
n <- res_1p$n[["subjects"]]
k <- clinss:::exact_crit_upper(n, 0.20, 0.05)
# k must reject (p-value <= alpha) while k - 1 must not
stopifnot(stats::binom.test(k, n, 0.20,
                            alternative = "greater")$p.value <= 0.05)
stopifnot(stats::binom.test(k - 1, n, 0.20,
                            alternative = "greater")$p.value > 0.05)
# achieved power matches direct binomial tail computation
stopifnot(abs(res_1p$power_achieved -
              stats::pbinom(k - 1, n, 0.35, lower.tail = FALSE)) < 1e-12)
# first-n property: no smaller n reaches the target
for (m in 2:(n - 1)) {
  stopifnot(power_one_proportion(m, p = 0.35, p0 = 0.20, alpha = 0.05,
                                 hypothesis = "superiority", sides = 1) < 0.80)
}

# ---------------------------------------------------------------------------
# 16. One proportion, z-test: closed-form check
#     n from  ( z_a * sqrt(p0 q0) + z_b * sqrt(p q) )^2 / (p - p0)^2
# ---------------------------------------------------------------------------
res_1z <- ss_one_proportion(p = 0.35, p0 = 0.20, power = 0.80, alpha = 0.05,
                            sides = 1, test = "z")
n_cf <- (qnorm(0.95) * sqrt(0.20 * 0.80) +
         qnorm(0.80) * sqrt(0.35 * 0.65))^2 / 0.15^2
stopifnot(res_1z$n[["subjects"]] == ceiling(n_cf))

# ---------------------------------------------------------------------------
# 17. Binary results: object integrity and report wording
# ---------------------------------------------------------------------------
stopifnot(res_2p$endpoint == "two_proportions",
          res_1p$endpoint == "one_proportion")
txt2 <- report(res_ni2)
stopifnot(grepl("non-inferiority", txt2), grepl("percentage points", txt2),
          grepl("85.0%", txt2))
txt1 <- report(res_1p)
stopifnot(grepl("reference value", txt1), grepl("35.0%", txt1))
df2 <- as.data.frame(res_2p)
stopifnot(df2$n_total == res_2p$n_total)

# ---------------------------------------------------------------------------
# 18. Binary input validation
# ---------------------------------------------------------------------------
expect_error(ss_two_proportions(p1 = 0.5, p2 = 0.5))              # no diff
expect_error(ss_two_proportions(p1 = 1.2, p2 = 0.5))              # bad p
expect_error(ss_two_proportions(p1 = 0.60, p2 = 0.80, margin = 0.10,
                                hypothesis = "noninferiority"))   # wrong side
expect_error(ss_one_proportion(p = 0.50, p0 = 0.10, margin = 0.20,
                               hypothesis = "noninferiority"))    # p0-m < 0
expect_error(power_one_proportion(50, p = 0.5, p0 = 0.95, margin = 0.10,
                                  hypothesis = "superiority_margin"))

cat("All clinss v0.2.0 tests passed.\n")
