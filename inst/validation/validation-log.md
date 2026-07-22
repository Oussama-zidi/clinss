# clinss validation log

Each entry records how a procedure was validated. Tests enforcing these
comparisons run on every `R CMD check` (see `tests/tests.R`).

## ss_two_means / power_two_means (v0.1.0)

Formulas: Julious (2004), Stat Med 23:1921-1986, Sections 2, 4, 6;
Chow, Shao & Wang (2008), Chapter 3.

| Scenario | Independent benchmark | Result |
|---|---|---|
| Superiority, t-test, effect 5, SD 10, alpha 0.05 two-sided, power 0.90 | `stats::power.t.test` (n = 85.03, round up) and Julious (2004) worked value of 86/group | 86/group, agreement |
| Power at n = 10..400 grid | `stats::power.t.test(strict = TRUE)` | agreement < 1e-6 |
| Superiority, z-test | closed form (1+1/r) s^2 (z_a + z_b)^2 / d^2 | exact match (85/group) |
| Unequal allocation r = 2 | same closed form | exact match |
| Non-inferiority, z-test, margin 10, SD 20 | closed form with shifted effect (Julious Sec. 4) | exact match |
| Superiority by margin m | equals superiority with effect - m (one-sided) | exact match |
| Equivalence (TOST), effect 0, z-test | closed form with z_{1-beta/2} (Julious Sec. 6) | exact match |
| TOST properties | power symmetric in sign of effect; power at boundary <= alpha | verified |
| Smallest-n property | power(n) >= target, power(n-1) < target | verified |

## ss_one_mean / ss_paired_means (v0.1.0)

| Scenario | Independent benchmark | Result |
|---|---|---|
| One mean, t-test | `stats::power.t.test(type = "one.sample")` | agreement |
| Paired means, t-test | `stats::power.t.test(type = "paired")` | agreement |

## Notes

- Our two-sided powers include both rejection tails (exact); this matches
  `power.t.test(strict = TRUE)`. The difference from `strict = FALSE` is
  negligible at all practical powers.

## ss_two_proportions / power_two_proportions (v0.2.0)

Formulas: Farrington & Manning (1990), Stat Med 9:1447-1454; Fleiss, Levin
& Paik (2003), Ch. 4; Chow, Shao & Wang (2008), Ch. 4.

| Scenario | Independent benchmark | Result |
|---|---|---|
| Superiority, pooled z-test, power grid over (p1, p2, n) | `stats::power.prop.test(strict = TRUE)` | agreement < 1e-6 |
| Superiority sample size (0.60 vs 0.40, 90% power) | `stats::power.prop.test(strict = TRUE)` | agreement |
| FM restricted MLE closed-form cubic, 50 random scenarios | direct numerical maximization of the restricted binomial likelihood (`stats::optimize`) | max abs. difference < 1e-6; constraint p1t - p2t = d0 holds to 1e-10 |
| FM restricted MLE at d0 = 0 | pooled proportion (n1 p1 + n2 p2)/(n1 + n2) | exact match |
| Non-inferiority smallest-n property | power(n) >= target, power(n-1) < target | verified |
| Equivalence | power symmetric in sign of the true difference; power at the margin boundary <= alpha | verified |

Monte Carlo simulation (40,000 replicates per condition, seed 2026,
performed at release; script not shipped in tests due to runtime):

| Design | n/group (clinss) | Nominal power | Empirical power (MC SE) | Empirical T1E at boundary (nominal 0.025) |
|---|---|---|---|---|
| p1 = p2 = 0.85, margin 0.10, alpha 0.025 | 276 | 0.900 | 0.9048 (0.0015) | 0.0260 |
| p1 = 0.70, p2 = 0.65, margin 0.15, alpha 0.025 | 86 | 0.802 | 0.7986 (0.0020) | 0.0261 |

The slight liberality of the FM score test in finite samples (~0.026 vs
0.025) is a documented property of the asymptotic score test, not an
implementation artefact.

## ss_one_proportion / power_one_proportion (v0.2.0)

| Scenario | Independent benchmark | Result |
|---|---|---|
| Exact test rejection threshold | `stats::binom.test` p-values: threshold k rejects (p <= alpha), k - 1 does not | verified |
| Exact achieved power | direct binomial tail sum `pbinom(k - 1, n, p, lower.tail = FALSE)` | exact match |
| Exact first-n property | power(m) < target for every m < n | verified exhaustively |
| z-test sample size | closed form (z_a sqrt(p0 q0) + z_b sqrt(p q))^2 / (p - p0)^2 | exact match |

Note: exact binomial power is a non-monotone step function of n
("sawtooth"); the returned n is the smallest reaching the target, which is
standard for exact procedures and documented in the help page.
