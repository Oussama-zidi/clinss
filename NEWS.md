# clinss 0.2.0

* Binary endpoints:
  - `ss_two_proportions()` / `power_two_proportions()`: difference of two
    independent proportions. Superiority uses the pooled z-test (uncorrected
    chi-square); non-inferiority, superiority by a margin, and equivalence
    use the Farrington-Manning score test with the closed-form restricted
    MLE null variance (Farrington & Manning 1990).
  - `ss_one_proportion()` / `power_one_proportion()`: one proportion against
    a reference value, exact binomial test (default) or z approximation,
    all four hypothesis types.
* `clinss_result` objects gained an `endpoint` field; `report()` now
  produces endpoint-specific protocol wording for proportions.
* Validation: the FM restricted MLE closed form is tested against direct
  numerical maximization of the restricted likelihood; superiority powers
  are cross-validated against `stats::power.prop.test()`; exact-test
  rejection thresholds are cross-validated against `stats::binom.test()`;
  operating characteristics confirmed by Monte Carlo simulation (see
  `inst/validation/validation-log.md`).

# clinss 0.1.0

* Initial release: core utilities (normal quantiles, dropout and allocation
  adjustments, validation, rich result objects, protocol text via report()).
* Continuous endpoints: ss_two_means(), ss_one_mean(), ss_paired_means(),
  each supporting superiority, non-inferiority, superiority by a margin, and
  equivalence (TOST), with exact noncentral-t or normal-approximation
  calculations, unequal allocation, and dropout inflation.
* All procedures validated against stats::power.t.test() and closed-form
  published formulas (Julious 2004; Chow, Shao & Wang 2008).
