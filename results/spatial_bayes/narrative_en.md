# Bayesian unexposed spatial reference


## Methods


The model fits natural-log1p imaging counts with anatomical and anisotropic spatial-basis predictors and proper physical-section nuisance effects. Exact Normal–Inverse-Gamma conjugate updates produce Student-t latent-mean and new-cell intervals. All priors, kernel lengths and candidate genes are fixed; whole-section blocks and their neighboring guards are excluded before learning feature scales, knots or coefficients. Expression-derived cell labels are excluded.


## Results


HY: all 37,811 real cells across 20 physical sections and 10 fixed measured genes were evaluated. Median gene RMSE = 0.140; median empirical 95% predictive coverage = 97.3%, versus 0.0% among positive cells. Median Gaussian predictive mass below zero = 41.1%. Spatial RMSE improves over the anatomy-only baseline for 1/10 genes. These numerical outcomes are retained even when calibration is poor.

ACA: all 14,582 real cells across 18 physical sections and 10 fixed measured genes were evaluated. Median gene RMSE = 0.373; median empirical 95% predictive coverage = 94.0%, versus 66.5% among positive cells. Median Gaussian predictive mass below zero = 38.2%. Spatial RMSE improves over the anatomy-only baseline for 4/10 genes. These numerical outcomes are retained even when calibration is poor.


## Interpretation and limits


The Gaussian likelihood is an inadequate working baseline for positive-cell prediction in these data, not a validated reconstruction or count-generating model. Its fitted means and latent intervals can also be negative. Impossible predictive mass below zero, weak RMSE improvement and poor positive-tail recovery together limit its usefulness, even where overall coverage looks close to 95%. Positive-only interval inclusion is an outcome-conditioned tail diagnostic, not a requirement to achieve nominal 95% coverage within that subset. Credible intervals are conditional on the specified basis/priors and do not establish uncertainty across animals or registration accuracy. Held-out sections quantify within-donor spatial transport only; they do not validate a PM response. Full-reference fitted arrays are distinct from all held-out arrays/metrics. No cell is generated, relocated or given a measured exposure response.


## Sources


[Original atlas](https://doi.org/10.1038/s41586-023-06808-9); [conjugate Gaussian analysis](https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf); [Student-t implementation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.t.html).
