"""Independent posterior-sampling check for analytic interval transport."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("spatial_bayes", ROOT / "scripts/analyze_spatial_bayes.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class PosteriorPredictiveSampling(unittest.TestCase):
    def test_analytic_intervals_match_independent_posterior_draws(self):
        rng = np.random.default_rng(20261020)
        coords = rng.normal(size=(300, 3))
        anatomy = np.array(["A", "B", "C"] * 100)
        sections = np.array(["S" + str(i % 6) for i in range(300)])
        response = np.log1p(rng.poisson(np.exp(.2 + .2 * coords[:, 0])))[:, None]
        design = m.Design(coords, anatomy, sections)
        x, _ = design.transform(coords, anatomy, sections, True)
        fitted = m.fit(x, response, design.prior)
        xx, unknown = design.transform(coords[:2], np.array(["A", "UNSEEN"]), np.array(["NEW", "NEW"]))
        pred = m.predict(fitted, xx, unknown)
        draws = 60000
        variance = fitted["b"][0] / rng.gamma(fitted["a"], 1, size=draws)
        beta = fitted["mean"][:, 0] + (rng.normal(size=(draws, len(design.prior))) @ np.linalg.cholesky(fitted["V"]).T) * np.sqrt(variance[:, None])
        latent = beta @ xx.T
        # Unknown anatomical intercept is integrated with its specified proper prior.
        latent += rng.normal(size=(draws, 2)) * np.sqrt(variance[:, None] * unknown[None, :])
        section_variance = m.methods()["prior"]["section_variance_multiplier"]
        predictive = latent + rng.normal(size=(draws, 2)) * np.sqrt(variance[:, None] * (1 + section_variance))
        for j in range(2):
            covered = ((predictive[:, j] >= pred["predictive_lower"][j, 0]) & (predictive[:, j] <= pred["predictive_upper"][j, 0])).mean()
            latent_covered = ((latent[:, j] >= pred["ci_lower"][j, 0]) & (latent[:, j] <= pred["ci_upper"][j, 0])).mean()
            self.assertAlmostEqual(covered, .95, delta=.008)
            self.assertAlmostEqual(latent_covered, .95, delta=.008)
            self.assertAlmostEqual(float((predictive[:, j] < 0).mean()), pred["negative_mass"][j, 0], delta=.008)

if __name__ == "__main__":
    unittest.main()
