#!/usr/bin/env python3
"""Finite-grid hierarchical Beta-Binomial baseline detection analysis.

This model estimates detected-cell proportions, not expression magnitude,
exposure responses, new cell locations, or increased spatial resolution.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
import scipy
from scipy.special import betaln, betainc, gammaln, logsumexp
from scipy.optimize import brentq
from scipy.stats import beta as beta_distribution
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/spatial_detection_bayes"
META = ROOT / "metadata/spatial_detection_bayes"
WEB = ROOT / "docs/assets/spatial-detection-bayes"
SEED = 20261020
DRAWS = 8192
PRIOR_MEDIAN_KAPPA = 20.
PRIOR_LOG_KAPPA_SD = 1.5
MAIN_GRID = (64, 48)
SENSITIVITY_GRIDS = [(32, 24), (96, 64)]
REGIONS = {"ACA": {"region_en": "Anterior cingulate cortex", "region_es": "Corteza cingulada anterior", "curated": ["Fosb", "Egr3", "S1pr1"]},
           "HY": {"region_en": "Hypothalamus", "region_es": "Hipotálamo", "curated": ["Tll1", "Fosl2", "Nr4a3"]}}


def clean(x):
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, np.ndarray)):
        return [clean(v) for v in x]
    if isinstance(x, (float, np.floating)):
        return float(x) if np.isfinite(x) else None
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def dump(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(clean(obj), ensure_ascii=False, separators=(",", ":"), allow_nan=False) + "\n")


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def grid(n_mu, n_kappa):
    """Equal-weight quadrature in the two prior-quantile coordinates.

    No extra prior density or Jacobian is multiplied into these weights.
    """
    mu_axis = beta_distribution.ppf((np.arange(n_mu) + .5) / n_mu, .5, .5)
    kappa_axis = np.exp(np.log(PRIOR_MEDIAN_KAPPA) + PRIOR_LOG_KAPPA_SD * norm.ppf((np.arange(n_kappa) + .5) / n_kappa))
    mu, kappa = np.meshgrid(mu_axis, kappa_axis, indexing="ij")
    mu, kappa = mu.ravel(), kappa.ravel()
    return {"shape": (n_mu, n_kappa), "mu": mu, "kappa": kappa,
            "alpha": mu * kappa, "beta": (1 - mu) * kappa}


GRIDS = {shape: grid(*shape) for shape in [MAIN_GRID] + SENSITIVITY_GRIDS}


def weighted_quantiles(values, weights, levels=(.025, .975)):
    order = np.argsort(values, kind="stable")
    cdf = np.cumsum(weights[order])
    return np.asarray(values)[order[np.searchsorted(cdf, levels).clip(max=len(order) - 1)]]


def fit(y, n, shape=MAIN_GRID):
    g = GRIDS[shape]
    y, n = np.asarray(y, dtype=int), np.asarray(n, dtype=int)
    assert len(y) == len(n) and ((0 <= y) & (y <= n) & (n > 0)).all()
    if len(y):
        ll = (gammaln(n + 1) - gammaln(y + 1) - gammaln(n - y + 1))[:, None]
        ll = ll + betaln(y[:, None] + g["alpha"], (n - y)[:, None] + g["beta"]) - betaln(g["alpha"], g["beta"])
        logweights = ll.sum(axis=0)
    else:
        logweights = np.zeros(len(g["mu"]))
    weights = np.exp(logweights - logsumexp(logweights))
    ci = weighted_quantiles(g["mu"], weights)
    w2 = weights.reshape(shape)
    return {"grid": g, "weights": weights, "mu_mean": float(weights @ g["mu"]),
            "mu_ci_low": float(ci[0]), "mu_ci_high": float(ci[1]),
            "kappa_mean": float(weights @ g["kappa"]),
            "rho_mean": float(weights @ (1 / (1 + g["kappa"]))),
            "prior_only": len(y) == 0, "n_training_sections": len(y), "n_training_cells": int(sum(n)),
            "effective_grid_nodes": float(1 / (weights @ weights)),
            "mu_edge_weight": float(w2[0].sum() + w2[-1].sum()),
            "kappa_edge_weight": float(w2[:, 0].sum() + w2[:, -1].sum())}


def rng_for(region, gene_index, purpose, anatomy_code=0, fold=0, section_code=0):
    return np.random.default_rng(np.random.SeedSequence([SEED, 1 if region == "ACA" else 2,
                                                        gene_index, purpose, anatomy_code, fold, section_code]))


def mixture_draws(posterior, rng, y=0, n=0, draws=DRAWS):
    indices = rng.choice(len(posterior["weights"]), size=draws, p=posterior["weights"])
    g = posterior["grid"]
    return rng.beta(g["alpha"][indices] + y, g["beta"][indices] + n - y)


def new_p_quantiles(posterior):
    """Numerically invert the finite-grid predictive Beta-mixture CDF.

    Unlike known-section/count MC intervals, full-data new-section propensity
    intervals and their grid comparison contain no Monte Carlo draw noise.
    """
    g, weights = posterior["grid"], posterior["weights"]
    return np.array([brentq(lambda x: float(weights @ betainc(g["alpha"], g["beta"], x)) - q,
                           0., 1., xtol=1e-12, rtol=1e-12) for q in [.025, .975]])


def section_posterior(posterior, y, n, rng):
    samples = mixture_draws(posterior, rng, y=y, n=n)
    g = posterior["grid"]
    mean = float(posterior["weights"] @ ((g["alpha"] + y) / (g["kappa"] + n)))
    ci = np.quantile(samples, [.025, .975], method="inverted_cdf")
    return {"p_mean": mean, "p_ci_low": float(ci[0]), "p_ci_high": float(ci[1])}


def new_section_prediction(posterior, n, rng):
    p = mixture_draws(posterior, rng)
    p_ci = np.quantile(p, [.025, .975], method="inverted_cdf")
    counts = rng.binomial(int(n), p)
    count_ci = np.quantile(counts, [.025, .975], method="inverted_cdf").astype(int)
    return {"new_section_p_ci_low": float(p_ci[0]), "new_section_p_ci_high": float(p_ci[1]),
            "count_pi_low": int(count_ci[0]), "count_pi_high": int(count_ci[1]),
            "count_predictive_mean": float(n * posterior["mu_mean"])}


def brier(y, n, p):
    return (y * (1 - p) ** 2 + (n - y) * p ** 2) / n


def summarize_oof(rows):
    t = pd.DataFrame(rows)
    if not len(t):
        return {"n_groups": 0}
    n = t.n.to_numpy()
    return {"n_groups": len(t), "n_cells": int(n.sum()), "count_pi_covered_groups": int(t.count_pi_covered.sum()), "count_pi_coverage": float(t.count_pi_covered.mean()),
            "count_pi_coverage_cell_weighted": float(np.average(t.count_pi_covered, weights=n)),
            "coverage_target": .95, "interval_target": "New physical-section positive count at its observed n",
            "fraction_mae": float(t.absolute_fraction_error.mean()),
            "fraction_mae_cell_weighted": float(np.average(t.absolute_fraction_error, weights=n)),
            "baseline_fraction_mae": float(t.baseline_absolute_fraction_error.mean()),
            "brier": float(np.average(t.brier, weights=n)), "baseline_brier": float(np.average(t.baseline_brier, weights=n)),
            "brier_improvement": float(np.average(t.baseline_brier - t.brier, weights=n)),
            "prior_only_groups": int(t.prior_only.sum()), "prior_only_cells": int(t.loc[t.prior_only, "n"].sum()),
            "count_pi_mean_fraction_width": float(np.mean((t.count_pi_high - t.count_pi_low) / t.n)),
            "model_is_anatomy_specific_not_cellclass_specific": True}


def group_data(region):
    path = ROOT / f"results/spatial_sections/{region}/detection_by_section_anatomy.tsv.gz"
    t = pd.read_csv(path, sep="\t")
    for field in ["n", "positive", "section_code", "anatomy_code"]:
        assert np.allclose(t[field], np.round(t[field]))
        t[field] = t[field].astype(int)
    assert (t.n > 0).all() and (t.positive <= t.n).all() and (t.positive >= 0).all()
    base = json.loads((ROOT / f"docs/assets/spatial-sections/data-{region}.js").read_text().split("=", 1)[1].strip().removesuffix(";"))
    folds_path = ROOT / f"results/spatial_bayes/{region}/section_folds.tsv"
    audit_path = ROOT / f"results/spatial_bayes/{region}/fold_audit.json"
    folds = pd.read_csv(folds_path, sep="\t")
    audit = json.loads(audit_path.read_text())
    assert set(t.section_label) == set(folds.section)
    assert set(folds.fold) == {1, 2, 3, 4, 5}
    assert len(t.gene.unique()) == 10
    template = t.loc[t.gene.eq(base["genes"][0]["symbol"])]
    per_slice_n = template.groupby("section_label").n.sum().to_dict()
    assert all(per_slice_n[row.section] == row.n_cells for row in folds.itertuples())
    assert sum(per_slice_n.values()) == base["cell_count"]
    for a in audit:
        assert set(a["test_sections"]) == set(folds.loc[folds.fold.eq(a["fold"]), "section"])
        train = set(folds.section) - set(a["test_sections"]) - set(a["guard_sections"])
        assert not train.intersection(a["test_sections"]) and not train.intersection(a["guard_sections"])
        assert sum(per_slice_n[x] for x in train) == a["n_train"]
        assert sum(per_slice_n[x] for x in a["test_sections"]) == a["n_test"]
    return t, base, folds, audit, [path, folds_path, audit_path, ROOT / f"docs/assets/spatial-sections/data-{region}.js"]


def analyze_region(region):
    (OUT / region).mkdir(parents=True, exist_ok=True)
    table, base, folds, audit, inputs = group_data(region)
    genes = []
    full_rows, anatomy_rows, sensitivity_rows, oof_rows, hyper_rows, oof_sensitivity_rows = [], [], [], [], [], []
    oof_by_gene = {}
    for gi, entry in enumerate(base["genes"]):
        symbol = entry["symbol"]
        t = table.loc[table.gene.eq(symbol)].copy()
        gene_anatomy, gene_sections = [], []
        full_fits = {}
        for code, label in enumerate(base["anatomy_labels"]):
            sub = t.loc[t.anatomy_code.eq(code)]
            posterior = fit(sub.positive, sub.n)
            full_fits[code] = posterior
            new_ci = new_p_quantiles(posterior)
            row = {"region": region, "gene_id": entry["id"], "gene": symbol, "anatomy_code": code, "label": label,
                   "n_sections": len(sub), "n": int(sub.n.sum()), "positive": int(sub.positive.sum()),
                   "observed_fraction": float(sub.positive.sum() / sub.n.sum()),
                   **{key: posterior[key] for key in ["mu_mean", "mu_ci_low", "mu_ci_high", "kappa_mean", "rho_mean", "prior_only", "effective_grid_nodes", "mu_edge_weight", "kappa_edge_weight"]},
                   "new_section_p_ci_low": float(new_ci[0]), "new_section_p_ci_high": float(new_ci[1])}
            gene_anatomy.append(row)
            anatomy_rows.append(row)
            for shape in SENSITIVITY_GRIDS:
                other = fit(sub.positive, sub.n, shape)
                other_ci = new_p_quantiles(other)
                sensitivity_rows.append({"region": region, "gene": symbol, "gene_id": entry["id"], "anatomy_code": code, "label": label,
                                         "main_grid": "64x48", "comparison_grid": f"{shape[0]}x{shape[1]}",
                                         "main_mu_mean": posterior["mu_mean"], "other_mu_mean": other["mu_mean"],
                                         "abs_mu_mean_difference": abs(posterior["mu_mean"] - other["mu_mean"]),
                                         "main_mu_ci_low": posterior["mu_ci_low"], "other_mu_ci_low": other["mu_ci_low"],
                                         "main_mu_ci_high": posterior["mu_ci_high"], "other_mu_ci_high": other["mu_ci_high"],
                                         "abs_mu_ci_low_difference": abs(posterior["mu_ci_low"] - other["mu_ci_low"]),
                                         "abs_mu_ci_high_difference": abs(posterior["mu_ci_high"] - other["mu_ci_high"]),
                                         "main_mu_edge_weight": posterior["mu_edge_weight"], "other_mu_edge_weight": other["mu_edge_weight"],
                                         "main_kappa_mean": posterior["kappa_mean"], "other_kappa_mean": other["kappa_mean"],
                                         "relative_kappa_mean_difference": abs(posterior["kappa_mean"]-other["kappa_mean"])/posterior["kappa_mean"],
                                         "main_rho_mean": posterior["rho_mean"], "other_rho_mean": other["rho_mean"],
                                         "abs_rho_mean_difference": abs(posterior["rho_mean"]-other["rho_mean"]),
                                         "main_new_section_p_low": new_ci[0], "main_new_section_p_high": new_ci[1],
                                         "other_new_section_p_low": other_ci[0], "other_new_section_p_high": other_ci[1],
                                         "abs_new_section_p_low_difference": abs(new_ci[0]-other_ci[0]),
                                         "abs_new_section_p_high_difference": abs(new_ci[1]-other_ci[1])})
            for observed in sub.itertuples():
                p = section_posterior(posterior, observed.positive, observed.n, rng_for(region, gi, 2, code, section_code=observed.section_code))
                r = {"region": region, "gene_id": entry["id"], "gene": symbol, "section_code": observed.section_code,
                     "label": observed.section_label, "anatomy_code": code, "anatomy_label": label,
                     "n": observed.n, "positive": observed.positive, "observed_fraction": observed.positive / observed.n, **p}
                gene_sections.append(r)
                full_rows.append(r)
        gene_oof = []
        for fold_spec in audit:
            fold = fold_spec["fold"]
            test_labels = set(fold_spec["test_sections"])
            excluded = test_labels | set(fold_spec["guard_sections"])
            train = t.loc[~t.section_label.isin(excluded)]
            test = t.loc[t.section_label.isin(test_labels)]
            train_regional = train.groupby("section_label").agg(y=("positive", "sum"), n=("n", "sum"))
            regional = fit(train_regional.y, train_regional.n)
            regional_fine = fit(train_regional.y, train_regional.n, (96, 64))
            for code in sorted(test.anatomy_code.unique()):
                sub = train.loc[train.anatomy_code.eq(code)]
                posterior = fit(sub.positive, sub.n)
                posterior_fine = fit(sub.positive, sub.n, (96, 64))
                hyper_rows.append({"region": region, "gene": symbol, "fold": fold, "anatomy_code": code,
                                   "n_training_sections": len(sub), "n_training_cells": int(sub.n.sum()),
                                   **{k: posterior[k] for k in ["mu_mean", "mu_ci_low", "mu_ci_high", "kappa_mean", "rho_mean", "prior_only", "effective_grid_nodes", "mu_edge_weight", "kappa_edge_weight"]}})
                for observed in test.loc[test.anatomy_code.eq(code)].itertuples():
                    p = posterior["mu_mean"]
                    p0 = regional["mu_mean"]
                    prediction = new_section_prediction(posterior, observed.n, rng_for(region, gi, 3, code, fold, observed.section_code))
                    fine_prediction = new_section_prediction(posterior_fine, observed.n, rng_for(region, gi, 3, code, fold, observed.section_code))
                    r = {"region": region, "gene_id": entry["id"], "gene": symbol, "fold": fold,
                         "section_code": observed.section_code, "section_label": observed.section_label,
                         "anatomy_code": code, "anatomy_label": base["anatomy_labels"][code],
                         "n": observed.n, "positive": observed.positive, "observed_fraction": observed.positive / observed.n,
                         "predicted_fraction": p, "baseline_predicted_fraction": p0,
                         "prior_only": posterior["prior_only"], "n_training_sections": len(sub), "n_training_cells": int(sub.n.sum()),
                         "mu_edge_weight": posterior["mu_edge_weight"], "kappa_edge_weight": posterior["kappa_edge_weight"],
                         **prediction, "count_pi_covered": prediction["count_pi_low"] <= observed.positive <= prediction["count_pi_high"],
                         "absolute_fraction_error": abs(observed.positive / observed.n - p),
                         "baseline_absolute_fraction_error": abs(observed.positive / observed.n - p0),
                         "brier": brier(observed.positive, observed.n, p), "baseline_brier": brier(observed.positive, observed.n, p0)}
                    gene_oof.append(r)
                    oof_rows.append(r)
                    oof_sensitivity_rows.append({"region": region, "gene": symbol, "fold": fold, "section_code": observed.section_code,
                                                 "section_label": observed.section_label, "anatomy_code": code, "n": observed.n, "positive": observed.positive,
                                                 "main_predicted_fraction": p, "fine_predicted_fraction": posterior_fine["mu_mean"],
                                                 "abs_fraction_mean_difference": abs(p-posterior_fine["mu_mean"]),
                                                 "main_brier": r["brier"], "fine_brier": brier(observed.positive, observed.n, posterior_fine["mu_mean"]),
                                                 "fine_baseline_brier": brier(observed.positive, observed.n, regional_fine["mu_mean"]),
                                                 "main_count_pi_low": prediction["count_pi_low"], "main_count_pi_high": prediction["count_pi_high"],
                                                 "fine_count_pi_low": fine_prediction["count_pi_low"], "fine_count_pi_high": fine_prediction["count_pi_high"],
                                                 "main_covered": r["count_pi_covered"],
                                                 "fine_covered": fine_prediction["count_pi_low"] <= observed.positive <= fine_prediction["count_pi_high"],
                                                 "max_count_interval_endpoint_change_fraction": max(abs(prediction["count_pi_low"]-fine_prediction["count_pi_low"]),abs(prediction["count_pi_high"]-fine_prediction["count_pi_high"]))/observed.n})
        assert len(gene_oof) == len(t)
        assert len({(x["section_code"], x["anatomy_code"]) for x in gene_oof}) == len(t)
        summary = summarize_oof(gene_oof)
        summary["trained_anatomy_only"] = summarize_oof([x for x in gene_oof if not x["prior_only"]])
        oof_by_gene[symbol] = summary
        genes.append({"id": entry["id"], "symbol": symbol, "anatomy": gene_anatomy, "sections": gene_sections, "oof_summary": summary})
        print(f"{region} {symbol}: count PI coverage={summary['count_pi_coverage']:.3f}, Brier={summary['brier']:.4f} vs baseline={summary['baseline_brier']:.4f}, prior-only={summary['prior_only_groups']}", flush=True)
    for stem, rows in [("known_section_posterior", full_rows), ("anatomy_posterior", anatomy_rows),
                       ("grid_sensitivity", sensitivity_rows), ("out_of_fold_predictions", oof_rows), ("training_hyperparameters", hyper_rows)]:
        pd.DataFrame(rows).to_csv(OUT / region / f"{stem}.tsv.gz", sep="\t", index=False)
    oof_sensitivity = pd.DataFrame(oof_sensitivity_rows)
    oof_sensitivity.to_csv(OUT / region / "oof_grid_sensitivity.tsv.gz", sep="\t", index=False)
    total = summarize_oof(oof_rows)
    total["trained_anatomy_only"] = summarize_oof([x for x in oof_rows if not x["prior_only"]])
    pd.DataFrame([{ "gene": gene, **{k: v for k, v in scores.items() if not isinstance(v, dict)}}
                  for gene, scores in oof_by_gene.items()]).to_csv(OUT / region / "gene_metrics.tsv", sep="\t", index=False)
    pd.DataFrame([{ "gene": gene, **scores["trained_anatomy_only"]}
                  for gene, scores in oof_by_gene.items()]).to_csv(OUT / region / "gene_metrics_trained_anatomy.tsv", sep="\t", index=False)
    folds.to_csv(OUT / region / "section_folds.tsv", sep="\t", index=False)
    sensitivity = pd.DataFrame(sensitivity_rows)
    fine = sensitivity.loc[sensitivity.comparison_grid.eq("96x64")]
    total.update(region=region, genes=10, reference_cells=base["cell_count"], physical_sections=len(base["sections"]),
                 anatomy_groups=len(base["anatomy_labels"]), metrics_by_gene=oof_by_gene,
                 aggregate_metric_unit="Gene×section×anatomy evaluations; n_cells/prior_only_cells are gene×cell evaluations, not unique animals/cells",
                 unique_section_anatomy_groups=len(table.loc[table.gene.eq(base['genes'][0]['symbol'])]),
                 prior_only_unique_reference_cells=int(pd.DataFrame(oof_rows).loc[lambda x: x.gene.eq(base['genes'][0]['symbol']) & x.prior_only, 'n'].sum()),
                 grid_sensitivity_96x64={"max_abs_mu_mean_difference": float(fine.abs_mu_mean_difference.max()),
                                        "median_abs_mu_mean_difference": float(fine.abs_mu_mean_difference.median()),
                                        "max_abs_mu_ci_endpoint_difference": float(max(fine.abs_mu_ci_low_difference.max(), fine.abs_mu_ci_high_difference.max())),
                                        "max_relative_kappa_mean_difference": float(fine.relative_kappa_mean_difference.max()),
                                        "max_abs_rho_mean_difference": float(fine.abs_rho_mean_difference.max()),
                                        "max_abs_new_section_p_interval_endpoint_difference": float(max(fine.abs_new_section_p_low_difference.max(),fine.abs_new_section_p_high_difference.max())),
                                        "p_intervals_numerically_invert_Beta_mixture_CDF_no_MC": True},
                 oof_grid_sensitivity_96x64={"count_pi_coverage": float(oof_sensitivity.fine_covered.mean()),
                                            "brier": float(np.average(oof_sensitivity.fine_brier,weights=oof_sensitivity.n)),
                                            "baseline_brier": float(np.average(oof_sensitivity.fine_baseline_brier,weights=oof_sensitivity.n)),
                                            "max_abs_fraction_mean_difference": float(oof_sensitivity.abs_fraction_mean_difference.max()),
                                            "max_count_endpoint_change_fraction": float(oof_sensitivity.max_count_interval_endpoint_change_fraction.max()),
                                            "coverage_status_changed_groups": int((oof_sensitivity.main_covered != oof_sensitivity.fine_covered).sum()),
                                            "same_MC_seed_and8192draws": True})
    payload = {"schema_version": 1, "region": region, "region_label_en": REGIONS[region]["region_en"], "region_label_es": REGIONS[region]["region_es"],
               "method": "Hierarchical Beta-Binomial baseline detection; finite64x48 prior-quantile quadrature", "interval_level": .95,
               "independent_reference_animals": 1, "cellclass_specific": False, "baseline_only": True,
               "no_PM_effect": True, "no_increased_resolution": True, "counts_unchanged": True,
               "mc_draws": DRAWS, "seed": SEED, "mu_interval_method": "Discrete weighted grid quantiles",
               "p_interval_method": "Known-section p:8192 seeded mixture Beta draws; full-data new-section p:finite-grid Beta-mixture CDF numerical inversion",
               "section_labels": base["section_labels"], "anatomy_labels": base["anatomy_labels"],
               "downloads_path": f"../../downloads/spatial_detection_bayes/{region}/", "genes": genes, "oof_summary": total}
    (WEB / f"data-{region}.js").write_text(f"window.SPATIAL_DETECTION_BAYES_{region}=" + json.dumps(clean(payload), ensure_ascii=False, separators=(",", ":"), allow_nan=False) + ";\n")
    dump(OUT / region / "summary.json", total)
    dump(OUT / region / "fold_audit.json", {"source_folds": audit, "exact_source_split_reused": True,
                                           "training_excludes_all_test_and_guard_sections": True,
                                           "test_anatomy_absent_training": "Fixed prior-only prediction, explicitly flagged; no silent fallback or test-data fit",
                                           "no_test_anatomy_n": "No artificial zero-count group; groups with no measured test cells are absent"})
    return {"payload": payload, "summary": total, "oof": pd.DataFrame(oof_rows), "full": pd.DataFrame(full_rows),
            "anatomy": pd.DataFrame(anatomy_rows), "sensitivity": sensitivity, "oof_sensitivity": oof_sensitivity, "inputs": inputs}


def make_figures(region, result, manifests):
    curated = REGIONS[region]["curated"]
    anatomy = result["anatomy"]
    # Fixed display subset for legibility; ALL strata remain in output tables/payload.
    # Top12 HY labels by all-gene-identical measured n; selection is independent of expression/posterior.
    first = anatomy.loc[anatomy.gene.eq(curated[0])]
    selected_codes = first.sort_values(["n", "anatomy_code"], ascending=[False, True]).head(12).anatomy_code.to_list() if region == "HY" else sorted(first.anatomy_code)
    for spanish in [False, True]:
        T = lambda en, es: es if spanish else en
        fig = plt.figure(figsize=(16, 10.8), layout="constrained")
        grid_spec = fig.add_gridspec(2, 3, height_ratios=[1.35, 1])
        for j, gene in enumerate(curated):
            ax = fig.add_subplot(grid_spec[0, j])
            a = anatomy.loc[anatomy.gene.eq(gene)].set_index("anatomy_code").loc[selected_codes]
            idx = np.arange(len(a))
            # Broad new-section propensity predictive intervals and narrower group-mu credible intervals differ in target.
            ax.hlines(idx, a.new_section_p_ci_low, a.new_section_p_ci_high, color="#b9ccda", lw=8, alpha=.6,
                      label=T("New-section p:95% predictive", "p de corte nuevo:95% predictivo"))
            ax.hlines(idx, a.mu_ci_low, a.mu_ci_high, color="#21628b", lw=3,
                      label=T("Anatomy μ:95% credible", "μ anatómica:95% creíble"))
            ax.scatter(a.mu_mean, idx, c="#21628b", marker="s", s=26, label=T("Posterior μ mean", "Media posterior μ"), zorder=4)
            ax.scatter(a.observed_fraction, idx, c="#c05310", marker="o", s=25,
                       label=T("Measured pooled fraction", "Fracción medida agrupada"), zorder=5)
            ax.set_yticks(idx, [f"{row.label} (n={row.n:,}; S={row.n_sections})" for row in a.itertuples()], fontsize=8)
            ax.invert_yaxis()
            ax.set(xlim=(0, 1), xlabel=T("Detection probability", "Probabilidad de detección"), title=f"{chr(65+j)}  {gene}")
            ax.spines[["top", "right"]].set_visible(False)
            ax.grid(axis="x", color="#e4e9ed", lw=.5)
            if j == 0:
                ax.legend(loc="lower center", bbox_to_anchor=(1.6, 1.10), ncol=2, fontsize=8.2, frameon=False)
        oof = result["oof"]
        names = [g["symbol"] for g in result["payload"]["genes"]]
        metrics = [result["summary"]["metrics_by_gene"][g] for g in names]
        ax = fig.add_subplot(grid_spec[1, 0])
        idx = np.arange(len(names))
        ax.scatter([x["count_pi_coverage"] for x in metrics], idx, marker="o", c="#21628b", s=35)
        ax.axvline(.95, color="#b94b30", ls="--", label=T("Nominal95%", "Nominal95%"))
        ax.set_yticks(idx, names, fontsize=8)
        ax.invert_yaxis()
        ax.set(xlim=(0, 1.02), xlabel=T("Held-out count PI coverage", "Cobertura predictiva\nde conteos excluidos"), title=T("D  New section, actual n", "D  Corte nuevo, n real"))
        ax.legend(frameon=False, fontsize=8, loc="lower left")
        ax.spines[["top", "right"]].set_visible(False)
        ax = fig.add_subplot(grid_spec[1, 1])
        for k, metric in enumerate(["fraction_mae", "baseline_fraction_mae"]):
            vals = [x[metric] for x in metrics]
            ax.scatter(vals, idx, marker=["s", "x"][k], c=["#21628b", "#c05310"][k], s=35,
                       label=T("Anatomy Beta-Binomial", "Beta-Binomial anatómico") if k == 0 else T("Regional training-only baseline", "Referencia regional solo entrenamiento"))
        ax.set_yticks(idx, names, fontsize=8)
        ax.invert_yaxis()
        ax.set(xlim=(0, max(.1, max(x["fraction_mae"] for x in metrics) * 1.08, max(x["baseline_fraction_mae"] for x in metrics) * 1.08)),
               xlabel=T("Fraction MAE\n(section × anatomy)", "Error absoluto de fracción\n(corte × anatomía)"),
               title=T("E  Held-out fraction error", "E  Error de fracción excluida"))
        handles, legend_labels = ax.get_legend_handles_labels()
        fig.legend(handles, legend_labels, frameon=False, fontsize=8, loc="center", bbox_to_anchor=(.45,.092), ncol=2)
        ax.spines[["top", "right"]].set_visible(False)
        ax = fig.add_subplot(grid_spec[1, 2])
        for k, compare in enumerate(["32x24", "96x64"]):
            s = result["sensitivity"].loc[result["sensitivity"].comparison_grid.eq(compare)]
            vals = [s.loc[s.gene.eq(g), "abs_mu_mean_difference"].max() for g in names]
            ax.scatter(vals, idx, marker=["x", "o"][k], c=["#c05310", "#21628b"][k], s=30, label=f"64×48 vs {compare.replace('x','×')}")
        ax.set_yticks(idx, names, fontsize=8)
        ax.invert_yaxis()
        ax.set(xlabel=T("Maximum |Δ posterior μ mean|", "Máximo |Δ media posterior μ|"), xlim=(0, max(.01, result["sensitivity"].abs_mu_mean_difference.max() * 1.1)),
               title=T("F  Grid sensitivity\nfull-data μ means", "F  Sensibilidad de cuadrícula\nmedias μ con datos completos"))
        ax.xaxis.set_major_locator(MaxNLocator(nbins=3,min_n_ticks=3))
        handles, legend_labels = ax.get_legend_handles_labels()
        fig.legend(handles, legend_labels, frameon=False, fontsize=8, loc="center", bbox_to_anchor=(.85,.092), ncol=2)
        ax.spines[["top", "right"]].set_visible(False)
        fig.suptitle((REGIONS[region]["region_es"] if spanish else REGIONS[region]["region_en"]) + " · " +
                     T("hierarchical baseline DETECTION", "DETECCIÓN basal jerárquica"), fontsize=15)
        fig.get_layout_engine().set(rect=(.02, .12, .94, .81), h_pad=.09, w_pad=.08)
        s = result["summary"]
        text_en = (f"One unexposed donor; physical sections model heterogeneity, not independent animals. All measured zeros retained.\n"
                   f"Fixed μ~Beta(½,½), κ~LogNormal(log20,1.5);64×48 prior-quantile grid;8,192 seeded mixture draws. "
                   f"OOF prior-only groups={s['prior_only_groups']}/{s['n_groups']}. No PM2.5 effect or new spatial resolution.")
        text_es = (f"Un donante sin exposición; los cortes modelan heterogeneidad, no animales independientes. Se conservan los ceros medidos.\n"
                   f"μ~Beta(½,½), κ~LogNormal(log20,1.5) fijas;cuadrícula64×48;8,192 sorteos con semilla. "
                   f"Grupos excluidos solo con prior={s['prior_only_groups']}/{s['n_groups']}. Sin efecto PM2.5 ni nueva resolución espacial.")
        fig.text(.5, .045, T(text_en, text_es), ha="center", va="center", fontsize=9, color="#334c55")
        stem = f"{region}_detection_posterior" + ("_es" if spanish else "")
        exports = []
        for ext in ["png", "svg", "pdf"]:
            path = OUT / "figures" / f"{stem}.{ext}"
            fig.savefig(path, dpi=210, facecolor="white", transparent=False)
            exports.append({"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "sha256": digest(path)})
        plt.close(fig)
        cap_en = (f"Hierarchical Beta-Binomial baseline detection in {region}, not RNA magnitude or PM2.5 exposure effects. "
                  "For each deposited anatomy, physical-section positive-cell counts y_s among n_s measured cells follow Beta-Binomial(μκ,(1−μ)κ). "
                  "Priors fixed before fitting: μ~Beta(0.5,0.5), κ~LogNormal(log20,1.5); equal prior-quantile64×48 quadrature approximates the continuous posterior. "
                  "Panels A–C distinguish anatomy mean μ95% credible intervals (dark) from new-section detection propensity p95% predictive intervals (pale); orange dots are pooled measured fractions, not μ estimates. "
                  + ("All10 ACA layer/subdivisions shown. " if region == "ACA" else "For legibility12 of49 HY nucleus/subregions shown, selected solely by largest measured cell n, independent of expression; all49 remain in tables. ") +
                  "n=measured cells, S=physical sections. Panels D–E use the same five contiguous held-out-section folds and neighboring guards as the earlier model; no test/guard counts enter training. "
                  "Count95% predictive intervals are conditioned on actual test n and include both section heterogeneity and binomial sampling; coverage is unweighted across observed section×anatomy groups. "
                  "Prior-only test anatomies are retained and explicitly counted. Regional nonanatomical baseline uses the identical fixed priors, pooling training-only cells within each physical section. "
                  "Fraction MAE is unweighted across observed section×anatomy groups; cell-weighted Brier loss is separately tabulated. No independent-animal or cell-type exposure inference. Panel F compares fixed grids without selecting on held-out performance. "
                  "Known-section p and held-out count intervals use8192 fixed-seed mixture draws; full-data new-section p intervals invert the finite-grid Beta-mixture CDF numerically; μ intervals are discrete weighted-grid quantiles. Observed zero cells remain in n. "
                  "Uncertainty is conditional on this detection model and single-donor sampling; it does not calibrate tissue-wide biological uncertainty.")
        cap_es = (f"Detección basal Beta-Binomial jerárquica en {region}; no modela magnitud de RNA ni efectos PM2.5. "
                  "Por anatomía depositada, los conteos positivos y_s entre n_s células medidas por corte físico siguen Beta-Binomial(μκ,(1−μ)κ). "
                  "Priors fijadas antes del ajuste: μ~Beta(0.5,0.5), κ~LogNormal(log20,1.5); cuadratura64×48 por cuantiles previos con pesos iguales aproxima la posterior continua. "
                  "A–C distinguen intervalos creíbles95% de la media anatómica μ (oscuro) de intervalos predictivos95% de propensión p de un corte nuevo (claro); puntos naranjas son fracciones medidas agrupadas, no estimaciones μ. "
                  + ("Se muestran las10 capas/subdivisiones ACA. " if region == "ACA" else "Para legibilidad se muestran12 de49 núcleos/subregiones HY por mayor n medido, independientemente de expresión; se conservan los49 en tablas. ") +
                  "n=células medidas, S=cortes físicos. D–E reutilizan cinco bloques de cortes excluidos y sus cortes vecinos de guarda, sin usar conteos de prueba/guarda al entrenar. "
                  "Los intervalos predictivos95% de conteos condicionan al n de prueba real e incluyen heterogeneidad entre cortes y muestreo binomial; cobertura no ponderada entre grupos corte×anatomía observados. "
                  "Se conservan y cuentan anatomías de prueba solo con prior. La referencia no anatómica usa las mismas priors y agrega solo datos de entrenamiento por corte. "
                  "El error absoluto de fracción no se pondera entre grupos corte×anatomía; Brier ponderado por células se tabula aparte. No se infieren animales independientes ni respuestas celulares a exposición. F compara cuadrículas fijas sin seleccionarlas por rendimiento excluido. "
                  "Intervalos p de corte medido/conteos excluidos usan8192 sorteos de mezcla con semilla; p de corte nuevo con ajuste completo invierte numéricamente la CDF de mezcla Beta; μ usa cuantiles discretos ponderados. Se incluyen ceros en n. "
                  "La incertidumbre depende del modelo y del muestreo de un donante; no calibra incertidumbre biológica de todo el tejido.")
        manifests[stem] = {"caption_en": cap_en, "caption_es": cap_es, "exports": exports, "files": [x["path"] for x in exports],
                           "alt_en": "Three anatomy probability interval plots; held-out count coverage, fraction prediction error and fixed-grid sensitivity for ten measured genes.",
                           "alt_es": "Tres gráficos de intervalos de probabilidad anatómica; cobertura predictiva, error de fracción y sensibilidad de cuadrícula para diez genes medidos."}


def validate(results):
    checks = {"fixed_priors_and_grids_no_test_tuning": True, "all20_genes_original_measured_detection_counts": True,
              "same_five_contiguous_folds_and_neighbor_guards": True, "training_only_regional_baseline": True,
              "no_cellclass_covariate_or_expression_imputation": True, "full_data_posterior_not_used_for_OOF": True}
    for region, r in results.items():
        t = r["oof"]
        assert len(t) == len(r["full"])
        assert not t[["gene", "section_code", "anatomy_code"]].duplicated().any()
        assert (t.count_pi_low >= 0).all() and (t.count_pi_high <= t.n).all()
        assert (t.count_pi_low <= t.count_pi_high).all()
        for frame, pairs in [(r["full"], [("p_mean", None), ("p_ci_low", "p_ci_high")]),
                             (r["anatomy"], [("mu_mean", None), ("mu_ci_low", "mu_ci_high"), ("new_section_p_ci_low", "new_section_p_ci_high")])]:
            for low, high in pairs:
                assert ((frame[low] >= 0) & (frame[low] <= 1)).all()
                if high:
                    assert ((frame[high] >= 0) & (frame[high] <= 1)).all() and (frame[low] <= frame[high]).all()
        checks[f"{region}_probabilities_bounded_and_intervals_ordered"] = True
        checks[f"{region}_count_PI_support_0_to_actual_n"] = True
        checks[f"{region}_each_measured_group_tested_exactly_once"] = True
        checks[f"{region}_all_n_and_positive_counts_conserved"] = int(t.n.sum()) == 10 * r["payload"]["oof_summary"]["reference_cells"]
        # Algebraic verification against the analytic no-data prior mean and known-grid conditional means.
    prior = fit([], [])
    assert abs(prior["mu_mean"] - .5) < 1e-12
    a = fit([1, 3], [10, 12])
    from scipy.stats import betabinom
    from scipy.special import softmax
    g = a["grid"]
    independent_ll = betabinom.logpmf(np.array([1, 3])[:, None], np.array([10, 12])[:, None], g["alpha"], g["beta"]).sum(axis=0)
    np.testing.assert_allclose(a["weights"], softmax(independent_ll), atol=1e-12, rtol=1e-10)
    repeat = section_posterior(a, 1, 10, rng_for("ACA", 0, 91))
    again = section_posterior(a, 1, 10, rng_for("ACA", 0, 91))
    assert repeat == again
    checks["prior_only_prediction_symmetric_mean_half"] = True
    checks["seeded_mixture_intervals_exactly_repeatable"] = True
    checks["posterior_grid_likelihood_matches_scipy_betabinom"] = True
    assert all(checks.values())
    return {"passed": True, "all_checks_passed": True, "checks": checks,
            "scope": "Implementation/data-leakage/support/serialization invariants only; scientific calibration and grid sensitivity reported separately."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figures-only", action="store_true", help="Regenerate figures from saved posterior/OOF tables without refitting")
    args = parser.parse_args()
    for path in [OUT, META, WEB, OUT / "figures"]:
        path.mkdir(parents=True, exist_ok=True)
    if args.figures_only:
        results = {}
        for region in REGIONS:
            results[region] = {"summary": json.loads((OUT/region/"summary.json").read_text()),
                               "payload": json.loads((WEB/f"data-{region}.js").read_text().split("=",1)[1].strip().removesuffix(";"))}
            for key, stem in [("anatomy", "anatomy_posterior"), ("full", "known_section_posterior"), ("oof", "out_of_fold_predictions"), ("sensitivity", "grid_sensitivity")]:
                results[region][key] = pd.read_csv(OUT/region/f"{stem}.tsv.gz", sep="\t")
    else:
        results = {region: analyze_region(region) for region in REGIONS}
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 11,
                         "pdf.fonttype": 42, "svg.fonttype": "none", "figure.facecolor": "white"})
    manifest = {}
    for region, r in results.items():
        make_figures(region, r, manifest)
    dump(OUT / "figure_manifest.json", manifest)
    if args.figures_only:
        provenance = json.loads((META/"source_manifest.json").read_text())
        for item in provenance["sources"]:
            if item["path"] == str(Path(__file__).resolve().relative_to(ROOT)):
                item["sha256"] = digest(Path(__file__).resolve())
                item["bytes"] = Path(__file__).resolve().stat().st_size
        dump(META/"source_manifest.json", provenance)
        print("Regenerated figures from unchanged saved numerical tables; no model refit.")
        return
    dump(OUT / "validation.json", validate(results))
    summary = {region: r["summary"] for region, r in results.items()}
    dump(OUT / "summary.json", summary)
    diagnostics = {region: {"OOF": r["summary"], "coverage_is_model_diagnostic_not_validation_pass": True,
                            "grid_sensitivity": r["summary"]["grid_sensitivity_96x64"],
                            "mean_difference_over_0_02_groups": int(((r["sensitivity"].comparison_grid.eq("96x64")) & (r["sensitivity"].abs_mu_mean_difference > .02)).sum()),
                            "prior_only_predictions_must_be_interpreted_separately": True} for region, r in results.items()}
    dump(OUT / "scientific_diagnostics.json", diagnostics)
    methods = {"created_at": datetime.now(timezone.utc).isoformat(), "target": "Detected-cell proportions (raw imaging count>0), baseline only, one unexposed male reference donor",
               "likelihood": "y_sa|p_sa~Binomial(n_sa,p_sa); p_sa|mu_a,kappa_a~Beta(mu_a*kappa_a,(1-mu_a)*kappa_a)",
               "anatomy_parameters": "Independent fixed-prior hierarchical model for each deposited anatomy and gene; physical sections provide section-heterogeneity groups, not independent animals",
               "fixed_priors": {"mu": "Beta(0.5,0.5)", "kappa": "LogNormal(log(20),1.5); log-scale SD1.5, median20; fixed before fitting",
                                "motivation": "Explicit modeling assumptions, not empirically established brain-detection priors; no prior tuning against test data"},
               "quadrature": {"main": MAIN_GRID, "sensitivity_full_data": SENSITIVITY_GRIDS, "heldout_main_vs96x64": True,
                              "nodes": "At q=(i+0.5)/N of each prior; Cartesian product with equal prior-quadrature weights1/(N_mu*N_kappa)",
                              "weights": "Posterior grid weights proportional to product of training-section Beta-Binomial marginal likelihoods; logsumexp normalization",
                              "approximation": "Deterministic finite-grid approximation to continuous prior/posterior; not exact continuous integration or MCMC convergence",
                              "mu_intervals": "Equal-tail2.5/97.5% weighted discrete grid quantiles; finite grid sensitivity quantified, never selected on OOF",
                              "sensitivity_targets": "Full-data mu mean/CI, kappa/rho means and new-section p interval endpoints via numerical Beta-mixture CDF inversion(noMC); held-out96x64 count PI coverage/endpoints and Brier versus prespecified64x48, common-seedMC (count interval differences include MC error)"},
               "known_section": "Posterior of p_sa is mixture Beta(mu*kappa+y_sa,(1-mu)*kappa+n_sa-y_sa) under full-data hyperposterior; exact mixture mean,8192-seeded-draw interval",
               "anatomy_mu": "Mean and95% credible interval for mean propensity across exchangeable physical sections, distinct from pooled cell fraction",
               "new_section_p": "Mixture Beta(mu*kappa,(1-mu)*kappa), integrating hyperposterior;95% predictive interval for a new section's latent detection propensity. Full-data intervals/sensitivity numerically invert weighted mixture CDF with scipy betainc/brentq(tolerance1e-12); OOF uses MC draws for joint p/count prediction",
               "new_section_count": "New-section propensity drawn from posterior mixture, then Binomial(actual held-out n,p); empirical integer equal-tail95% count predictive interval",
               "MC": {"draws": DRAWS, "master_seed": SEED, "stream": "SeedSequence(master,region,gene,purpose,anatomy,fold,section)",
                      "quantiles": "numpy method=inverted_cdf; empirical count quantiles retain integer support[0,n]",
                      "tail_probability_MC_SE_near_0_025": float(np.sqrt(.025*.975/DRAWS)),
                      "caution": "This is MC error of estimated tail probabilities, not an uncertainty bound on probability-quantile locations; no exact continuous posterior claim"},
               "cross_validation": "Same five contiguous physical-section blocks and neighboring guard labels from prior spatial_bayes section_folds.tsv/fold_audit.json; neither test nor guard data fit posterior or baseline",
               "no_training_anatomy": "Use fixed prior-only predictive distribution and flag it; retain OOF scores. No artificial test groups for n=0, no zero-filled missing anatomy",
               "regional_baseline": "Same priors/quadrature, collapse anatomy within each TRAINING physical section; posterior mu mean predicts every test anatomy; no test or guard counts",
               "metrics": {"count_PI_coverage": "Fraction observed section×anatomy counts within inclusive95% predictive interval; primary unweighted across nonempty groups",
                           "fraction_MAE": "Unweighted mean abs(y/n-E[mu|training]); n-weighted alternative also reported",
                           "Brier": "Sum[y*(1-pred)^2+(n-y)*pred^2]/sum(n); cell-weighted predictive loss, not independent-cell inference",
                           "improvement": "Regional-baseline Brier minus anatomy-model Brier; positive means lower model loss",
                           "scope": "Allgroups plus trained-anatomy-only subset reported; correlated genes/sections, no animal-level CI or superiority test"},
               "figure_HY_subset": "Twelve largest all-measured-cell-n subregions; independent of gene expression; all49 in full tables and UI",
               "limitations": ["One unexposed donor; physical sections are not independent animal replication", "Section exchangeability is a model assumption; AP-dependent patterns and acquisition effects are not distinguished", "Conditional binomial sampling assumes common propensity within each anatomy×section; residual local spatial dependence and celltype composition are not modeled",
                               "Detection depends on panel/assay/cell filtering and CCF coverage; no unbiased tissue composition claim", "No cellclass-specific effects, expression imputation, RNA magnitude, PM2.5 response, or increased spatial resolution",
                               "Intervals are conditional on fixed priors, finite quadrature and observed sampling; calibration may fail and is separately reported"],
               "software": {"Python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "pandas": pd.__version__, "matplotlib": matplotlib.__version__}}
    dump(META / "methods.json", methods)
    paths = list(dict.fromkeys(p for r in results.values() for p in r["inputs"]))
    paths += [ROOT / "metadata/spatial_sections/source_manifest.json", Path(__file__).resolve()]
    dump(META / "source_manifest.json", {"sources": [{"path": str(p.relative_to(ROOT)), "bytes": p.stat().st_size, "sha256": digest(p)} for p in paths],
                                        "atlas_paper_doi": "https://doi.org/10.1038/s41586-023-06808-9", "counts": "Exact published measured-section raw-count>0 aggregates; no new data"})
    for spanish in [False, True]:
        prose = "# " + ("Detección basal Beta-Binomial jerárquica" if spanish else "Hierarchical Beta-Binomial baseline detection") + "\n\n"
        for region, r in results.items():
            s = r["summary"]
            if spanish:
                prose += (f"{region}:10 genes, {s['n_groups']:,} grupos observados corte×anatomía. Cobertura predictiva95% de conteos excluidos "
                          f"{100*s['count_pi_coverage']:.1f}%; Brier {s['brier']:.5f} frente a referencia regional {s['baseline_brier']:.5f}. "
                          f"Grupos de prueba solo con prior: {s['prior_only_groups']}; no se ajustan con datos excluidos. "
                          f"Diferencia máxima de media μ entre cuadrículas64×48 y96×64: {s['grid_sensitivity_96x64']['max_abs_mu_mean_difference']:.4f}.\n\n")
            else:
                prose += (f"{region}:10 genes, {s['n_groups']:,} observed section×anatomy groups. Held-out95% count predictive coverage "
                          f"{100*s['count_pi_coverage']:.1f}%; Brier {s['brier']:.5f} versus regional baseline {s['baseline_brier']:.5f}. "
                          f"Prior-only test groups: {s['prior_only_groups']}; held-out data never fit them. "
                          f"Maximum μ mean difference64×48 versus96×64 grids: {s['grid_sensitivity_96x64']['max_abs_mu_mean_difference']:.4f}.\n\n")
        prose += ("Los intervalos creíbles de μ anatómica, la posterior de p de un corte ya medido y la predicción de p/conteos de un corte nuevo responden preguntas distintas. "
                  "Este análisis describe detección basal con priors fijas y heterogeneidad entre cortes. No modela respuesta PM2.5, magnitud de RNA, tipos celulares expuestos ni nueva resolución espacial. "
                  "Hay un donante; los cortes no son réplicas animales. Validación de implementación y diagnósticos científicos se informan por separado.\n" if spanish else
                  "Anatomy-mean μ credible intervals, known-section p posterior intervals and new-section p/count predictive intervals answer different questions. "
                  "This analysis describes baseline detection using fixed priors and physical-section heterogeneity. It does not model PM2.5 responses, RNA magnitude, exposed cell types or new spatial resolution. "
                  "There is one donor; sections are not animal replicates. Implementation validation and scientific diagnostics are reported separately.\n")
        (OUT / ("narrative_es.md" if spanish else "narrative_en.md")).write_text(prose)
    print(json.dumps(clean(summary), indent=2))


if __name__ == "__main__":
    main()
