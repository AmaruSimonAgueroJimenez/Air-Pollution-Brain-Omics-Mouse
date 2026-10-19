#!/usr/bin/env python3
"""Check the compiled local report and required completed-study artifacts.

This complements actual browser inspection; it cannot certify visual appearance.
"""
from __future__ import annotations

import argparse
import ast
from datetime import datetime, timezone
from html.parser import HTMLParser
import json
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "docs/_site"


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
        self.ids = set()

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        if "id" in attributes:
            self.ids.add(attributes["id"])
        for name in ("href", "src"):
            if attributes.get(name):
                self.links.append(attributes[name])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--allow-incomplete-ozone", action="store_true")
    args = parser.parse_args()
    errors, checks = [], {}
    pages = sorted((ROOT / "docs").glob("*.qmd")) + sorted((ROOT / "docs/es").glob("*.qmd"))
    for page in pages:
        for index, chunk in enumerate(page.read_text().split("```{python}")[1:], 1):
            try:
                ast.parse(chunk.split("```")[0], filename=f"{page.name}:chunk{index}")
            except SyntaxError as error:
                errors.append(str(error))
        if not (SITE / page.relative_to(ROOT / "docs").with_suffix(".html")).exists():
            errors.append("Missing rendered page: " + str(page.relative_to(ROOT / "docs")))
    checks["qmd_pages"] = len(pages)
    documents = {}
    widget_dependencies = []
    for path in SITE.rglob("*.html"):
        # Downloaded QC reports are provenance, not authored report pages.
        if "downloads" in path.relative_to(SITE).parts:
            continue
        source = path.read_text()
        if "requirejs@" in source or "jquery@" in source or "define('jquery'" in source:
            widget_dependencies.append(str(path.relative_to(SITE)))
        parsed = Links()
        parsed.feed(source)
        documents[path.resolve()] = parsed
    link_count = 0
    for path, parsed in documents.items():
        for link in parsed.links:
            parts = urlsplit(link)
            if parts.scheme or parts.netloc or link.startswith("data:"):
                continue
            target = (SITE / unquote(parts.path).lstrip("/")) if parts.path.startswith("/") else path.parent / unquote(parts.path)
            if not parts.path:
                target = path
            target = target.resolve()
            link_count += 1
            if not target.exists():
                errors.append(f"Missing target in {path.relative_to(SITE)}: {link}")
            elif parts.fragment and target in documents and unquote(parts.fragment) not in documents[target].ids:
                # Empty href="#" is a legitimate toggle; nonempty anchors must resolve.
                errors.append(f"Missing anchor in {path.relative_to(SITE)}: {link}")
    checks.update(html_documents=len(documents), local_links_checked=link_count)
    checks["unexpected_static_widget_dependencies"] = widget_dependencies
    if widget_dependencies:
        errors.append("Unexpected widget loader in static report pages: " + ", ".join(widget_dependencies))
    expected = [
        "results/hypothalamus/validation_checks.json", "results/cortex/validation.json",
        "results/mirna/validation.json", "results/spatial_viewer/validation_checks.json",
        "results/hypothalamus/independent_deseq2_check/summary.json",
        "results/spatial_reference/validation_checks.json",
        "results/enrichment/validation_checks.json",
    ]
    for relative in expected:
        if not (ROOT / relative).exists():
            errors.append("Missing analysis audit: " + relative)
    for relative in (
        "results/hypothalamus/validation_checks.json",
        "results/spatial_reference/validation_checks.json",
        "results/enrichment/validation_checks.json",
    ):
        path = ROOT / relative
        if path.exists() and not json.loads(path.read_text()).get("all_checks_passed", False):
            errors.append("Analysis audit reports unsuccessful checks: " + relative)
    r_check = ROOT / "results/hypothalamus/independent_deseq2_check/summary.json"
    if r_check.exists():
        implementation_check = json.loads(r_check.read_text())
        if not all(implementation_check.get(k, False) for k in ("primary_analysis_unchanged", "same_gene_universe", "all_descriptive_tolerances_met")):
            errors.append("Independent implementation comparison has an unresolved check")
    ozone_file = ROOT / "results/ozone/summary.json"
    ozone = json.loads(ozone_file.read_text()) if ozone_file.exists() else {}
    ozone_complete = ozone.get("complete", False)
    checks["ozone_inference_complete"] = ozone_complete
    checks["ozone_quantification_checkpoints"] = len(list((ROOT / "results/ozone/quant").glob("*/ozone_quant_complete.json")))
    if not args.allow_incomplete_ozone:
        if not ozone_complete or checks["ozone_quantification_checkpoints"] != 24:
            errors.append("Ozone requires 24 verified libraries and completed gene-level inference")
        validation_path = ROOT / "results/ozone/validation.json"
        if not validation_path.exists() or not json.loads(validation_path.read_text()).get("passed", False):
            errors.append("Ozone numerical validation has not passed")
        enrichment = json.loads((ROOT / "results/enrichment/summary.json").read_text())["datasets"]
        for key in ("M3_E3", "M3_E4", "M17_E3", "M17_E4"):
            if not (ROOT / f"results/ozone/de_{key}_O3_vs_FA.csv").exists():
                errors.append("Missing ozone contrast: " + key)
            if "ozone_" + key not in enrichment:
                errors.append("Missing ozone enrichment: " + key)
        for key in ("hypothalamus", "cortex", "ozone_M3_E3", "ozone_M3_E4", "ozone_M17_E3", "ozone_M17_E4"):
            for name in ("GSEA_all_terms.tsv", "GSEA_ranked_genes.tsv", "ORA_up_all_terms.tsv", "ORA_down_all_terms.tsv"):
                if not (ROOT / "results/enrichment" / key / name).exists():
                    errors.append(f"Missing enrichment artifact: {key}/{name}")
    result = {
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "scope": "Compiled local targets, source syntax and required artifacts; visual review is separate",
        "allows_incomplete_ozone": args.allow_incomplete_ozone,
        "checks": checks, "errors": errors, "passed_within_checked_scope": not errors,
    }
    destination = ROOT / "results/report_validation/compiled_report_checks.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
