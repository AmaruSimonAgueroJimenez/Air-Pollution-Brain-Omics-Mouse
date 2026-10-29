#!/usr/bin/env python3
"""Validate bilingual presentation structure and provenance, apart from visual QA."""
from pathlib import Path
import base64
import hashlib
import json
import re
from urllib.parse import unquote, urlsplit, parse_qs
import zipfile

from bs4 import BeautifulSoup

ROOT=Path(__file__).resolve().parents[1]
BASE="https://amarusimonaguerojimenez.github.io/Air-Pollution-Brain-Omics-Mouse/"
OUT=ROOT/"results/defense_slides"

def main():
    errors=[]; checks={}
    manifest=json.loads((OUT/"presentation_manifest.json").read_text())
    expected=["title-slide",*manifest["explicit_ids"]]
    for lang,relative in (("en","defense.html"),("es","es/defense.html")):
        p=ROOT/"docs"/relative
        soup=BeautifulSoup(p.read_text(),"html.parser")
        slides=soup.select(".reveal .slides > section")
        ids=[s.get("id") for s in slides]
        if ids!=expected: errors.append(f"{lang}: slide IDs/order mismatch")
        notes=[s.select_one("aside.notes") for s in slides[1:]]
        if any(n is None or len(n.get_text(" ",strip=True))<100 for n in notes):
            errors.append(f"{lang}: absent/trivial speaker notes")
        math=soup.find_all("math")
        qmd=p.with_suffix(".qmd").read_text()
        expected_math=len(re.findall(r"\$\$[\s\S]*?\$\$|(?<!\$)\$(?!\$)[^\n$]*?\$(?!\$)",qmd))
        if len(math)!=expected_math or expected_math<1:
            errors.append(f"{lang}: rendered equations do not match the source")
        images=soup.select(".reveal .slides img")
        if any(not (i.get("src") or i.get("data-src","")).startswith("data:") for i in images):
            errors.append(f"{lang}: figure requires an external asset")
        figure_specs=[s for s in manifest["slides"] if s.get("figure")]
        if len(images)!=len(figure_specs):
            errors.append(f"{lang}: figure count differs from the slide specification")
        for element,spec in zip(images,figure_specs):
            uri=element.get("src") or element.get("data-src", "")
            stem=spec["figure"]+("_es" if lang=="es" else "")
            original=ROOT/"results/defense_slides/figures"/(stem+".png")
            if not uri.startswith("data:image/png;base64,") or hashlib.sha256(base64.b64decode(uri.split(",",1)[1])).digest()!=hashlib.sha256(original.read_bytes()).digest():
                errors.append(f"{lang}: embedded scientific figure differs from saved output: {stem}")
        explanations=[s for s in slides if (s.get("id") or "").startswith("explain-")]
        if len(explanations)!=12:
            errors.append(f"{lang}: expected twelve visible explanation slides")
        for slide in explanations:
            visible=" ".join(p.get_text(" ",strip=True) for p in slide.find_all("p",recursive=False))
            if len(visible.split())<45:
                errors.append(f"{lang}: explanation lacks visible detail: {slide['id']}")
        for sid,region in (("hy-3d","HY"),("aca-3d","ACA")):
            slide=soup.find("section",id=sid)
            frame=slide.select_one("iframe[data-map-src]") if slide else None
            if frame is None:
                errors.append(f"{lang}: missing in-slide interactive map: {sid}")
                continue
            component=urlsplit(frame["data-map-src"])
            query=parse_qs(component.query)
            local_component=(p.parent/unquote(component.path)).resolve()
            if local_component!=ROOT/"docs/assets/defense-slides/inline-map.html" or query.get("region")!=[region] or query.get("lang")!=[lang] or frame.get("src"):
                errors.append(f"{lang}: invalid lazy map source or language/region: {sid}")
            if not slide.select_one(".map-poster img") or not slide.select_one(".live-map-caption") or not slide.select_one(".map-enlarge"):
                errors.append(f"{lang}: map fallback, explanation or enlargement absent: {sid}")
        for anchor in soup.select("a[data-map-slide]"):
            if anchor["data-map-slide"] not in ids:
                errors.append(f"{lang}: broken map slide jump")
        resources=[s.get("src","") for s in soup.find_all("script",src=True)]
        resources += [s.get("href","") for s in soup.select("link[rel=stylesheet]")]
        external=[u for u in resources if not u.startswith("data:")]
        if external: errors.append(f"{lang}: non-embedded runtime/CSS: {external}")
        if "defense-notes-toggle" not in p.read_text():
            errors.append(f"{lang}: offline notes runtime missing")
        citations=soup.select(".reveal .source a[href]")
        if any(".json" in a.get("href", "").lower() or ".json" in a.get_text().lower() for a in citations):
            errors.append(f"{lang}: internal JSON presented as a scientific citation")
        scholarly=[a for a in citations if not a["href"].startswith(BASE)]
        if len(scholarly)<80:
            errors.append(f"{lang}: missing scientific references")
        source_links=0
        for a in soup.select("a[href]"):
            href=a["href"]
            if href.startswith(BASE):
                target=ROOT/"docs"/unquote(urlsplit(href[len(BASE):]).path)
                if not target.exists(): errors.append(f"{lang}: unavailable published source target {href}")
                source_links+=1
        text=" ".join(s.get_text(" ",strip=True) for s in slides)
        longwords=sorted(set(re.findall(r"[A-Za-zÁÉÍÓÚáéíóúñÑ]{30,}",text)))
        if longwords: errors.append(f"{lang}: possible compacted prose {longwords}")
        if re.search(r"and three\?|por confirmar",text,re.I) or re.search(r"\b(?:TODO|FIXME)\b",text):
            errors.append(f"{lang}: unresolved draft wording")
        for sid,tokens in {
            "cortex-reproduction":["109","1,205","1,222","639"],
            "hy-primary":["11"],
            "beta-performance":["0.13725","0.14051","0.09246","0.03621"],
            "numerics":["0.000192","0.000537","0.0203","0.0181"],
            "scenario-method":["0.42808","1.99600"],
        }.items():
            found=soup.find("section",id=sid).get_text(" ",strip=True)
            if not all(t in found for t in tokens): errors.append(f"{lang}: missing result anchors in {sid}")
        checks[lang]={"slides":len(slides),"notes":len(notes),"native_mathml_expressions":len(math),
                      "embedded_figures":len(images),"embedded_runtime":not external,
                      "visible_explanation_slides":len(explanations),"in_slide_interactive_gene_maps":2,
                      "published_evidence_links":source_links,"scholarly_citations":len(scholarly),"bytes":p.stat().st_size,
                      "sha256":hashlib.sha256(p.read_bytes()).hexdigest()}
    for item in manifest["source_files"]:
        if hashlib.sha256((ROOT/item["path"]).read_bytes()).hexdigest()!=item["sha256"]:
            errors.append("Changed scientific source: "+item["path"])
    for path,expected_hash in manifest["presentation_code_sha256"].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=expected_hash:
            errors.append("Changed presentation builder: "+path)
    archive=ROOT/"docs/downloads/defense-slides-bilingual.zip"
    if archive.exists():
        with zipfile.ZipFile(archive) as z:
            for relative in ("defense.html","es/defense.html"):
                if z.read(relative)!=(ROOT/"docs"/relative).read_bytes():
                    errors.append("Stale portable archive: "+relative)
            required_assets=["assets/defense-slides/"+n for n in ("inline-map.html","inline-map.js","inline-map.css")]
            required_assets += ["assets/viewer/"+n for n in ("plotly.min.js","annotations.js","data.js")]
            required_assets += ["assets/viewer-cortex/data.js"]
            for directory in ("viewer","viewer-cortex"):
                required_assets += [str(p.relative_to(ROOT/"docs")) for p in sorted((ROOT/"docs/assets"/directory/"genes").glob("*.js"))]
            for relative in required_assets:
                if relative not in z.namelist() or z.read(relative)!=(ROOT/"docs"/relative).read_bytes():
                    errors.append("Missing/stale portable map asset: "+relative)
            checks["portable_interactive_map_assets"]=len(required_assets)
        checks["portable_zip_bytes"]=archive.stat().st_size
    else:
        errors.append("Missing portable bilingual archive")
    report={"scope":"Structure, embedded resources, paired IDs, result anchors and frozen-source hashes. Visual/browser review is separate; no model validity certificate.",
            "checks":checks,"errors":errors,"passed":not errors}
    (OUT/"validation.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))
    raise SystemExit(bool(errors))

if __name__=="__main__": main()
