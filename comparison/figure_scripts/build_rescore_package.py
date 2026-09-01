#!/usr/bin/env python3
"""Build the SECOND-RATER swipe package (review Ombeline 3, point 13).

The original swipe tools iterate over ALL candidates of their dataset (240 +
406 + 361 + 145), with no subset mechanism: asking a second rater to locate
the 173 sampled candidates by hand is not workable. This script produces
comparison/rescore_package/ — four self-contained mini-tools, one per decision
file, each restricted to the sampled candidates in the sample's shuffled
presentation order, with:
  * a SEPARATE localStorage key (never collides with the original tools);
  * the export renamed rater2_<key>.csv and its dataset column forced to the
    sample key, i.e. exactly the format compute_rescore_kappa.py expects;
  * only the map sections actually referenced by the sampled candidates.
The rater sees the same geometry-only rendering as the original review; no
decision or expression value is present anywhere in the package.
"""
import json
import re
import shutil

import pandas as pd

from project_paths import ROOT, ZHUANG, HAO, SUYAMA

TABLES = ROOT / "comparison" / "tables"
OUT = ROOT / "comparison" / "rescore_package"
OUT.mkdir(exist_ok=True)

TOOLS = {
    "zhuang": (ZHUANG / "swipe_validation" / "swipe_zhuang.json",
               ZHUANG / "swipe_validation" / "swipe.html"),
    "zhuang_abca1": (ZHUANG / "swipe_validation" / "swipe_abca1.json",
                     ZHUANG / "swipe_validation" / "swipe_abca1.html"),
    "hao": (HAO / "swipe_validation" / "swipe_data.json",
            HAO / "swipe_validation" / "swipe.html"),
    "suyama": (None,  # data inline in the HTML
               SUYAMA / "results" / "swipe_validation" / "swipe.html"),
}

sample = pd.read_csv(TABLES / "swipe_rescore_sample.csv", dtype={"cell_label": str})

for key, (json_path, html_path) in TOOLS.items():
    want = sample.loc[sample.dataset == key, "cell_label"].tolist()
    html = html_path.read_text()

    if json_path is not None:
        data = json.loads(json_path.read_text())
    else:  # suyama: extract the inline window.SWIPE_DATA={...}
        m = re.search(r"window\.SWIPE_DATA=(\{.*\})\s*(?:;)?\s*</script>", html)
        assert m, "inline SWIPE_DATA not found in the Suyama tool"
        data = json.loads(m.group(1))

    by_id = {str(c["id"]): c for c in data["cards"]}
    missing = [w for w in want if w not in by_id]
    assert not missing, f"{key}: {len(missing)} sampled ids absent from the tool data"
    cards = [by_id[w] for w in want]                      # sample (shuffled) order
    secs = {c["section"] for c in cards}
    out = {"meta": {**data["meta"],
                    "dataset": f"SECOND RATER · {key} (rescore)",
                    "n": len(cards),
                    "source": data["meta"].get("source", "") + " | seed-fixed 15% rescore sample"},
           "sections": {s: data["sections"][s] for s in sorted(secs)},
           "cards": cards}
    (OUT / f"rescore_{key}.js").write_text("window.SWIPE_DATA=" + json.dumps(out))

    # ---- patch the HTML copy ----
    # 1. data source -> our filtered js
    if json_path is not None:
        html, n = re.subn(r'<script src="[^"]*"></script>',
                          f'<script src="rescore_{key}.js"></script>', html, count=1)
        assert n == 1, f"{key}: data <script src> not found"
    else:
        html, n = re.subn(r"<script>window\.SWIPE_DATA=\{.*?</script>",
                          f'<script src="rescore_{key}.js"></script>', html, count=1)
        assert n == 1, "suyama: inline data block not replaced"
    # 2. separate localStorage key
    html, n = re.subn(r'const KEY="[^"]+"', f'const KEY="rescore_v1_{key}"', html, count=1)
    assert n == 1, f"{key}: localStorage KEY not patched"
    # 3. export: filename + dataset column forced to the sample key
    html, n = re.subn(r'a\.download="[^"]+\.csv"', f'a.download="rater2_{key}.csv"', html, count=1)
    assert n == 1, f"{key}: export filename not patched"
    html, n = re.subn(r"rows\.push\(\[c\.id,\s*c\.dataset,",
                      f'rows.push([c.id,"{key}",', html, count=1)
    assert n == 1, f"{key}: export dataset column not patched"
    (OUT / f"rescore_{key}.html").write_text(html)
    print(f"  {key}: {len(cards)} candidates, {len(secs)} sections -> rescore_{key}.html")

(OUT / "README.md").write_text(f"""# Second-rater swipe package ({len(sample)} candidates)

Open each `rescore_*.html` in a browser (double-click; no server needed) and score every
candidate — same rules as the original review: geometry only, judge whether the starred soma
is a clear, isolated ML candidate. Keys: y / → = yes · n / ← = no · m / ↑ = maybe.
**Resolve every "maybe" to yes or no before exporting** (the original protocol resolved all
final decisions to yes/no; the kappa script rejects anything else). Progress is saved locally
in the browser (separate storage from the original tools).

Then export each tool ("⬇ Exporter les décisions") and concatenate the four files:

    head -n1 rater2_zhuang.csv > rater2_all.csv
    for f in rater2_zhuang.csv rater2_zhuang_abca1.csv rater2_hao.csv rater2_suyama.csv; do
      tail -n +2 "$f" >> rater2_all.csv
    done
    python3 comparison/figure_scripts/compute_rescore_kappa.py rater2_all.csv

The script prints per-dataset and pooled Cohen's kappa + raw agreement, to be inserted in
Methods, "Manual validation". Do not consult the original decision files while scoring.
""")
print("package:", OUT)
