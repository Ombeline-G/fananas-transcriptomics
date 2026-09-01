#!/usr/bin/env python3
"""FC share by cerebellar lobule / territory in MERFISH (review point 18).

Every censused Zhuang soma (190 FC + 3,570 confident BG) is assigned to the
Allen CCF 2020 parcellation through the ABC-atlas per-cell CCF registration
(cell_label -> parcellation_index -> structure acronym). Structure-level
acronyms map vermal lobules I-X and the hemispheric territories:

  LING I · CENT II-III · CUL IV-V · DEC VI · FOTU VII · PYR VIII · UVU IX ·
  NOD X · hemispheric: SIM, AN (crus I/II), PRM, COPY, FL, PFL.

Shares are FC/(FC + confident BG) per structure. Cells are not independent
within animals (Statistics): shares are descriptive, reported with per-animal
counts and no cell-level intervals; the dense-sampling animal (ABCA-3)
carries the profile and the two sparse animals are shown as counts.

Output: tables/fc_share_by_lobule.csv (per structure x animal + 'all' rows,
plus aggregate territories).
"""
import pandas as pd

from project_paths import ROOT, ZHUANG

TABLES = ROOT / "comparison" / "tables"
CCF_META = ZHUANG / "abc_cache" / "metadata"

LOBULE = {"LING": "I", "CENT": "II-III", "CUL": "IV-V", "DEC": "VI", "FOTU": "VII",
          "PYR": "VIII", "UVU": "IX", "NOD": "X",
          "SIM": "hemis.", "AN": "hemis.", "PRM": "hemis.", "COPY": "hemis.",
          "FL": "flocc.", "PFL": "flocc."}
TERRITORY = {"LING": "anterior vermis (I-V)", "CENT": "anterior vermis (I-V)",
             "CUL": "anterior vermis (I-V)",
             "DEC": "posterior vermis (VI-IX)", "FOTU": "posterior vermis (VI-IX)",
             "PYR": "posterior vermis (VI-IX)", "UVU": "posterior vermis (VI-IX)",
             "NOD": "nodulus (X)",
             "SIM": "hemispheres", "AN": "hemispheres", "PRM": "hemispheres",
             "COPY": "hemispheres",
             "FL": "flocculus/paraflocculus", "PFL": "flocculus/paraflocculus"}
ORDER = ["LING", "CENT", "CUL", "DEC", "FOTU", "PYR", "UVU", "NOD",
         "SIM", "AN", "PRM", "COPY", "FL", "PFL"]

meta = pd.read_parquet(TABLES / "zhuang_astro_cells.parquet")[["cell", "group", "animal", "section"]]
meta = meta[meta.group.isin(["FC", "BG"])]
assert meta.group.value_counts().to_dict() == {"BG": 6084, "FC": 362}   # 4 animals since ABCA-1 integration (2026-08-18)

membership = pd.read_csv(CCF_META / "Allen-CCF-2020" / "20230630" / "views"
                         / "parcellation_to_parcellation_term_membership_acronym.csv")
ccf = pd.concat([
    pd.read_csv(CCF_META / f"{animal}-CCF" / "20230830" / "ccf_coordinates.csv",
                dtype={"cell_label": str})[["cell_label", "parcellation_index"]]
    for animal in ["Zhuang-ABCA-1", "Zhuang-ABCA-2", "Zhuang-ABCA-3", "Zhuang-ABCA-4"]
]).set_index("cell_label")

j = (meta.set_index(meta.cell.astype(str)).join(ccf, how="left")
     .merge(membership, on="parcellation_index", how="left"))
assert j.parcellation_index.notna().all(), "unmatched cells in CCF registration"
assert (j.division == "CB").all() or set(j[j.division != "CB"].structure) <= {"CB-unassigned"}
unassigned = j[~j.structure.isin(ORDER)]
print(f"CB-unassigned (excluded from shares): {unassigned.group.value_counts().to_dict()}")
j = j[j.structure.isin(ORDER)]

rows = []
def add_rows(frame, animal_label):
    for structure, sub in frame.groupby("structure"):
        n_fc = int((sub.group == "FC").sum())
        n_bg = int((sub.group == "BG").sum())
        rows.append(dict(level="structure", unit=structure, lobule=LOBULE[structure],
                         territory=TERRITORY[structure], animal=animal_label,
                         n_FC=n_fc, n_BG=n_bg,
                         share_pct=round(100 * n_fc / (n_fc + n_bg), 2) if n_fc + n_bg else float("nan")))
    for territory, sub in frame.groupby(frame.structure.map(TERRITORY)):
        n_fc = int((sub.group == "FC").sum())
        n_bg = int((sub.group == "BG").sum())
        rows.append(dict(level="territory", unit=territory, lobule="", territory=territory,
                         animal=animal_label, n_FC=n_fc, n_BG=n_bg,
                         share_pct=round(100 * n_fc / (n_fc + n_bg), 2) if n_fc + n_bg else float("nan")))

add_rows(j, "all")
for animal, sub in j.groupby("animal"):
    add_rows(sub, animal)

out = pd.DataFrame(rows)
out.to_csv(TABLES / "fc_share_by_lobule.csv", index=False)

view = out[(out.level == "structure") & (out.animal == "all")].set_index("unit").loc[
    [s for s in ORDER]]
print(view[["lobule", "territory", "n_FC", "n_BG", "share_pct"]].to_string())
print()
terr = out[(out.level == "territory") & (out.animal == "all")]
print(terr[["unit", "n_FC", "n_BG", "share_pct"]].to_string(index=False))
print()
for animal in ["Zhuang-ABCA-2", "Zhuang-ABCA-4"]:   # sparse animals only; ABCA-1/3 shown in the main view
    sub = out[(out.level == "structure") & (out.animal == animal) & (out.n_FC > 0)]
    print(f"{animal} FC locations:", dict(zip(sub.unit, sub.n_FC)))
