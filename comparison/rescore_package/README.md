# Second-rater swipe package (173 candidates)

Open each `rescore_*.html` in a browser (double-click; no server needed) and score every
candidate — same rules as the original review: geometry only, judge whether the starred soma
is a clear, isolated ML candidate. Keys: y / → = yes · n / ← = no · m / ↑ = maybe.
**Resolve every "maybe" to yes or no before exporting** (the original protocol resolved all
final decisions to yes/no; the kappa script rejects anything else). Progress is saved locally
in the browser (separate storage from the original tools).

If a tool opens with pre-filled decisions (same browser as a previous scoring session), click
« Effacer / réinitialiser » before starting — otherwise you inherit someone else's decisions.

Then export each tool ("⬇ Exporter les décisions") and concatenate the four files (awk, not
`tail`: browser exports carry no final newline, which a `tail` merge would silently corrupt):

    awk 'FNR==1 && NR!=1 {next} {print}' rater2_zhuang.csv rater2_zhuang_abca1.csv \
      rater2_hao.csv rater2_suyama.csv > rater2_all.csv
    python3 comparison/figure_scripts/compute_rescore_kappa.py rater2_all.csv

The script prints per-dataset and pooled Cohen's kappa + raw agreement, to be inserted in
Methods, "Manual validation". Do not consult the original decision files while scoring.

Status 2026-08-31: scoring COMPLETE. All four rater-2 batches are A.S.'s (zhuang scored
2026-08-20; the other three 2026-08-31). `rater2_all.csv` is the official merge (zhuang = first,
blind pass); `rater2_zhuang_retest.csv` is A.S.'s unplanned second pass on zhuang, kept for the
record but excluded from the official κ (pooled κ = 0.463 — in Methods "Manual validation" and
addendum 54). The tools remain usable for full re-review.
