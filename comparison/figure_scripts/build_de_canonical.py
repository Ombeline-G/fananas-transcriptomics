#!/usr/bin/env python3
"""Build canonical DE counts plus a transparent audit of every statistical hit.

The audit gene sets are defined below, in this file only, and exported gene by
gene with their provenance to TableS_DE_audit_gene_sets.csv (review point 10).
None of the sets is constructed from the FC-versus-BG comparison: they come
from cell-type marker literature/atlas taxonomies (Purkinje, granule/synaptic),
standard technical-artifact lists, the Stereo-seq selection pipeline (gate
genes), canonical pan-astrocyte markers overlapping the
harmonized astro module (astro lean), and the Bergmann signature derived from the
BG-versus-parenchymal contrast without FC (Bergmann program). The labels
classify hits; they do not remove them.
"""
import os

import numpy as np
import pandas as pd

from project_paths import ROOT as PROJECT_ROOT

TAB = os.path.join(str(PROJECT_ROOT), "comparison", "tables")

PURKINJE = {"Calb1", "Car8", "Pcp2", "Pcp4", "Itpr1", "Grid2", "Aldoc", "Slc1a6", "Homer3"}
GRANULE = {"Gabra6", "Gng13", "Neurod1", "Cbln1", "Cbln3", "Snap25", "Stmn2", "Syn1", "Sncb",
           "Slc17a7", "Fat2", "Rbfox3", "Adcy1", "Nrxn3"}
TECH = {"Malat1", "AY036118", "Bc1", "Neat1"}
HAO_GATES = {"Aqp4", "Gja1", "Glul", "Aldh1l1"}
LEAN = {"Aqp4", "Ntsr2", "Gfap", "Stk32a", "Gja1", "Agt", "Aldh1l1", "Fgfr3", "S100b",
        "Gjb6", "Slc1a2", "Slc1a3", "Ndrg2"}
BERG = {"Gdf10", "Hopx", "A2m", "Timp4", "Slc6a11", "Gpr37l1", "Tnc", "Fabp7"}

FAMILY_TABLE = [
    ("contam_Purkinje", PURKINJE,
     "Canonical Purkinje-cell markers (calbindin, carbonic anhydrase 8, Pcp genes, "
     "Purkinje-enriched signaling); cerebellar marker literature and atlas cluster markers.",
     "yes - cell-type markers, no FC/BG input"),
    ("contam_granule/synaptic", GRANULE,
     "Canonical granule-cell and pan-neuronal/synaptic markers; cerebellar marker "
     "literature and atlas cluster markers.",
     "yes - cell-type markers, no FC/BG input"),
    ("technical", TECH,
     "Nuclear-retained lncRNAs (Malat1, Neat1), rRNA-associated repeat AY036118, dendritic "
     "BC1 RNA; all mt- transcripts are additionally classified technical by rule.",
     "yes - standard artifact genes"),
    ("gate_gene(Hao)", HAO_GATES,
     "The four genes of the astrocyte bin-selection module of the Stereo-seq pipeline; "
     "flagged only in that dataset, where selection renders them uninterpretable "
     "(depth-matched analysis README).",
     "yes - selection-pipeline genes, fixed before DE"),
    ("astro_lean", LEAN,
     "Canonical pan-astrocyte markers (Gfap, Aqp4, Agt, Gja1, Gjb6, Aldh1l1, Fgfr3, S100b, "
     "Slc1a2, Slc1a3, Ndrg2) plus the astrocyte-enriched genes Ntsr2 and Stk32a; overlaps "
     "the harmonized astro module, defined before that comparison "
     "(HARMONIZATION_protocol.md).",
     "yes - pan-astrocyte markers; module defined before the comparison"),
    ("bergmann_program", BERG,
     "The Bergmann signature of this study, defined from the BG-versus-parenchymal "
     "contrast without any use of FC (Results section 3), plus literature Bergmann markers.",
     "yes - defined from BG vs parenchymal, FC never used"),
]


def classify(gene, dataset):
    if gene.startswith("mt-") or gene in TECH:
        return "technical"
    if gene in PURKINJE:
        return "contam_Purkinje"
    if gene in GRANULE:
        return "contam_granule/synaptic"
    if dataset == "hao" and gene in HAO_GATES:
        return "gate_gene(Hao)"
    if gene in LEAN:
        return "astro_lean"
    if gene in BERG:
        return "bergmann_program"
    return "other"


audit, canonical, counts = [], [], []
for dataset in ["zhuang", "hao", "suyama"]:
    df = pd.read_csv(os.path.join(TAB, f"fig4_DE_{dataset}_depthmatched.csv"))
    sig = df[(df.padj < 0.05) & (df.log2FC.abs() > 0.25)].copy()
    sig["dataset"] = dataset
    sig["class"] = [classify(g, dataset) for g in sig.gene]
    sig["direction"] = np.where(sig.log2FC > 0, "up_FC", "up_BG")
    # Hao depth-matched hits are retained in the audit but are not interpreted as a coherent
    # FC identity: they are technical, contamination-related or isolated residence effects.
    sig["canonical_included"] = dataset == "zhuang"
    sig["interpretation"] = np.where(
        sig.canonical_included, "canonical matched DE",
        "audit-only: no coherent FC identity",
    )
    audit.append(sig[["dataset", "gene", "log2FC", "padj", "direction", "class",
                      "canonical_included", "interpretation"]])
    kept = sig[sig.canonical_included]
    canonical.append(kept[["dataset", "gene", "log2FC", "padj", "direction", "class"]])
    counts.append(dict(dataset=dataset, n_statistical_hits=len(sig), n_canonical_DE=len(kept)))

audit_df = pd.concat(audit, ignore_index=True)
canonical_df = pd.concat(canonical, ignore_index=True)
audit_df.to_csv(os.path.join(TAB, "DE_audit_candidates.csv"), index=False)
canonical_df.to_csv(os.path.join(TAB, "DE_canonical_summary.csv"), index=False)
pd.DataFrame(counts).to_csv(os.path.join(TAB, "DE_canonical_counts.csv"), index=False)

families_df = pd.DataFrame(
    [dict(family=name, n_genes=len(genes), genes=";".join(sorted(genes)),
          provenance=prov, defined_without_FC_vs_BG=indep)
     for name, genes, prov, indep in FAMILY_TABLE])
families_df.to_csv(os.path.join(TAB, "TableS_DE_audit_gene_sets.csv"), index=False)

print(pd.DataFrame(counts).to_string(index=False))
print(canonical_df.groupby(["dataset", "class"]).size())
print(families_df[["family", "n_genes"]].to_string(index=False))
