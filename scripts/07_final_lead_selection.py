"""
Stage 7 — Final lead selection: combine QSAR score, docking affinity, and
confirmed zinc coordination into one defensible ranking.

This is the direct answer to "does docking score alone matter?" — it
doesn't. A compound with the single best raw docking affinity of the
whole set (CHEMBL323371, -11.58 kcal/mol) has only borderline zinc
coordination and drops out of the top 5 once that's weighed in, while
several compounds with more modest affinities but confirmed catalytic
engagement rank higher.

final_score = 0.4 * normalized QSAR composite score
            + 0.4 * normalized docking affinity
            + 0.2 * confirmed zinc-coordination bonus
"""

import os

import pandas as pd

from config import PROJECT_DIR


def main():
    composite_df = pd.read_csv(os.path.join(PROJECT_DIR, "MMP9_composite_ranking_full.csv"))
    docking_df = pd.read_csv(os.path.join(PROJECT_DIR, "MMP9_docking_results.csv"))
    zinc_df = pd.read_csv(os.path.join(PROJECT_DIR, "MMP9_zinc_interaction_check.csv"))

    merged = composite_df.merge(docking_df, on="compound_id").merge(zinc_df, on="compound_id")

    merged["norm_composite"] = (
        (merged["composite_score"] - merged["composite_score"].min())
        / (merged["composite_score"].max() - merged["composite_score"].min())
    )
    merged["norm_affinity"] = (
        (-merged["best_affinity_kcal_mol"] - (-merged["best_affinity_kcal_mol"]).min())
        / ((-merged["best_affinity_kcal_mol"]).max() - (-merged["best_affinity_kcal_mol"]).min())
    )
    merged["zn_bonus"] = merged["zinc_coordinating"].astype(float)

    merged["final_score"] = (
        0.4 * merged["norm_composite"]
        + 0.4 * merged["norm_affinity"]
        + 0.2 * merged["zn_bonus"]
    )

    final = merged.sort_values("final_score", ascending=False).reset_index(drop=True)
    final.insert(0, "final_rank", range(1, len(final) + 1))

    out_path = os.path.join(PROJECT_DIR, "MMP9_FINAL_lead_selection.csv")
    final.to_csv(out_path, index=False)
    print(f"Saved: {out_path}\n")

    print(final[[
        "final_rank", "compound_id", "composite_rank", "predicted_pIC50",
        "best_affinity_kcal_mol", "zinc_coordinating", "dist_to_catalytic_zn", "final_score"
    ]].head(10).to_string(index=False))

    return final


if __name__ == "__main__":
    main()
