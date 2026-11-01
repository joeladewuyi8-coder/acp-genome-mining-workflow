"""
rank_peptides.py

Screen and rank candidate peptides using predicted properties.

Input: a CSV file with one row per peptide and these columns:
    peptide_id, net_charge, hydrophobic_ratio, boman_index,
    gravy, instability_index, half_life_h, anticp_score

Steps:
    1. Screen: keep peptides with net charge > 0, instability index < 40
       and Boman index < 1 kcal/mol.
    2. Min-max scale each property to 0-1 across the peptides that passed.
       Boman and instability indices are inverted and half-life is
       log-transformed, so that higher is always better.
    3. Composite score = 0.7 * membrane-activity score
                       + 0.3 * stability score
       Membrane activity = mean of scaled net charge, hydrophobic ratio,
       GRAVY and inverted Boman index.
       Stability = mean of inverted instability index and log half-life.
    4. Report the Spearman correlation between the composite score and
       the AntiCP score (the AntiCP score is NOT part of the composite).

Usage:
    python scripts/rank_peptides.py input.csv output.csv

All values are computational predictions, not experimental results.
"""

import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REQUIRED = [
    "peptide_id", "net_charge", "hydrophobic_ratio", "boman_index",
    "gravy", "instability_index", "half_life_h", "anticp_score",
]
W_MEMBRANE, W_STABILITY = 0.7, 0.3


def min_max(series):
    """Scale to 0-1. If all values are equal, return 0.5 for every row."""
    span = series.max() - series.min()
    if span == 0:
        return pd.Series(0.5, index=series.index)
    return (series - series.min()) / span


def screen(df):
    keep = (
        (df["net_charge"] > 0)
        & (df["instability_index"] < 40)
        & (df["boman_index"] < 1)
    )
    return df[keep].copy()


def score(df):
    membrane = (
        min_max(df["net_charge"])
        + min_max(df["hydrophobic_ratio"])
        + min_max(df["gravy"])
        + (1 - min_max(df["boman_index"]))
    ) / 4
    stability = (
        (1 - min_max(df["instability_index"]))
        + min_max(np.log(df["half_life_h"]))
    ) / 2
    df["composite_score"] = W_MEMBRANE * membrane + W_STABILITY * stability
    return df.sort_values("composite_score", ascending=False)


def main(in_path, out_path):
    df = pd.read_csv(in_path)
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        sys.exit(f"Missing columns: {', '.join(missing)}")

    passed = screen(df)
    print(f"{len(df)} peptides read, {len(passed)} passed the screen.")
    if len(passed) < 3:
        sys.exit("Too few peptides passed the screen to rank.")

    ranked = score(passed)
    ranked["rank"] = range(1, len(ranked) + 1)
    ranked.round(2).to_csv(out_path, index=False)

    rho, p = spearmanr(ranked["composite_score"], ranked["anticp_score"])
    print(f"Spearman rho (composite vs AntiCP) = {rho:.2f}, p = {p:.2f}")
    print("Note: tied AntiCP scores weaken this correlation.")
    print(f"Ranked table written to {out_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: python scripts/rank_peptides.py input.csv output.csv")
    main(sys.argv[1], sys.argv[2])
