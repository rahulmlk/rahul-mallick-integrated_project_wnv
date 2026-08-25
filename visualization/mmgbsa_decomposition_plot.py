"""
mmgbsa_decomposition_plot.py

Standalone visualization utility: plots MM/GBSA binding free energy
decomposition for lead compounds, in the style of the figures generated
in my thesis (Figure 3.7/3.8, MM/GBSA binding free energy decomposition
for NS4B lead compounds).

Uses the thesis-reported total binding free energies as example data
(no confidential intermediate/raw data), decomposed into representative
energy terms for illustration.

Usage:
    python mmgbsa_decomposition_plot.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Reported total MM/GBSA binding free energies (kcal/mol) for the two
# thesis lead compounds against WNV NS4B.
COMPOUNDS = {
    "Compound 4366": -23.88,
    "Compound 4558": -19.34,
}

# Illustrative energy-term decomposition (van der Waals / electrostatic /
# polar solvation / nonpolar solvation) summing to each reported total.
# Proportions are representative, not the thesis's exact raw values.
DECOMPOSITION_FRACTIONS = {
    "van der Waals": 0.55,
    "Electrostatic": 0.30,
    "Polar solvation": -0.20,   # unfavorable term, opposes binding
    "Nonpolar solvation": 0.10,
}


def build_decomposition(total: float) -> dict[str, float]:
    return {term: round(total * frac, 2) for term, frac in DECOMPOSITION_FRACTIONS.items()}


def plot(out_path: Path) -> None:
    terms = list(DECOMPOSITION_FRACTIONS.keys())
    x = np.arange(len(terms))
    width = 0.35

    fig, ax = plt.subplots(figsize=(7, 5))
    for i, (compound, total) in enumerate(COMPOUNDS.items()):
        decomposition = build_decomposition(total)
        values = [decomposition[t] for t in terms]
        ax.bar(x + i * width, values, width, label=f"{compound} (ΔG = {total} kcal/mol)")

    ax.set_xticks(x + width / 2)
    ax.set_xticklabels(terms, rotation=20, ha="right")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Energy contribution (kcal/mol)")
    ax.set_title("MM/GBSA binding free energy decomposition: NS4B lead compounds")
    ax.legend()
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    plot(Path(__file__).resolve().parent.parent / "figures" / "mmgbsa_decomposition.png")
