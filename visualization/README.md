# Visualization: MM/GBSA Decomposition Plot

Standalone plotting utility reflecting the MM/GBSA binding free energy
decomposition figures from my thesis (NS4B lead compound binding).

## Files

- `mmgbsa_decomposition_plot.py` - generates a grouped bar chart of
  energy-term contributions for the two thesis lead compounds, using the
  reported total binding free energies (ΔG = −23.88 and −19.34 kcal/mol)
  decomposed into representative van der Waals / electrostatic / polar
  solvation / nonpolar solvation terms

The decomposition proportions are illustrative (for plotting purposes
only) — the exact per-term raw values from the thesis are not included
here; only the reported totals are real.

## Running

```bash
pip install matplotlib numpy
python mmgbsa_decomposition_plot.py
```

Writes `../figures/mmgbsa_decomposition.png`.
