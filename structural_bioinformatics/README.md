# Structural Bioinformatics: Protein Structure QC

Representative structural QC checks reflecting the kind of validation used
before committing a modeled/docked protein structure to molecular dynamics
in my thesis (e.g. comparing glycosylated vs. non-glycosylated WNV NS1
homodimer models: radius of gyration, solvent accessible surface area,
disulfide connectivity).

## Files

- `protein_structure_analysis.py` - loads a PDB structure and reports
  residue/atom counts, radius of gyration, B-factor distribution, total
  SASA (Shrake-Rupley algorithm), and candidate disulfide bonds
- `data/sample_structure.pdb` - a small synthetic 20-residue toy helix
  (not a real viral protein) generated purely to exercise the code

## Real pipeline (thesis)

- 3D structure modeling: **trRosetta** + **GalaxyRefine2**
- Structure validation: **ProSA-web** (Z-score)
- Glycosylation: **CHARMM-GUI Glycan Reader**
- Force field parameterization: **OPLS4** / **CHARMM36m**

Those steps depend on external web servers/GPU resources and aren't
reproducible offline, so this script focuses on the structure-level QC
metrics (Rg, SASA, B-factor, disulfides) that are computable directly from
any PDB file with Biopython.

## Running

```bash
pip install biopython numpy
python protein_structure_analysis.py                       # uses the toy structure
python protein_structure_analysis.py /path/to/real.pdb      # or any real PDB
```

Note: the toy helix's two cysteines are positioned ~19 A apart along the
helix axis, so the disulfide check correctly reports none found — the
script isn't tuned to force a hit, it reports what the geometry says.
