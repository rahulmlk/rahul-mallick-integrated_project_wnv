# Molecular Dynamics: Trajectory Analysis

Representative MD trajectory analysis reflecting the GROMACS + MDAnalysis
workflow used in my thesis for membrane-embedded WNV NS4B simulations
(lead compound binding, POPC bilayer, CHARMM36m) and glycosylated vs.
non-glycosylated NS1 homodimer comparisons: backbone RMSD, per-residue
RMSF (over the equilibrated portion of the trajectory), radius of
gyration, and ligand hydrogen-bond occupancy.

## Real pipeline (thesis)

- Simulations run in **GROMACS** with **CHARMM36m** (protein/membrane)
  and **CGenFF** (ligand) force fields, **OPLS4**/**OPLS3** used
  elsewhere in the docking/refinement stage
- Membrane systems built with **CHARMM-GUI Membrane Builder**
- Trajectory analysis with **GROMACS analysis tools** and **MDAnalysis**
- Metrics: backbone RMSD, per-residue RMSF (post-equilibration),
  radius of gyration, intermolecular hydrogen bond occupancy, PCA of
  C-alpha displacement
- Binding free energy: **MM/GBSA**

## Files

- `trajectory_analysis.py` - RMSD, RMSF, Rg, and a simplified ligand
  contact-occupancy proxy for hydrogen bonding
- `data/toy_trajectory.pdb` - a synthetic 30-frame multi-MODEL PDB
  standing in for a GROMACS `.xtc`/`.tpr` pair, with a small drift-then-
  equilibrate pattern and a persistent nearby "ligand" pseudo-atom

This uses a synthetic trajectory (not real GROMACS output) so the
analysis code is fully runnable without a multi-gigabyte trajectory file.
The same MDAnalysis calls work directly against real GROMACS output —
just point `Universe()` at your `.tpr`/`.xtc` files instead.

## Running

```bash
pip install MDAnalysis numpy
python trajectory_analysis.py                          # uses the toy trajectory
python trajectory_analysis.py /path/to/real_traj.pdb    # or real MD output
```

## Note on scope

MM/GBSA binding free energy calculation and PCA of collective motions
(both used in the thesis) aren't included here to keep this example
focused and dependency-light; `trajectory_analysis.py` covers the
RMSD/RMSF/Rg/H-bond metrics that make up the bulk of day-to-day
trajectory QC.
