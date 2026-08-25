"""
trajectory_analysis.py

Representative MD trajectory analysis reflecting the GROMACS + MDAnalysis
workflow used in my thesis for membrane-embedded NS4B simulations and
glycosylated/non-glycosylated NS1 comparisons: backbone RMSD, per-residue
RMSF, radius of gyration over the trajectory, and ligand hydrogen-bond
occupancy.

Runs against data/toy_trajectory.pdb, a synthetic 30-frame multi-model PDB
standing in for a GROMACS .xtc/.tpr pair (not a real simulation) so the
script is self-contained and dependency-light. Point it at a real
topology + trajectory to use it for actual MD output.

Usage:
    python trajectory_analysis.py [path/to/trajectory.pdb]
"""
from __future__ import annotations

import sys
from pathlib import Path

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import align, rms

DEFAULT_TRAJECTORY = Path(__file__).resolve().parent / "data" / "toy_trajectory.pdb"
HBOND_DISTANCE_CUTOFF_A = 3.5  # donor/acceptor heavy-atom distance, GROMACS-style cutoff
EQUILIBRATION_FRAMES = 10      # frames treated as pre-equilibration, excluded from RMSF


def load_universe(path: Path) -> mda.Universe:
    return mda.Universe(str(path))


def backbone_rmsd(u: mda.Universe, selection: str = "name CA") -> np.ndarray:
    ref = u.copy()
    ref.trajectory[0]
    r = rms.RMSD(u, ref, select=selection)
    r.run()
    return r.results.rmsd[:, 2]  # time, frame, RMSD columns -> take RMSD


def per_residue_rmsf(u: mda.Universe, selection: str = "name CA") -> np.ndarray:
    ca = u.select_atoms(selection)
    aligner = align.AlignTraj(u, u, select=selection, in_memory=True)
    aligner.run()
    from MDAnalysis.analysis.rms import RMSF

    # skip early "equilibration" frames, matching thesis convention of
    # computing RMSF over the equilibrated portion of the trajectory
    rmsf_calc = RMSF(ca).run(start=EQUILIBRATION_FRAMES)
    return rmsf_calc.results.rmsf


def radius_of_gyration_over_time(u: mda.Universe) -> np.ndarray:
    protein = u.select_atoms("name CA")
    rg_values = []
    for _ in u.trajectory:
        rg_values.append(protein.radius_of_gyration())
    return np.array(rg_values)


def ligand_hbond_occupancy(u: mda.Universe) -> float:
    """Fraction of frames where the ligand pseudo-atom is within
    HBOND_DISTANCE_CUTOFF_A of any CA atom -- a simplified stand-in for
    the donor/acceptor hydrogen-bond occupancy analysis used in the
    thesis (there, real polar atoms and GROMACS-style geometric criteria
    were used; here a single heavy-atom distance is used since the toy
    trajectory only has CA-level detail)."""
    ligand = u.select_atoms("resname LIG")
    protein = u.select_atoms("name CA")
    if len(ligand) == 0:
        return 0.0

    contacts = 0
    n_frames = len(u.trajectory)
    for _ in u.trajectory:
        dists = np.linalg.norm(protein.positions - ligand.positions[0], axis=1)
        if dists.min() <= HBOND_DISTANCE_CUTOFF_A:
            contacts += 1
    return round(100 * contacts / n_frames, 1)


def main(argv: list[str]) -> int:
    traj_path = Path(argv[1]) if len(argv) > 1 else DEFAULT_TRAJECTORY
    u = load_universe(traj_path)

    print(f"Trajectory: {traj_path.name}")
    print(f"  Atoms: {len(u.atoms)}, Frames: {len(u.trajectory)}\n")

    rmsd = backbone_rmsd(u)
    print(f"Backbone RMSD: mean={rmsd.mean():.2f} A, "
          f"final={rmsd[-1]:.2f} A, "
          f"first-{EQUILIBRATION_FRAMES}-frame mean={rmsd[:EQUILIBRATION_FRAMES].mean():.2f} A "
          f"vs. post-equilibration mean={rmsd[EQUILIBRATION_FRAMES:].mean():.2f} A")

    rmsf = per_residue_rmsf(u)
    print(f"\nPer-residue RMSF (equilibrated portion, residues 1-{len(rmsf)}):")
    print(f"  mean={rmsf.mean():.2f} A, max={rmsf.max():.2f} A "
          f"at residue {int(np.argmax(rmsf)) + 1}")

    rg = radius_of_gyration_over_time(u)
    print(f"\nRadius of gyration: mean={rg.mean():.2f} A, "
          f"std={rg.std():.2f} A (stable if std is small relative to mean)")

    occupancy = ligand_hbond_occupancy(u)
    print(f"\nLigand-protein contact occupancy (proxy for H-bond occupancy): "
          f"{occupancy}% of frames")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
