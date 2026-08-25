"""
protein_structure_analysis.py

Representative structural QC analysis reflecting the structural
bioinformatics checks used in my thesis when comparing glycosylated vs.
non-glycosylated WNV NS1 homodimer models (radius of gyration, solvent
accessibility, disulfide connectivity) before committing structures to
molecular dynamics.

Runs against data/sample_structure.pdb, a small synthetic toy helix (not
a real viral protein) included so the script is fully self-contained.
Point PARSER at any real PDB/mmCIF file to use it for real work.

Usage:
    python protein_structure_analysis.py [path/to/structure.pdb]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

DEFAULT_STRUCTURE = Path(__file__).resolve().parent / "data" / "sample_structure.pdb"
DISULFIDE_SG_SG_CUTOFF_A = 2.5  # typical S-S bond length ~2.05 A; allow some slack


def load_structure(pdb_path: Path):
    parser = PDBParser(QUIET=True)
    return parser.get_structure("model", str(pdb_path))


def radius_of_gyration(structure) -> float:
    coords = np.array([atom.coord for atom in structure.get_atoms()])
    center = coords.mean(axis=0)
    diffs = coords - center
    rg = np.sqrt((diffs ** 2).sum(axis=1).mean())
    return round(float(rg), 2)


def bfactor_stats(structure) -> dict:
    bfactors = np.array([atom.get_bfactor() for atom in structure.get_atoms()])
    return {
        "mean": round(float(bfactors.mean()), 2),
        "min": round(float(bfactors.min()), 2),
        "max": round(float(bfactors.max()), 2),
    }


def total_sasa(structure) -> float:
    sr = ShrakeRupley()
    sr.compute(structure, level="S")
    return round(float(structure.sasa), 2)


def find_disulfides(structure) -> list[tuple[int, int, float]]:
    sg_atoms = [
        (res.id[1], atom)
        for res in structure.get_residues()
        if res.get_resname() == "CYS"
        for atom in res
        if atom.get_name() == "SG"
    ]
    pairs = []
    for i in range(len(sg_atoms)):
        for j in range(i + 1, len(sg_atoms)):
            resi_i, atom_i = sg_atoms[i]
            resi_j, atom_j = sg_atoms[j]
            dist = atom_i - atom_j
            if dist <= DISULFIDE_SG_SG_CUTOFF_A:
                pairs.append((resi_i, resi_j, round(float(dist), 2)))
    return pairs


def main(argv: list[str]) -> int:
    pdb_path = Path(argv[1]) if len(argv) > 1 else DEFAULT_STRUCTURE
    structure = load_structure(pdb_path)

    n_residues = sum(1 for _ in structure.get_residues())
    n_atoms = sum(1 for _ in structure.get_atoms())

    print(f"Structure: {pdb_path.name}")
    print(f"  Residues: {n_residues}")
    print(f"  Atoms:    {n_atoms}")
    print(f"  Radius of gyration: {radius_of_gyration(structure)} A")

    bf = bfactor_stats(structure)
    print(f"  B-factor (mean/min/max): {bf['mean']} / {bf['min']} / {bf['max']}")

    print(f"  Total SASA: {total_sasa(structure)} A^2")

    disulfides = find_disulfides(structure)
    if disulfides:
        print(f"  Candidate disulfide bonds (SG-SG <= {DISULFIDE_SG_SG_CUTOFF_A} A):")
        for resi_i, resi_j, dist in disulfides:
            print(f"    CYS{resi_i} - CYS{resi_j}: {dist} A")
    else:
        print("  No candidate disulfide bonds found within cutoff.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
