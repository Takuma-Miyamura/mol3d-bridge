"""
Mol3D Bridge Python API – the same chemistry as the desktop app, without the window.

    from mol3d_bridge import Molecule3D

    m = Molecule3D.from_smiles("CC(=O)N[C@@H](Cc1ccccc1)C(=O)OC")
    m.minimize("MMFF94")                 # or "MMFF94s", "UFF", "GFN2-xTB", "GFN-FF"
    print(m.energy)                      # kcal/mol
    m.save("ac-phe-ome.xyz")             # .xyz .mol .sdf .pdb .inp (ORCA) .gjf/.com (Gaussian)
    print(m.orca_input("! r2SCAN-3c Opt Freq"))

    best = m.conformer_search(n=30)[0]   # list of Molecule3D, lowest energy first

Everything runs locally.  GFN2-xTB / GFN-FF need the free `xtb` program (on PATH, or pass
xtb_path=...).
"""
from __future__ import annotations

import os
from typing import List, Optional, Sequence

import numpy as np
from rdkit import Chem

from . import core

__all__ = ["Molecule3D", "METHODS"]

METHODS = core.METHODS


class Molecule3D:
    """A molecule with hydrogens and 3D coordinates (one conformer)."""

    def __init__(self, mol: Chem.Mol, source: Optional[Chem.Mol] = None,
                 energy: float = float("nan"), method: str = ""):
        if mol.GetNumConformers() == 0:
            raise ValueError("Molecule3D needs 3D coordinates – use Molecule3D.from_*().")
        if mol.GetNumConformers() > 1:
            mol = Chem.Mol(mol, confId=mol.GetConformers()[0].GetId())
        mol.GetConformer().SetId(0)
        self.mol = mol
        self.source = source if source is not None else Chem.Mol(mol)
        self.energy = energy
        self.method = method
        self.messages: List[str] = []

    # ------------------------------------------------------------------ construction
    @classmethod
    def _from_parsed(cls, m: Optional[Chem.Mol], what: str, seed: int = 0xF00D) -> "Molecule3D":
        if m is None or m.GetNumAtoms() == 0:
            raise ValueError(f"Could not read a structure from {what}.")
        m = Chem.Mol(m)
        Chem.SanitizeMol(m)
        core.normalize_scale(m)
        if core.is_3d(m):
            obj = cls(Chem.AddHs(m, addCoords=True), source=m)
        else:
            r = core.build_3d(m, None, 1, seed)
            obj = cls(r.mol, source=m)
            obj.messages = list(r.messages)
        return obj

    @classmethod
    def from_smiles(cls, smiles: str, seed: int = 0xF00D) -> "Molecule3D":
        """SMILES -> 3D model (ETKDG, not minimized)."""
        return cls._from_parsed(Chem.MolFromSmiles(smiles), f"SMILES {smiles!r}", seed)

    @classmethod
    def from_molblock(cls, molblock: str, seed: int = 0xF00D) -> "Molecule3D":
        """MDL molfile text (2D drawing or 3D) -> 3D model."""
        return cls._from_parsed(Chem.MolFromMolBlock(molblock, removeHs=False), "the molfile", seed)

    @classmethod
    def from_text(cls, text: str, seed: int = 0xF00D) -> "Molecule3D":
        """Any of: SMILES, molfile, CDXML, InChI, XYZ."""
        return cls._from_parsed(core.parse_text(text), "the text", seed)

    @classmethod
    def from_file(cls, path: str, seed: int = 0xF00D) -> "Molecule3D":
        """.cdx .cdxml .mol .sdf .xyz .pdb .smi"""
        return cls._from_parsed(core.load_file(path), os.path.basename(path), seed)

    @classmethod
    def from_rdkit(cls, mol: Chem.Mol, seed: int = 0xF00D) -> "Molecule3D":
        return cls._from_parsed(mol, "the RDKit molecule", seed)

    # ------------------------------------------------------------------ calculations
    def minimize(self, method: str = "MMFF94", fixed: Optional[Sequence[int]] = None,
                 xtb_path: str = "", solvent: str = "none", threads: int = 4,
                 charge: Optional[int] = None, uhf: Optional[int] = None) -> float:
        """Minimize from the current geometry. Returns the energy in kcal/mol."""
        if method not in METHODS:
            raise ValueError(f"method must be one of {METHODS}")
        fixed = list(fixed or [])
        if core.is_xtb(method):
            st = core.XtbSettings(path=core.find_xtb(xtb_path) or "", solvent=solvent,
                                  threads=threads, charge=charge, uhf=uhf)
            r = core.xtb_optimize(self.mol, 0, method, st, fixed)
        else:
            r = core.ff_trajectory(self.mol, 0, method, fixed)
        core.apply_positions(self.mol, 0, r.final)
        self.energy, self.method, self.messages = r.energy, r.method, list(r.messages)
        return self.energy

    def single_point(self, method: str = "MMFF94") -> float:
        """Force-field energy of the current geometry (kcal/mol), no optimisation."""
        return core.single_point(self.mol, 0, method)

    def trajectory(self, method: str = "MMFF94", fixed: Optional[Sequence[int]] = None,
                   max_step: float = 0.02) -> List[np.ndarray]:
        """Minimize and return every step's coordinates (N x 3 arrays) – e.g. for movies."""
        r = core.ff_trajectory(self.mol, 0, method, list(fixed or []), max_step=max_step)
        core.apply_positions(self.mol, 0, r.final)
        self.energy, self.method = r.energy, r.method
        return [f[0] for f in r.frames]

    def conformer_search(self, n: int = 30, method: str = "MMFF94",
                         rms_prune: float = 0.25, seed: int = 0xF00D) -> List["Molecule3D"]:
        """n starting geometries, each minimized; unique conformers sorted by energy."""
        ff = method if not core.is_xtb(method) else "MMFF94"
        r = core.build_3d(self.source, ff, n, seed, rms_prune=rms_prune, force_embed=True)
        out = []
        for c, e in zip(r.mol.GetConformers(), r.energies):
            m = Chem.Mol(r.mol, confId=c.GetId())
            out.append(Molecule3D(m, self.source, e, r.forcefield))
        return out

    def stereo_ok(self) -> bool:
        """True if the 3D model has the same stereochemistry as the input."""
        return core.stereo_changes(self.source, self.mol, self.mol.GetConformer().GetId()) == ""

    # ------------------------------------------------------------------ editing
    @property
    def coordinates(self) -> np.ndarray:
        return self.mol.GetConformer().GetPositions()

    @coordinates.setter
    def coordinates(self, xyz) -> None:
        core.apply_positions(self.mol, 0, np.asarray(xyz, float))
        self.energy, self.method = float("nan"), ""

    def set_internal(self, atoms: Sequence[int], value: float) -> None:
        """Distance (2 atoms, Å), angle (3, °) or dihedral (4, °); the last atom's side moves."""
        core.set_internal(self.mol, 0, list(atoms), value)
        self.energy, self.method = float("nan"), ""

    # ------------------------------------------------------------------ output
    @property
    def smiles(self) -> str:
        m = Chem.Mol(self.mol)
        Chem.AssignStereochemistryFrom3D(m)
        return core.stereo_smiles(m)

    @property
    def formula(self) -> str:
        return core.formula_and_mass(self.mol)[0]

    def charge_multiplicity(self):
        return core.default_charge_mult(self.mol)

    def to_xyz(self, comment: str = "") -> str:
        return core.xyz_block(self.mol, 0, comment or self._comment())

    def to_molblock(self) -> str:
        return core.molblock_3d(self.mol, 0)

    def to_pdb(self) -> str:
        return Chem.MolToPDBBlock(self.mol)

    def to_2d_molblock(self, rotation=None, clean: bool = False) -> str:
        """2D drawing (wedges from the 3D stereo): projection of `rotation` (3x3) or a clean layout."""
        R = np.eye(3) if rotation is None else np.asarray(rotation, float)
        return core.drawing_from_3d(self.mol, 0, R, "clean" if clean else "view").molblock

    def orca_input(self, keywords: str = "! r2SCAN-3c Opt Freq TightSCF", charge=None, mult=None,
                   nprocs: int = 4, maxcore: int = 2000, extra: str = "", title: str = "") -> str:
        ch, mu = self.charge_multiplicity()
        return core.orca_input(self.mol, 0, keywords, ch if charge is None else charge,
                               mu if mult is None else mult, nprocs, maxcore, extra, title)

    def gaussian_input(self, route: str = "#p opt freq b3lyp/6-31g(d) empiricaldispersion=gd3bj",
                       charge=None, mult=None, nprocs: int = 4, mem: str = "4GB",
                       chk: str = "", title: str = "", extra: str = "") -> str:
        ch, mu = self.charge_multiplicity()
        return core.gaussian_input(self.mol, 0, route, ch if charge is None else charge,
                                   mu if mult is None else mult, nprocs, mem, chk, title, extra)

    def save(self, path: str, **kw) -> str:
        """Write by extension: .xyz .mol .sdf .pdb .inp (ORCA) .gjf/.com (Gaussian)."""
        ext = os.path.splitext(path)[1].lower()
        if ext == ".xyz":
            text = self.to_xyz()
        elif ext == ".mol":
            text = self.to_molblock()
        elif ext == ".sdf":
            text = self.to_molblock() + "$$$$\n"
        elif ext == ".pdb":
            text = self.to_pdb()
        elif ext == ".inp":
            text = self.orca_input(**kw)
        elif ext in (".gjf", ".com"):
            kw.setdefault("chk", os.path.splitext(os.path.basename(path))[0] + ".chk")
            text = self.gaussian_input(**kw)
        else:
            raise ValueError(f"Unknown file type {ext!r}")
        with open(path, "w", newline="\n") as fh:
            fh.write(text)
        return path

    def _comment(self):
        e = "" if np.isnan(self.energy) else f"  {self.method} E={core.format_energy(self.energy, self.method)}"
        return self.formula + e

    def __repr__(self):
        e = "not minimized" if np.isnan(self.energy) else f"{self.method} {self.energy:.3f} kcal/mol"
        return f"<Molecule3D {self.formula}, {self.mol.GetNumAtoms()} atoms, {e}>"
