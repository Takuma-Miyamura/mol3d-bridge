"""
m3b – Mol3D Bridge from the command line.

  m3b "CC(=O)N[C@@H](Cc1ccccc1)C(=O)OC" -o phe.xyz
  m3b drawing.cdx --method GFN2-xTB --orca job.inp --keywords "! r2SCAN-3c Opt Freq"
  m3b input.mol --conformers 30 -o confs.sdf
  m3b "c1ccccc1O" --gaussian phenol.gjf --route "#p opt freq m062x/6-311+g(d,p)"
"""
from __future__ import annotations

import argparse
import os
import sys

from . import __version__


def main(argv=None) -> int:
    from rdkit import RDLogger

    from . import core
    from .api import METHODS, Molecule3D

    RDLogger.DisableLog("rdApp.*")      # clean terminal output; errors are reported below

    ap = argparse.ArgumentParser(prog="m3b", description="Build molecules in 3D and minimize them "
                                 "(MMFF94, MMFF94s, UFF, GFN2-xTB, GFN-FF).",
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("structure", help="SMILES, or a file: .cdx .cdxml .mol .sdf .xyz .pdb .smi")
    ap.add_argument("--method", "-m", default="MMFF94", choices=METHODS)
    ap.add_argument("--no-minimize", action="store_true", help="only build the 3D model")
    ap.add_argument("--conformers", "-c", type=int, default=0, metavar="N",
                    help="conformer search with N starts (MMFF); with -o .sdf all are written")
    ap.add_argument("--fixed", default="", help="atoms to hold fixed, 0-based, e.g. 0,1,5")
    ap.add_argument("--xtb", default="", help="path to the xtb program")
    ap.add_argument("--solvent", default="none", help="ALPB solvent for xtb (e.g. water)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("-o", "--output", help="write .xyz .mol .sdf .pdb")
    ap.add_argument("--orca", metavar="FILE.inp", help="write an ORCA input")
    ap.add_argument("--keywords", default="! r2SCAN-3c Opt Freq TightSCF")
    ap.add_argument("--gaussian", metavar="FILE.gjf", help="write a Gaussian input")
    ap.add_argument("--route", default="#p opt freq b3lyp/6-31g(d) empiricaldispersion=gd3bj")
    ap.add_argument("--charge", type=int)
    ap.add_argument("--mult", type=int)
    ap.add_argument("--nprocs", type=int, default=4)
    ap.add_argument("--version", action="version", version=f"Mol3D Bridge {__version__}")
    a = ap.parse_args(argv)

    try:
        if os.path.exists(a.structure):
            m = Molecule3D.from_file(a.structure)
        else:
            m = Molecule3D.from_text(a.structure)
        confs = None
        if a.conformers > 0:
            confs = m.conformer_search(a.conformers, a.method)
            m = confs[0]
            print(f"{len(confs)} conformers; lowest {m.energy:.3f} kcal/mol ({m.method})", file=sys.stderr)
        if not a.no_minimize and (confs is None or core.is_xtb(a.method)):
            fixed = [int(x) for x in a.fixed.split(",") if x.strip()]
            m.minimize(a.method, fixed, xtb_path=a.xtb, solvent=a.solvent, threads=a.threads)
        for msg in m.messages:
            print("note:", msg, file=sys.stderr)
        if not m.stereo_ok():
            print("warning: stereochemistry differs from the input", file=sys.stderr)
        print(repr(m), file=sys.stderr)
        if a.output:
            if a.output.lower().endswith(".sdf") and confs:
                from rdkit import Chem
                w = Chem.SDWriter(a.output)
                for c in confs:
                    c.mol.SetProp("Energy_kcal_mol", f"{c.energy:.4f}")
                    w.write(c.mol)
                w.close()
            else:
                m.save(a.output)
        qm = dict(charge=a.charge, mult=a.mult, nprocs=a.nprocs)
        if a.orca:
            m.save(a.orca, keywords=a.keywords, **qm)
        if a.gaussian:
            m.save(a.gaussian, route=a.route, **qm)
        if not (a.output or a.orca or a.gaussian):
            sys.stdout.write(m.to_xyz())
        return 0
    except Exception as ex:  # noqa: BLE001
        print(f"m3b: error: {ex}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
