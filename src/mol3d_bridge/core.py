"""
Mol3D Bridge - core chemistry (no GUI code here, so it can be tested/scripted).

Responsibilities
  * parse whatever ChemDraw / ReShiki put on the clipboard (CDX, CDXML, MOL, SMILES, ...)
  * build 3D coordinates (ETKDG) and pre-optimize them with MMFF94 / MMFF94s / UFF
  * project the current 3D view back to a 2D drawing with wedge/hash bonds
  * export coordinates (ORCA input, XYZ, MOL/SDF, PDB)
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdDepictor, rdMolDescriptors

RDLogger.DisableLog("rdApp.warning")

# Writing binary CDX needs Chem.MolToCDXMLBlock (newer RDKit, Python 3.10+).  Older RDKit
# (e.g. 2025.09 on Python 3.9) still reads ChemDraw data; copies to ChemDraw then use MDL molfiles.
HAS_CDX = bool(getattr(Chem, "HasChemDrawCDXSupport", lambda: False)()) and hasattr(Chem, "MolToCDXMLBlock")
HAS_CDXML_READ = hasattr(Chem, "MolsFromCDXML")


# --------------------------------------------------------------------------------------
# MDLCT (the "MDL Clipboard Text" format ChemDraw/ISIS use on Windows):
# every molfile line is prefixed with one length byte, no line terminators.
# --------------------------------------------------------------------------------------
def mdlct_encode(molblock: str) -> bytes:
    out = bytearray()
    for line in molblock.replace("\r\n", "\n").split("\n"):
        raw = line.encode("ascii", "replace")[:255]
        out.append(len(raw))
        out += raw
    return bytes(out)


def molfile_text_from_bytes(data: bytes) -> str | None:
    """Accept plain molfile text or MDLCT length-prefixed lines."""
    try:
        text = data.decode("utf-8")
        if "M  END" in text and "\n" in text:
            return text.rstrip("\x00")
    except UnicodeDecodeError:
        pass
    lines, i = [], 0
    data = data.rstrip(b"\x00") + b""  # tolerate trailing NUL padding
    while i < len(data):
        n = data[i]
        line = data[i + 1 : i + 1 + n]
        if len(line) < n:
            break
        lines.append(line.decode("latin-1"))
        i += 1 + n
    text = "\n".join(lines)
    if "M  END" not in text:
        # trailing NUL stripping may have eaten zero-length lines at the end
        return None
    return text + "\n"


# --------------------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------------------
def _combine(mols) -> Chem.Mol | None:
    mols = [m for m in mols if m is not None and m.GetNumAtoms() > 0]
    if not mols:
        return None
    out = mols[0]
    for m in mols[1:]:
        out = Chem.CombineMols(out, m)
    out = Chem.Mol(out)
    try:
        Chem.SanitizeMol(out)
    except Exception:
        pass
    return out


def parse_cdx(data: bytes | str) -> Chem.Mol | None:
    if not HAS_CDXML_READ:
        raise RuntimeError("This RDKit build cannot read ChemDraw data (update RDKit).")
    # 3D positions are dropped on purpose: RDKit reads them with the wrong order/scale/handedness.
    # The drawing (2D positions + wedges) is exact; the 3D model is rebuilt from it.
    return _combine(Chem.MolsFromCDXML(strip_cdx_3d(data), sanitize=True, removeHs=False))


def parse_text(text: str) -> Chem.Mol | None:
    t = text.strip().strip("\x00")
    if not t:
        return None
    if "M  END" in t:
        return Chem.MolFromMolBlock(text, removeHs=False)
    if t.lstrip().startswith("<?xml") or "<CDXML" in t[:400]:
        return parse_cdx(t)
    if t.startswith("InChI="):
        return Chem.MolFromInchi(t.splitlines()[0])
    first = t.splitlines()[0].strip()
    if first.isdigit() and len(t.splitlines()) >= int(first) + 2:
        return mol_from_xyz(t)
    # SMILES (first token of first line)
    return Chem.MolFromSmiles(first.split()[0])


def mol_from_xyz(text: str, charge: int = 0) -> Chem.Mol | None:
    from rdkit.Chem import rdDetermineBonds

    m = Chem.MolFromXYZBlock(text)
    if m is None:
        return None
    try:
        rdDetermineBonds.DetermineBonds(m, charge=charge)
    except Exception:
        rdDetermineBonds.DetermineConnectivity(m)
    return m


def parse_clipboard_item(kind: str, data: bytes) -> Chem.Mol | None:
    if kind == "cdx":
        return parse_cdx(data)
    if kind == "cdxml":
        return parse_cdx(data.rstrip(b"\x00").decode("utf-8", "replace"))
    if kind == "molfile":
        txt = molfile_text_from_bytes(data)
        return Chem.MolFromMolBlock(txt, removeHs=False) if txt else None
    if kind in ("smiles", "text"):
        return parse_text(data.rstrip(b"\x00").decode("utf-8", "replace"))
    return None


def load_file(path: str) -> Chem.Mol | None:
    p = path.lower()
    if p.endswith((".cdx", ".cdxml")):
        with open(path, "rb") as fh:
            return parse_cdx(fh.read())
    if p.endswith((".mol", ".mdl")):
        return Chem.MolFromMolFile(path, removeHs=False)
    if p.endswith(".sdf"):
        sup = Chem.SDMolSupplier(path, removeHs=False)
        mols = [m for m in sup if m is not None]
        if not mols:
            return None
        base = Chem.Mol(mols[0])
        # several records of the same molecule -> keep as conformers
        same = all(m.GetNumAtoms() == base.GetNumAtoms() for m in mols)
        if same and len(mols) > 1:
            base.RemoveAllConformers()
            for m in mols:
                base.AddConformer(Chem.Conformer(m.GetConformer()), assignId=True)
        return base
    if p.endswith(".xyz"):
        with open(path, encoding="utf-8", errors="replace") as fh:
            return mol_from_xyz(fh.read())
    if p.endswith(".pdb"):
        return Chem.MolFromPDBFile(path, removeHs=False)
    with open(path, encoding="utf-8", errors="replace") as fh:
        return parse_text(fh.read())


def is_3d(mol: Chem.Mol) -> bool:
    """
    True only for real 3D coordinates.  A drawing that merely carries some depth values
    (perspective / projection hints, as ReShiki or ChemDraw may write) is treated as 2D:
    in real 3D all bond lengths are consistent, in a drawing with fake depth they are not.
    """
    if mol.GetNumConformers() == 0 or mol.GetNumAtoms() < 2:
        return False
    pos = mol.GetConformer().GetPositions()
    if np.ptp(pos[:, 2]) <= 0.05:
        return False
    if mol.GetNumBonds() == 0:
        return True
    bl = np.array([np.linalg.norm(pos[b.GetBeginAtomIdx()] - pos[b.GetEndAtomIdx()])
                   for b in mol.GetBonds()])
    med = float(np.median(bl))
    if med <= 1e-9:
        return False
    if np.ptp(pos[:, 2]) < 0.25 * med:
        return False
    zc = pos[:, 2] - np.median(pos[:, 2])
    if np.mean(np.abs(zc) > 0.05 * med) < 0.3:   # only a few atoms carry depth: a drawing
        return False
    # heavy-heavy bonds of a real 3D model vary by ~±25 % at most
    return float(np.percentile(np.abs(bl / med - 1), 90)) < 0.35


# --------------------------------------------------------------------------------------
# 3D building and force fields
# --------------------------------------------------------------------------------------
FORCE_FIELDS = ("MMFF94", "MMFF94s", "UFF")


@dataclass
class BuildResult:
    mol: Chem.Mol
    energies: list = field(default_factory=list)  # per conformer, kcal/mol
    forcefield: str = ""
    messages: list = field(default_factory=list)


def stereo_smiles(mol: Chem.Mol) -> str:
    m = Chem.Mol(mol)
    try:
        m = Chem.RemoveHs(m)
    except Exception:
        pass
    return Chem.MolToSmiles(m)


def _ff_available(mol, ff):
    if ff.startswith("MMFF"):
        return AllChem.MMFFHasAllMoleculeParams(mol)
    return AllChem.UFFHasAllMoleculeParams(mol)


def _make_ff(mol, ff, conf_id):
    if ff.startswith("MMFF"):
        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=ff)
        return AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
    return AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)


def choose_ff(mol, requested, messages):
    if _ff_available(mol, requested):
        return requested
    for alt in FORCE_FIELDS:
        if alt != requested and _ff_available(mol, alt):
            messages.append(f"{requested} has no parameters for this molecule; used {alt}.")
            return alt
    messages.append("No force-field parameters available; geometry is the raw embedding.")
    return ""


def minimize(mol, ff, conf_id=-1, max_iters=2000, fixed=None):
    """Minimize one conformer in place (optionally with frozen atoms). Returns (energy, converged)."""
    field_ = _make_ff(mol, ff, conf_id)
    if field_ is None:
        return float("nan"), False
    for i in fixed or ():
        field_.AddFixedPoint(int(i))
    field_.Initialize()
    not_conv = 1
    for _ in range(10):  # restart the BFGS a few times like Chem3D's iterative MM2
        not_conv = field_.Minimize(maxIts=max_iters)
        if not_conv == 0:
            break
    return field_.CalcEnergy(), not_conv == 0


def build_3d(
    mol_in: Chem.Mol,
    ff: str = "MMFF94",
    n_confs: int = 1,
    seed: int = 0xF00D,
    max_iters: int = 2000,
    rms_prune: float = 0.25,
    force_embed: bool = False,
) -> BuildResult:
    """Drawing -> 3D: add H, ETKDGv3 embedding, force-field pre-optimisation."""
    msgs = []
    ref = stereo_smiles(mol_in)
    has3d = is_3d(mol_in)
    mol = Chem.AddHs(mol_in, addCoords=has3d)

    if has3d and n_confs <= 1 and not force_embed:
        conf_ids = [mol.GetConformer().GetId()]
    else:
        conf_ids = embed_cascade(mol, n_confs, seed, rms_prune, msgs)
        if not conf_ids:
            # last resort: the drawing itself, lifted into 3D (always works)
            flat = from_drawing(mol_in, seed)
            mol = flat
            conf_ids = [0]
            msgs.append("3D embedding failed; started from the drawing instead (minimize it).")

    used = choose_ff(mol, ff, msgs) if ff else ""
    energies = []
    for cid in conf_ids:
        if used:
            e, conv = minimize(mol, used, cid, max_iters)
            if not conv:
                msgs.append(f"Conformer {cid}: not fully converged.")
        else:
            e = float("nan")
        energies.append(e)

    # sort conformers by energy, renumber 0..n-1
    order = sorted(range(len(conf_ids)), key=lambda i: (math.isnan(energies[i]), energies[i]))
    confs = [Chem.Conformer(mol.GetConformer(conf_ids[i])) for i in order]
    energies = [energies[i] for i in order]
    mol.RemoveAllConformers()
    for i, c in enumerate(confs):
        c.SetId(i)
        c.Set3D(True)
        mol.AddConformer(c, assignId=False)

    _separate_fragments(mol)
    if not has3d and mol_in.GetNumConformers():
        for c in mol.GetConformers():
            align_to_drawing(mol, c.GetId(), mol_in)

    # stereo sanity check (drawing vs. 3D)
    chk = Chem.Mol(mol, confId=0)
    Chem.AssignStereochemistryFrom3D(chk, confId=0)
    if ref and stereo_smiles(chk) != ref and "@" in ref + stereo_smiles(chk):
        msgs.append("Warning: stereochemistry of the 3D model differs from the drawing.")
    return BuildResult(mol, energies, used, msgs)


# --------------------------------------------------------------------------------------
# robust embedding and "as drawn" 3D
# --------------------------------------------------------------------------------------
def embed_cascade(mol, n_confs=1, seed=0xF00D, rms_prune=0.25, msgs=None, cancel=None):
    """
    ETKDG first, then progressively more tolerant distance geometry.  Each step has a time
    limit, so strained or unusual structures (cycloparaphenylenes, cages, metal complexes,
    pre-cleaned drawings with odd stereo flags) never hang the program.
    """
    msgs = [] if msgs is None else msgs
    n = max(1, n_confs)

    def etkdg(random_coords):
        p = AllChem.ETKDGv3()
        p.useSmallRingTorsions = True
        p.useMacrocycleTorsions = True
        p.useRandomCoords = random_coords
        p.timeout = (3 if random_coords else 4) if n == 1 else max(5, min(60, n))
        return p

    def plain(enforce):
        p = AllChem.EmbedParameters()
        p.useExpTorsionAnglePrefs = False
        p.useBasicKnowledge = False
        p.useRandomCoords = True
        p.ignoreSmoothingFailures = True
        p.enforceChirality = enforce
        p.timeout = 10
        return p

    attempts = [("ETKDG", etkdg(False)), ("ETKDG (random start)", etkdg(True)),
                ("plain distance geometry", plain(True)),
                ("distance geometry without stereo enforcement", plain(False))]
    for k, (name, p) in enumerate(attempts):
        if cancel is not None and cancel.is_set():
            break
        p.randomSeed = seed
        p.pruneRmsThresh = rms_prune if n > 1 else -1.0
        p.numThreads = 0
        try:
            ids = [i for i in AllChem.EmbedMultipleConfs(mol, numConfs=n, params=p) if i >= 0]
        except Exception:  # noqa: BLE001
            ids = []
        if ids:
            if k:
                msgs.append(f"Embedded with {name} (strained or unusual structure).")
            return ids
    return []


def align_to_drawing(mol3d: Chem.Mol, conf_id: int, mol2d: Chem.Mol):
    """Rotate a 3D model (proper rotation only) so that it looks like the 2D drawing."""
    n = mol2d.GetNumAtoms()
    if mol2d.GetNumConformers() == 0 or n < 3:
        return
    P = np.array(mol3d.GetConformer(conf_id).GetPositions())
    A = P[:n].copy()
    B = np.array(mol2d.GetConformer().GetPositions())[:n].copy()
    B[:, 2] = 0.0
    ca, cb = A.mean(axis=0), B.mean(axis=0)
    H = (A - ca).T @ (B - cb)
    U, _, Vt = np.linalg.svd(H)
    D = np.eye(3)
    D[2, 2] = np.sign(np.linalg.det(Vt.T @ U.T)) or 1.0
    Rm = Vt.T @ D @ U.T
    out = (P - ca) @ Rm.T
    conf = mol3d.GetConformer(conf_id)
    for i in range(len(out)):
        conf.SetAtomPosition(i, out[i].tolist())


def normalize_scale(mol: Chem.Mol, target=1.5):
    """Bring drawing / foreign coordinates to ~Å (ChemDraw points, 1/65536 units, ...)."""
    if mol.GetNumConformers() == 0 or mol.GetNumBonds() == 0:
        return False
    pos = mol.GetConformer().GetPositions()
    bl = [np.linalg.norm(pos[b.GetBeginAtomIdx()] - pos[b.GetEndAtomIdx()]) for b in mol.GetBonds()
          if mol.GetAtomWithIdx(b.GetBeginAtomIdx()).GetAtomicNum() > 1
          and mol.GetAtomWithIdx(b.GetEndAtomIdx()).GetAtomicNum() > 1] or \
         [np.linalg.norm(pos[b.GetBeginAtomIdx()] - pos[b.GetEndAtomIdx()]) for b in mol.GetBonds()]
    med = float(np.median(bl)) if bl else 0.0
    if med <= 1e-6 or not np.isfinite(med):
        return False
    if 0.9 <= med <= 2.2:
        return False
    conf = mol.GetConformer()
    c = pos.mean(axis=0)
    f = target / med
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, ((pos[i] - c) * f).tolist())
    return True


def from_drawing(mol2d: Chem.Mol, seed=0xF00D) -> Chem.Mol:
    """
    'Paste as drawn': keep the drawing's layout, lift it into 3D just enough that the
    stereochemistry is right (wedges -> depth, H placed on the correct side), add hydrogens.
    Nothing is optimised; the user minimizes it.
    """
    rng = np.random.default_rng(seed)
    m = Chem.Mol(mol2d)
    if m.GetNumConformers() == 0:
        rdDepictor.Compute2DCoords(m)
    normalize_scale(m)
    conf = m.GetConformer()
    pos = conf.GetPositions()
    pos[:, 2] = 0.0
    pos -= pos.mean(axis=0)
    # wedge / hash bonds of the drawing -> depth of the narrow-end's partner
    for b in m.GetBonds():
        d = b.GetBondDir()
        if d in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH):
            j = b.GetEndAtomIdx()
            pos[j, 2] += 0.75 if d == Chem.BondDir.BEGINWEDGE else -0.75
    pos[:, 2] += rng.uniform(-0.04, 0.04, len(pos))   # break perfect planarity
    for i in range(m.GetNumAtoms()):
        conf.SetAtomPosition(i, pos[i].tolist())
    conf.Set3D(True)
    mh = Chem.AddHs(m, addCoords=True)
    _fix_stereo_from_drawing(mol2d, mh)
    _separate_fragments(mh)
    mh.GetConformer().Set3D(True)
    return mh


def _cip_labels(mol):
    out = {}
    for idx, lab in Chem.FindMolChiralCenters(mol, includeUnassigned=False, useLegacyImplementation=False):
        out[idx] = lab
    return out


def _fix_stereo_from_drawing(mol2d, mh, rounds=4):
    """Flip H / substituent depth until each stereocentre matches the drawing."""
    want = _cip_labels(mol2d)
    if not want:
        return
    conf = mh.GetConformer()
    for _ in range(rounds):
        chk = Chem.Mol(mh)
        Chem.AssignStereochemistryFrom3D(chk)
        have = _cip_labels(chk)
        wrong = [c for c, lab in want.items() if have.get(c) != lab]
        if not wrong:
            return
        pos = conf.GetPositions()
        for c in wrong:
            nbrs = [n.GetIdx() for n in mh.GetAtomWithIdx(c).GetNeighbors()]
            hs = [n for n in nbrs if mh.GetAtomWithIdx(n).GetAtomicNum() == 1]
            pick = hs[0] if hs else max(nbrs, key=lambda n: abs(pos[n, 2] - pos[c, 2]))
            z = pos[pick, 2] - pos[c, 2]
            pos[pick, 2] = pos[c, 2] - (z if abs(z) > 0.2 else 0.8 * (1 if z <= 0 else -1))
            conf.SetAtomPosition(pick, pos[pick].tolist())


def _separate_fragments(mol, gap=3.0):
    """ETKDG embeds disconnected fragments (salts, reagents) on top of each other."""
    frags = Chem.GetMolFrags(mol)
    if len(frags) < 2:
        return
    for conf in mol.GetConformers():
        pos = conf.GetPositions()
        x = 0.0
        for idx in frags:
            idx = list(idx)
            p = pos[idx] - pos[idx].mean(axis=0)
            p[:, 0] += x - p[:, 0].min()
            x = p[:, 0].max() + gap
            for j, i in enumerate(idx):
                conf.SetAtomPosition(i, p[j].tolist())


def minimize_current(mol: Chem.Mol, conf_id: int, ff: str, max_iters=2000, fixed=None):
    """Minimize starting from the current (possibly hand-edited) coordinates."""
    msgs = []
    used = choose_ff(mol, ff, msgs)
    if not used:
        return float("nan"), "", msgs
    e, conv = minimize(mol, used, conf_id, max_iters, fixed)
    if fixed:
        msgs.append(f"{len(fixed)} atom(s) held fixed.")
    if not conv:
        msgs.append("Minimization did not fully converge (try again).")
    return e, used, msgs


def stereo_changes(reference: Chem.Mol | None, mol3d: Chem.Mol, conf_id: int) -> str:
    """Compare the stereo of the drawing with the (edited) 3D model."""
    if reference is None:
        return ""
    ref = stereo_smiles(reference)
    chk = Chem.Mol(mol3d, confId=conf_id)
    Chem.AssignStereochemistryFrom3D(chk, confId=chk.GetConformer().GetId())
    now = stereo_smiles(chk)
    if ref == now:
        return ""
    if Chem.CanonSmiles(ref, useChiral=0) != Chem.CanonSmiles(now, useChiral=0):
        return "Connectivity differs from the drawing."
    return ("Stereochemistry now differs from the drawing (an atom was moved through "
            "a stereocentre or C=C). Undo, or keep it if intended.")


def side_of_bond(mol: Chem.Mol, keep: int, move: int) -> list[int]:
    """Atoms on the `move` side of the keep-move bond (empty if the bond is in a ring)."""
    b = mol.GetBondBetweenAtoms(keep, move)
    if b is None or b.IsInRing():
        return []
    rw = Chem.RWMol(mol)
    rw.RemoveBond(keep, move)
    for frag in Chem.GetMolFrags(rw.GetMol(), sanitizeFrags=False):
        if move in frag:
            return list(frag)
    return [move]


def attached_hydrogens(mol: Chem.Mol, atoms) -> set[int]:
    out = set()
    for i in atoms:
        for n in mol.GetAtomWithIdx(int(i)).GetNeighbors():
            if n.GetAtomicNum() == 1:
                out.add(n.GetIdx())
    return out


def set_internal(mol: Chem.Mol, conf_id: int, atoms: list[int], value: float):
    """Set a distance (2 atoms), angle (3) or dihedral (4), moving the side of the last atom."""
    from rdkit.Chem import rdMolTransforms as T
    conf = mol.GetConformer(conf_id)
    if len(atoms) == 2:
        i, j = atoms
        side = side_of_bond(mol, i, j)
        if not side:  # not bonded / ring bond: move only the last atom (+ its H)
            side = [j] + sorted(attached_hydrogens(mol, [j]))
        pos = conf.GetPositions()
        d = pos[j] - pos[i]
        n = np.linalg.norm(d)
        if n < 1e-6:
            raise ValueError("Atoms coincide.")
        shift = d / n * (value - n)
        for k in side:
            conf.SetAtomPosition(k, (pos[k] + shift).tolist())
    elif len(atoms) == 3:
        T.SetAngleDeg(conf, *atoms, value)
    elif len(atoms) == 4:
        T.SetDihedralDeg(conf, *atoms, value)
    else:
        raise ValueError("Select 2, 3 or 4 atoms.")


def single_point(mol, conf_id, ff):
    if not ff or not _ff_available(mol, ff):
        return float("nan")
    f = _make_ff(mol, ff, conf_id)
    return f.CalcEnergy() if f else float("nan")


# --------------------------------------------------------------------------------------
# 3D -> 2D drawing (what Chem3D calls "copy to ChemDraw")
# --------------------------------------------------------------------------------------
@dataclass
class Drawing2D:
    mol: Chem.Mol
    molblock: str
    smiles: str
    cdx: bytes | None
    warnings: list
    depth: np.ndarray | None = None   # per atom of `mol`, same scale as x/y, + = toward viewer


def _single_conf(mol, conf_id):
    m = Chem.Mol(mol, confId=conf_id) if conf_id >= 0 else Chem.Mol(mol)
    return m


def _projection_wedges(m: Chem.Mol, z: np.ndarray):
    """Wedge (toward viewer) / hash (away) per stereocentre, chosen from real depth."""
    rw = Chem.RWMol(m)
    for b in rw.GetBonds():
        if b.GetBondType() == Chem.BondType.SINGLE:
            b.SetBondDir(Chem.BondDir.NONE)
    centers = [a.GetIdx() for a in rw.GetAtoms()
               if a.GetChiralTag() in (Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                                       Chem.ChiralType.CHI_TETRAHEDRAL_CCW)]
    cset = set(centers)
    used = set()
    pos = rw.GetConformer().GetPositions()
    plan = []
    for c in centers:
        best = None
        for b in rw.GetAtomWithIdx(c).GetBonds():
            if b.GetBondType() != Chem.BondType.SINGLE or b.GetIdx() in used:
                continue
            n = b.GetOtherAtomIdx(c)
            ln = np.linalg.norm(pos[n, :2] - pos[c, :2]) + 1e-6
            dz = z[n] - z[c]
            score = abs(dz) / (ln + abs(dz))
            score += 0.35 if n not in cset else 0.0
            score += 0.15 if not b.IsInRing() else 0.0
            score += 0.10 if rw.GetAtomWithIdx(n).GetAtomicNum() == 1 else 0.0
            if best is None or score > best[0]:
                best = (score, b.GetIdx(), n, dz)
        if best:
            used.add(best[1])
            plan.append((c, best[2], best[3]))
    for c, n, dz in plan:
        b = rw.GetBondBetweenAtoms(c, n)
        if b.GetBeginAtomIdx() != c:  # wedge must start at the stereocentre
            arom = b.GetIsAromatic()
            rw.RemoveBond(c, n)
            rw.AddBond(c, n, Chem.BondType.SINGLE)
            b = rw.GetBondBetweenAtoms(c, n)
            b.SetIsAromatic(arom)
        b.SetBondDir(Chem.BondDir.BEGINWEDGE if dz > 0 else Chem.BondDir.BEGINDASH)
    return rw.GetMol()


def _roundtrip_ok(molblock, ref):
    back = Chem.MolFromMolBlock(molblock)
    return back is not None and Chem.MolToSmiles(back) == ref


def drawing_from_3d(mol3d: Chem.Mol, conf_id: int, rotation: np.ndarray,
                    mode: str = "view", keep_h: bool = False,
                    bond_length: float = 1.5, cdx_bond_pt: float = 14.4) -> Drawing2D:
    """
    mode="view"  : orthographic projection of the current 3D view (Chem3D-style)
    mode="clean" : a standard clean 2D depiction (stereo taken from the 3D model)
    """
    warnings = []
    depth = None
    m = _single_conf(mol3d, conf_id)
    Chem.AssignStereochemistryFrom3D(m, confId=m.GetConformer().GetId())
    if not keep_h:
        m = Chem.RemoveHs(m)
    ref = Chem.MolToSmiles(Chem.RemoveHs(m))
    conf = m.GetConformer()

    if mode == "clean":
        rdDepictor.SetPreferCoordGen(True)
        m2 = Chem.Mol(m)
        m2.RemoveAllConformers()
        rdDepictor.Compute2DCoords(m2)
        rdDepictor.StraightenDepiction(m2)
        Chem.WedgeMolBonds(m2, m2.GetConformer())
        m = m2
    else:
        pos = conf.GetPositions()
        heavy = [a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum() > 1] or list(range(len(pos)))
        center = pos[heavy].mean(axis=0)
        v = (pos - center) @ np.asarray(rotation).T
        # scale so that the *3D* bond length maps to the requested 2D bond length
        bl = [np.linalg.norm(pos[b.GetBeginAtomIdx()] - pos[b.GetEndAtomIdx()])
              for b in m.GetBonds()
              if m.GetAtomWithIdx(b.GetBeginAtomIdx()).GetAtomicNum() > 1
              and m.GetAtomWithIdx(b.GetEndAtomIdx()).GetAtomicNum() > 1]
        scale = bond_length / (float(np.median(bl)) if bl else 1.5)
        v *= scale
        for i in range(m.GetNumAtoms()):
            conf.SetAtomPosition(i, (float(v[i, 0]), float(v[i, 1]), 0.0))
        conf.Set3D(False)
        depth = v[:, 2].copy()
        m = _projection_wedges(m, v[:, 2])
        # overlapping atoms make a bad drawing
        xy = v[:, :2]
        if len(xy) > 1:
            d = np.linalg.norm(xy[:, None, :] - xy[None, :, :], axis=-1)
            np.fill_diagonal(d, 1e9)
            if d.min() < 0.25 * bond_length:
                warnings.append("Some atoms overlap in this orientation - rotate the model, "
                                "use 'Best 2D orientation', or copy as clean 2D.")

    block = Chem.MolToMolBlock(m)
    if not _roundtrip_ok(block, ref):
        if mode == "view":
            m3 = Chem.Mol(m)
            for b in m3.GetBonds():
                b.SetBondDir(Chem.BondDir.NONE)
            Chem.WedgeMolBonds(m3, m3.GetConformer())
            b3 = Chem.MolToMolBlock(m3)
            if _roundtrip_ok(b3, ref):
                m, block = m3, b3
        if not _roundtrip_ok(block, ref):
            warnings.append("Stereo (likely a C=C geometry) is ambiguous in this projection; "
                            "check the drawing or use clean 2D.")
    cdx = None
    if HAS_CDX:
        try:
            cdx = chemdraw_cdx(m, cdx_bond_pt, None if mode == "clean" else bond_length)
        except Exception:  # noqa: BLE001
            cdx = None
    return Drawing2D(m, block, ref, cdx, warnings, depth)


def best_view_rotation(mol: Chem.Mol, conf_id: int) -> np.ndarray:
    """Rotation that puts the two largest principal axes in the screen plane."""
    pos = mol.GetConformer(conf_id).GetPositions()
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1] or list(range(len(pos)))
    p = pos[heavy] - pos[heavy].mean(axis=0)
    if len(p) < 3:
        return np.eye(3)
    _, _, vt = np.linalg.svd(p, full_matrices=False)
    r = vt.copy()
    if np.linalg.det(r) < 0:
        r[2] *= -1
    return r


# --------------------------------------------------------------------------------------
# Exports
# --------------------------------------------------------------------------------------
def default_charge_mult(mol: Chem.Mol):
    charge = Chem.GetFormalCharge(mol)
    rad = Descriptors.NumRadicalElectrons(mol)
    mult = rad + 1
    nelec = sum(a.GetAtomicNum() for a in mol.GetAtoms()) - charge
    if (nelec % 2 == 0) != (mult % 2 == 1):  # parity fix
        mult = 2 if nelec % 2 else 1
    return charge, mult


def xyz_lines(mol, conf_id):
    conf = mol.GetConformer(conf_id)
    out = []
    for a in mol.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        out.append(f"  {a.GetSymbol():<2s} {p.x:14.8f} {p.y:14.8f} {p.z:14.8f}")
    return out


def xyz_block(mol, conf_id, comment=""):
    lines = xyz_lines(mol, conf_id)
    return f"{len(lines)}\n{comment}\n" + "\n".join(lines) + "\n"


ORCA_PRESETS = [
    ("r2SCAN-3c  Opt + Freq (fast, robust)", "! r2SCAN-3c Opt Freq TightSCF"),
    ("B3LYP-D3BJ / def2-SVP  Opt", "! B3LYP D3BJ def2-SVP Opt TightSCF"),
    ("PBE0-D4 / def2-TZVP  Opt + Freq", "! PBE0 D4 def2-TZVP Opt Freq TightSCF"),
    ("wB97X-D3 / def2-TZVP  Single point", "! wB97X-D3 def2-TZVP TightSCF"),
    ("GFN2-xTB  Opt (pre-optimisation)", "! XTB2 Opt"),
    ("DLPNO-CCSD(T) / cc-pVTZ  Single point", "! DLPNO-CCSD(T) cc-pVTZ cc-pVTZ/C TightSCF"),
    ("Custom", "! "),
]


def orca_input(mol, conf_id, keywords="! r2SCAN-3c Opt Freq TightSCF", charge=0, mult=1,
               nprocs=4, maxcore=2000, extra="", title="", xyzfile=None):
    lines = []
    if title:
        lines.append(f"# {title}")
    lines.append(f"# formula: {rdMolDescriptors.CalcMolFormula(mol)}   generated by Mol3D Bridge")
    lines.append(keywords.strip())
    if nprocs and nprocs > 1:
        lines.append(f"%pal nprocs {nprocs} end")
    if maxcore:
        lines.append(f"%maxcore {maxcore}")
    if extra.strip():
        lines.append(extra.strip())
    lines.append("")
    if xyzfile:
        lines.append(f"* xyzfile {charge} {mult} {xyzfile}")
    else:
        lines.append(f"* xyz {charge} {mult}")
        lines.extend(xyz_lines(mol, conf_id))
        lines.append("*")
    return "\n".join(lines) + "\n"


def molblock_3d(mol, conf_id):
    return Chem.MolToMolBlock(mol, confId=conf_id)


def write_sdf(path, mol, energies=None, ff=""):
    w = Chem.SDWriter(path)
    for i, c in enumerate(mol.GetConformers()):
        if energies and i < len(energies):
            mol.SetProp("Energy_kcal_mol", f"{energies[i]:.4f}")
            mol.SetProp("ForceField", ff)
        w.write(mol, confId=c.GetId())
    w.close()


def formula_and_mass(mol):
    return rdMolDescriptors.CalcMolFormula(mol), Descriptors.MolWt(mol)


# ======================================================================================
# Step-by-step optimisation (trajectory for the "falling into the minimum" animation)
# ======================================================================================
KCAL_PER_EH = 627.509474
XTB_METHODS = ("GFN2-xTB", "GFN-FF")
METHODS = FORCE_FIELDS + XTB_METHODS


def is_xtb(method: str) -> bool:
    return method in XTB_METHODS


@dataclass
class OptResult:
    method: str
    frames: list            # [(positions Nx3, energy kcal/mol, gradient-RMS or None)]
    final: np.ndarray       # converged positions
    energy: float           # kcal/mol (xtb: total energy converted from Eh)
    converged: bool = True
    messages: list = field(default_factory=list)
    cancelled: bool = False


class Cancelled(Exception):
    pass


def ff_trajectory(mol: Chem.Mol, conf_id: int, ff: str, fixed=None, max_step=0.02,
                  max_iter=4000, grms_tol=0.05, cancel=None) -> OptResult:
    """
    FIRE minimisation with a small per-atom step cap (default 0.02 Å per step), so every
    recorded frame is a small, smooth move; finished by a BFGS polish to the true minimum.
    """
    msgs = []
    used = choose_ff(mol, ff, msgs)
    if not used:
        raise RuntimeError("No force-field parameters for this molecule. Try GFN2-xTB or GFN-FF.")
    work = Chem.Mol(mol)
    field_ = _make_ff(work, used, conf_id)
    x = np.array(work.GetConformer(conf_id).GetPositions(), dtype=float)
    n = len(x)
    mask = np.ones(n)
    if fixed:
        mask[list(fixed)] = 0.0
    nfree = max(1.0, mask.sum() * 3)
    v = np.zeros_like(x)
    dt, dt_max, alpha, npos = 0.02, 0.12, 0.1, 0
    frames = [(x.copy(), field_.CalcEnergy(x.ravel().tolist()), None)]
    for it in range(max_iter):
        if cancel is not None and cancel.is_set():
            break
        g = np.array(field_.CalcGrad(x.ravel().tolist())).reshape(-1, 3) * mask[:, None]
        grms = float(np.sqrt((g ** 2).sum() / nfree))
        if grms < grms_tol:
            break
        F = -g
        P = float((F * v).sum())
        if P > 0:
            v = (1 - alpha) * v + alpha * np.linalg.norm(v) * F / (np.linalg.norm(F) + 1e-12)
            npos += 1
            if npos > 5:
                dt, alpha = min(dt * 1.1, dt_max), alpha * 0.99
        else:
            v[:] = 0.0
            dt, alpha, npos = dt * 0.5, 0.1, 0
        v += F * dt
        dx = v * dt
        norms = np.linalg.norm(dx, axis=1)
        dx *= np.minimum(1.0, max_step / (norms + 1e-12))[:, None]
        x += dx * mask[:, None]
        frames.append((x.copy(), field_.CalcEnergy(x.ravel().tolist()), grms))

    cancelled = cancel is not None and cancel.is_set()
    conf = work.GetConformer(conf_id)
    for i in range(n):
        conf.SetAtomPosition(i, x[i].tolist())
    conv = True
    if not cancelled:  # BFGS polish -> exactly the same minimum as a normal minimisation
        e, conv = minimize(work, used, conf_id, 2000, fixed)
        fin = np.array(conf.GetPositions())
        frames.append((fin.copy(), e, 0.0))
    else:
        fin, e = x.copy(), frames[-1][1]
        msgs.append("Stopped by user.")
    if fixed:
        msgs.append(f"{len(fixed)} atom(s) held fixed.")
    if not conv:
        msgs.append("Not fully converged - minimize again.")
    return OptResult(used, frames, fin, e, conv, msgs, cancelled)


# --------------------------------------------------------------------------- xtb
XTB_VERSION = "6.7.1"
XTB_URLS = {
    "win": f"https://github.com/grimme-lab/xtb/releases/download/v{XTB_VERSION}/xtb-{XTB_VERSION}pre-windows-x86_64.zip",
    "linux": f"https://github.com/grimme-lab/xtb/releases/download/v{XTB_VERSION}/xtb-{XTB_VERSION}-linux-x86_64.tar.xz",
}
ALPB_SOLVENTS = ("none", "water", "methanol", "acetonitrile", "dmso", "dmf", "acetone", "chcl3",
                 "ch2cl2", "thf", "toluene", "benzene", "ether", "ethylacetate", "hexane", "octanol")


@dataclass
class XtbSettings:
    path: str = ""
    solvent: str = "none"
    level: str = "normal"      # crude sloppy loose normal tight vtight
    threads: int = 4
    charge: int | None = None  # None = from the structure
    uhf: int | None = None


def _exe_name():
    import sys
    return "xtb.exe" if sys.platform.startswith("win") else "xtb"


def find_xtb(hint: str = "", extra_dirs=()) -> str | None:
    import os, shutil
    if hint and os.path.isfile(hint):
        return hint
    w = shutil.which("xtb")
    if w:
        return w
    for d in extra_dirs:
        for root, _, files in os.walk(d):
            if _exe_name() in files and os.path.basename(root) == "bin":
                return os.path.join(root, _exe_name())
    return None


def download_xtb(dest_dir: str, progress=None) -> str:
    """Download the official xtb release (Windows/Linux x86-64), verify sha256, unpack."""
    import hashlib, os, sys, tarfile, urllib.request, zipfile
    key = "win" if sys.platform.startswith("win") else "linux" if sys.platform.startswith("linux") else None
    if key is None:
        raise RuntimeError("No official macOS binary. Install with:  conda install -c conda-forge xtb  "
                           "(or  brew install grimme-lab/qc/xtb), then pick the xtb file.")
    url = XTB_URLS[key]
    os.makedirs(dest_dir, exist_ok=True)
    arc = os.path.join(dest_dir, os.path.basename(url))

    def hook(blocks, bs, total):
        if progress and total > 0:
            progress(f"Downloading xtb… {min(100, blocks * bs * 100 // total)} %")

    urllib.request.urlretrieve(url, arc, hook)
    try:
        with urllib.request.urlopen(url + ".sha256", timeout=30) as r:
            want = r.read().decode().split()[0].strip().lower()
        h = hashlib.sha256()
        with open(arc, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        if h.hexdigest() != want:
            raise RuntimeError("Checksum mismatch – download corrupted, try again.")
    except RuntimeError:
        raise
    except Exception:
        pass  # checksum file unavailable: continue
    if progress:
        progress("Unpacking xtb…")
    if arc.endswith(".zip"):
        with zipfile.ZipFile(arc) as z:
            z.extractall(dest_dir)
    else:
        with tarfile.open(arc) as t:
            try:
                t.extractall(dest_dir, filter="data")
            except TypeError:  # Python < 3.12
                t.extractall(dest_dir)
    os.remove(arc)
    exe = find_xtb("", [dest_dir])
    if not exe:
        raise RuntimeError("Downloaded xtb but could not find the executable.")
    if not sys.platform.startswith("win"):
        os.chmod(exe, 0o755)
    return exe


def _parse_multixyz(text: str, natoms: int):
    """Complete frames of an xtb trajectory log -> [(positions, energy_Eh, gnorm)]."""
    lines = text.splitlines()
    out, i = [], 0
    while i + natoms + 2 <= len(lines):
        try:
            n = int(lines[i].split()[0])
        except (ValueError, IndexError):
            break
        if n != natoms:
            break
        com = lines[i + 1].split()
        e = g = float("nan")
        for k, tok in enumerate(com):
            if tok.startswith("energy") and k + 1 < len(com):
                e = float(com[k + 1])
            if tok.startswith("gnorm") and k + 1 < len(com):
                g = float(com[k + 1])
        pos = np.array([[float(t) for t in ln.split()[1:4]] for ln in lines[i + 2:i + 2 + natoms]])
        out.append((pos, e, g))
        i += natoms + 2
    return out


def xtb_optimize(mol: Chem.Mol, conf_id: int, method: str, st: XtbSettings, fixed=None,
                 progress=None, cancel=None, keep_dir=None) -> OptResult:
    """Run `xtb --opt` and stream every optimisation step (xtbopt.log) through `progress`."""
    import os, subprocess, sys, tempfile, time
    exe = find_xtb(st.path)
    if exe and not sys.platform.startswith("win") and not os.access(exe, os.X_OK):
        try:
            os.chmod(exe, 0o755)
        except OSError:
            pass
    if not exe:
        raise RuntimeError("xtb executable not found. Calculate ▸ xtb settings… lets you download "
                           "it or point to an installed xtb.")
    ch, mu = default_charge_mult(mol)
    charge = ch if st.charge is None else st.charge
    uhf = (mu - 1) if st.uhf is None else st.uhf
    natoms = mol.GetNumAtoms()
    wd = keep_dir or tempfile.mkdtemp(prefix="m3b_xtb_")
    with open(os.path.join(wd, "start.xyz"), "w") as fh:
        fh.write(xyz_block(mol, conf_id, "Mol3D Bridge"))
    cmd = [exe, "start.xyz", "--opt", st.level, "--chrg", str(charge), "--uhf", str(uhf)]
    cmd += ["--gfnff"] if method == "GFN-FF" else ["--gfn", "2"]
    if st.solvent and st.solvent != "none":
        cmd += ["--alpb", st.solvent]
    if fixed:
        with open(os.path.join(wd, "xcontrol"), "w") as fh:
            fh.write("$fix\n   atoms: " + ",".join(str(i + 1) for i in sorted(fixed)) + "\n$end\n")
        cmd += ["--input", "xcontrol"]
    env = dict(os.environ)
    env.update(OMP_NUM_THREADS=f"{max(1, st.threads)},1", MKL_NUM_THREADS=str(max(1, st.threads)),
               OMP_STACKSIZE="1G", OMP_MAX_ACTIVE_LEVELS="1")
    kw = {}
    if sys.platform.startswith("win"):
        kw["creationflags"] = 0x08000000  # CREATE_NO_WINDOW
    out = open(os.path.join(wd, "xtb.out"), "w")
    proc = subprocess.Popen(cmd, cwd=wd, stdout=out, stderr=subprocess.STDOUT, env=env, **kw)
    log = os.path.join(wd, "xtbopt.log")
    frames, sent = [], 0
    cancelled = False
    while True:
        rc = proc.poll()
        if cancel is not None and cancel.is_set():
            proc.kill()
            cancelled = True
            rc = proc.wait()
        if os.path.exists(log):
            with open(log, errors="replace") as fh:
                frames = _parse_multixyz(fh.read(), natoms)
            while sent < len(frames):
                p, e, g = frames[sent]
                if progress:
                    progress((p, e * KCAL_PER_EH, g))
                sent += 1
        if rc is not None:
            break
        time.sleep(0.05)
    out.close()
    with open(os.path.join(wd, "xtb.out"), errors="replace") as fh:
        text = fh.read()
    frames_k = [(p, e * KCAL_PER_EH, g) for p, e, g in frames]
    msgs = []
    if cancelled:
        if not frames_k:
            raise Cancelled()
        return OptResult(method, frames_k, frames_k[-1][0], frames_k[-1][1], False, ["Stopped by user."], True)
    final_file = None
    for name in ("xtbopt.xyz", "xtblast.xyz"):
        if os.path.exists(os.path.join(wd, name)):
            final_file = os.path.join(wd, name)
            break
    if rc != 0 and final_file is None:
        tail = "\n".join(l for l in text.splitlines()[-25:] if l.strip())
        raise RuntimeError(f"xtb failed (exit {rc}).\n\n{tail}")
    with open(final_file) as fh:
        fin = _parse_multixyz(fh.read(), natoms)
    pos, e, g = fin[-1]
    conv = "GEOMETRY OPTIMIZATION CONVERGED" in text
    if not conv:
        msgs.append("xtb optimisation did not converge.")
    if st.solvent and st.solvent != "none":
        msgs.append(f"ALPB({st.solvent}).")
    if fixed:
        msgs.append(f"{len(fixed)} atom(s) held fixed.")
    msgs.append(f"charge {charge}, unpaired {uhf}.")
    frames_k.append((pos, e * KCAL_PER_EH, g))
    return OptResult(method, frames_k, pos, e * KCAL_PER_EH, conv, msgs)


def apply_positions(mol: Chem.Mol, conf_id: int, pos: np.ndarray):
    conf = mol.GetConformer(conf_id)
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, pos[i].tolist())


def format_energy(e_kcal: float, method: str) -> str:
    if e_kcal is None or math.isnan(e_kcal):
        return "–"
    if is_xtb(method):
        return f"{e_kcal / KCAL_PER_EH:.6f} Eh"
    return f"{e_kcal:.3f} kcal/mol"


# --------------------------------------------------------------------------- Gaussian
GAUSSIAN_PRESETS = [
    ("B3LYP-D3BJ / 6-31G(d)  Opt + Freq", "#p opt freq b3lyp/6-31g(d) empiricaldispersion=gd3bj"),
    ("B3LYP-D3BJ / def2-SVP  Opt + Freq", "#p opt freq b3lyp/def2svp empiricaldispersion=gd3bj"),
    ("ωB97X-D / def2-TZVP  Opt + Freq", "#p opt freq wb97xd/def2tzvp"),
    ("M06-2X / 6-311+G(d,p)  Opt + Freq", "#p opt freq m062x/6-311+g(d,p) int=ultrafine"),
    ("PBE0-D3BJ / def2-SVP  Opt, SMD(water)",
     "#p opt pbe1pbe/def2svp empiricaldispersion=gd3bj scrf=(smd,solvent=water)"),
    ("B3LYP / 6-31G(d)  TS search (Opt=TS)", "#p opt=(ts,calcfc,noeigentest) freq b3lyp/6-31g(d)"),
    ("Custom", "#p "),
]


def gaussian_input(mol, conf_id, route="#p opt freq b3lyp/6-31g(d) empiricaldispersion=gd3bj",
                   charge=0, mult=1, nprocs=4, mem="4GB", chk="", title="", extra=""):
    lines = []
    if chk:
        lines.append(f"%chk={chk}")
    if nprocs and nprocs > 1:
        lines.append(f"%nprocshared={nprocs}")
    if mem:
        lines.append(f"%mem={mem}")
    lines.append(route.strip())
    lines.append("")
    lines.append(title.strip() or f"{rdMolDescriptors.CalcMolFormula(mol)} generated by Mol3D Bridge")
    lines.append("")
    lines.append(f"{charge} {mult}")
    conf = mol.GetConformer(conf_id)
    for a in mol.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        lines.append(f" {a.GetSymbol():<2s} {p.x:14.8f} {p.y:14.8f} {p.z:14.8f}")
    lines.append("")
    if extra.strip():
        lines.append(extra.strip())
        lines.append("")
    return "\n".join(lines) + "\n"


# ======================================================================================
# 3D-preserving copy (Chem3D-style "copy to ChemDraw" that can still be rotated in 3D)
# ======================================================================================
# ChemDraw stores 3D atoms as a CDX node property kCDXProp_3DPosition (0x0201, CDXML "xyz").
# CDX specification (CambridgeSoft/Revvity): CDXCoordinate = INT32 in 1/65536 pt; a CDXPoint2D
# is stored y, x and a CDXPoint3D is stored **z, y, x**; the 3D frame is right-handed with the
# 2D axes (x right, y down), so +z points *into* the screen; x/y of the 3D position should equal
# the 2D position.  (CDXML text writes "x y z".)
CDX_NODE, CDX_P2D, CDX_P3D = 0x8004, 0x0200, 0x0201


def _cdx_header_len(buf: bytes) -> int:
    if not buf.startswith(b"VjCD0100"):
        raise ValueError("not a CDX document")
    for off in (22, 28):
        if len(buf) > off + 1 and int.from_bytes(buf[off:off + 2], "little") == 0x8000:
            return off
    raise ValueError("unrecognised CDX header")


def _cdx_walk(buf: bytes, node_extra=None, drop_3d=False, found3d=None):
    """Walk CDX objects. Returns (node 2D positions in order, rewritten bytes)."""
    import struct
    i = _cdx_header_len(buf)
    out = bytearray(buf[:i])
    stack, nodes = [], []
    while i < len(buf):
        tag = int.from_bytes(buf[i:i + 2], "little")
        if tag == 0:
            out += buf[i:i + 2]
            i += 2
            if not stack:
                break
            stack.pop()
            if not stack:
                out += buf[i:]  # trailing bytes after the document
                break
            continue
        if tag & 0x8000:
            out += buf[i:i + 6]
            i += 6
            stack.append(tag)
            if tag == CDX_NODE:
                nodes.append(None)
                if node_extra is not None:
                    out += node_extra(len(nodes) - 1)
            continue
        ln = int.from_bytes(buf[i + 2:i + 4], "little")
        hl = 4
        if ln == 0xFFFF:
            ln = int.from_bytes(buf[i + 4:i + 8], "little")
            hl = 8
        data = buf[i + hl:i + hl + ln]
        if stack and stack[-1] == CDX_NODE:
            if tag == CDX_P2D and len(data) == 8:
                y, x = struct.unpack("<ii", data)
                nodes[-1] = (x / 65536.0, y / 65536.0)
            if tag == CDX_P3D:
                if found3d is not None and len(data) == 12:
                    zz, yy, xx = struct.unpack("<iii", data)      # spec order: z, y, x
                    found3d[len(nodes) - 1] = (xx / 65536.0, yy / 65536.0, zz / 65536.0)
                if node_extra is not None or drop_3d:
                    i += hl + ln      # drop an existing 3D position (ours replaces it)
                    continue
        out += buf[i:i + hl + ln]
        i += hl + ln
    return nodes, bytes(out)


CDX_BONDLENGTH = 0x0805          # document default bond length (CDXCoordinate, 1/65536 pt)
CDX_COORD_PROPS = {0x0200: 2, 0x0201: 3, 0x0204: 4}   # 2D point, 3D point, rectangle


def cdx_set_scale(cdx: bytes, factor: float, bond_length_pt: float) -> bytes:
    """
    Scale every coordinate of a CDX document by `factor` and set the document's default
    bond length.  RDKit writes BondLength = 14 *raw units* (0.0002 pt); ChemDraw then reports
    'different settings' and, if the user keeps the document's settings, scales the pasted
    structure up ~67,000x ("expand the document to 9212 pages").
    """
    import struct
    i = _cdx_header_len(cdx)
    out = bytearray(cdx[:i])
    stack, wrote_bl = [], False
    while i < len(cdx):
        tag = int.from_bytes(cdx[i:i + 2], "little")
        if tag == 0:
            if len(stack) == 1 and not wrote_bl:     # document had no BondLength: add one
                out += struct.pack("<HHi", CDX_BONDLENGTH, 4, int(round(bond_length_pt * 65536)))
                wrote_bl = True
            out += cdx[i:i + 2]
            i += 2
            if not stack:
                break
            stack.pop()
            if not stack:
                out += cdx[i:]
                break
            continue
        if tag & 0x8000:
            out += cdx[i:i + 6]
            i += 6
            stack.append(tag)
            continue
        ln = int.from_bytes(cdx[i + 2:i + 4], "little")
        hl = 4
        if ln == 0xFFFF:
            ln = int.from_bytes(cdx[i + 4:i + 8], "little")
            hl = 8
        data = cdx[i + hl:i + hl + ln]
        if tag == CDX_BONDLENGTH and len(stack) == 1 and ln == 4:
            data = struct.pack("<i", int(round(bond_length_pt * 65536)))
            wrote_bl = True
        elif tag in CDX_COORD_PROPS and ln == 4 * CDX_COORD_PROPS[tag]:
            vals = struct.unpack(f"<{ln // 4}i", data)
            data = struct.pack(f"<{ln // 4}i", *(int(round(v * factor)) for v in vals))
        out += cdx[i:i + hl] + data
        i += hl + ln
    return bytes(out)


def _cdx_fit_scale(cdx: bytes, xy: np.ndarray) -> float:
    """pt per drawing unit of an RDKit CDX (p = s*(x, -y) + t); also checks the atom order."""
    nodes, _ = _cdx_walk(cdx)
    if len(nodes) != len(xy) or any(n is None for n in nodes):
        raise ValueError("CDX nodes do not match the atoms")
    P = np.array(nodes)
    A = np.zeros((2 * len(xy), 3))
    A[0::2, 0], A[0::2, 1] = xy[:, 0], 1.0
    A[1::2, 0], A[1::2, 2] = -xy[:, 1], 1.0
    sol, *_ = np.linalg.lstsq(A, P.reshape(-1), rcond=None)
    if not sol[0] > 0 or np.abs(A @ sol - P.reshape(-1)).max() > 0.05 * max(1.0, sol[0]):
        raise ValueError("CDX atom order could not be matched")
    return float(sol[0])


def chemdraw_cdx(mol2d: Chem.Mol, bond_length_pt: float = 14.4, unit_bond: float | None = None) -> bytes | None:
    """RDKit CDX, rescaled so one bond (unit_bond drawing units; default the median 2D bond) is bond_length_pt."""
    if not HAS_CDX:
        return None
    if unit_bond is None:
        pos = mol2d.GetConformer().GetPositions()
        bl = [np.linalg.norm(pos[b.GetBeginAtomIdx(), :2] - pos[b.GetEndAtomIdx(), :2]) for b in mol2d.GetBonds()]
        unit_bond = float(np.median(bl)) if bl else 1.5
        unit_bond = unit_bond if unit_bond > 1e-6 else 1.5
    cdx = Chem.MolToCDXMLBlock(mol2d, Chem.CDXMLFormat.CDX)
    xy = mol2d.GetConformer().GetPositions()[:, :2]
    if mol2d.GetNumAtoms() >= 2:
        try:
            s_ = _cdx_fit_scale(cdx, xy)
            return cdx_set_scale(cdx, bond_length_pt / (unit_bond * s_), bond_length_pt)
        except ValueError:
            pass
    return cdx_set_scale(cdx, 1.0, bond_length_pt)


def strip_cdx_3d(data):
    """Remove 3D positions so readers use the 2D drawing + wedges (RDKit misreads CDX 3D)."""
    import re
    if isinstance(data, (bytes, bytearray)) and bytes(data[:8]) == b"VjCD0100":
        try:
            return _cdx_walk(bytes(data), drop_3d=True)[1]
        except Exception:  # noqa: BLE001
            return data
    if isinstance(data, str):
        return re.sub(r'\sxyz="[^"]*"', "", data)
    return data


def cdx_3d_positions(cdx: bytes):
    """Node 3D positions per the CDX spec, as right-handed y-up coordinates (pt)."""
    found = {}
    nodes, _ = _cdx_walk(cdx, found3d=found)
    return [None if k not in found else (found[k][0], -found[k][1], -found[k][2])
            for k in range(len(nodes))]


def cdx_with_depth(cdx: bytes, xy: np.ndarray, depth: np.ndarray, z_sign: float = 1.0) -> bytes:
    """
    Add a 3D position to every atom node of an RDKit-written CDX (atoms in mol order).
    depth: + = toward the viewer.  Written per spec: order z, y, x; z into the screen.
    z_sign = -1 inverts the depth (Settings, in case a ChemDraw version disagrees).
    """
    import struct
    s_ = _cdx_fit_scale(cdx, xy)
    P = np.array(_cdx_walk(cdx)[0])
    zc = -z_sign * s_ * depth            # toward viewer -> negative z (right-handed, y down)

    def extra(k):
        x, y, z = (int(round(c * 65536.0)) for c in (P[k, 0], P[k, 1], zc[k]))
        return struct.pack("<HH", CDX_P3D, 12) + struct.pack("<iii", z, y, x)

    _, out = _cdx_walk(cdx, extra)
    return out


@dataclass
class Drawing3D:
    smiles: str
    molblock3d: str          # x/y = current view, z = depth (toward viewer), "3D" molfile
    cdx3d: bytes | None      # ChemDraw CDX: 2D positions + wedges + 3D positions
    warnings: list
    molblock3d_wedged: str = ""   # same + wedge/hash flags (for ChemDraw's 2D display)


def drawing_3d_from_view(mol3d: Chem.Mol, conf_id: int, rotation: np.ndarray,
                         keep_h: bool = False, z_sign: float = 1.0,
                         cdx_bond_pt: float = 14.4) -> Drawing3D:
    d = drawing_from_3d(mol3d, conf_id, rotation, "view", keep_h, cdx_bond_pt=cdx_bond_pt)
    warnings = [w for w in d.warnings if "overlap" in w]
    m = Chem.Mol(d.mol)
    conf = m.GetConformer()
    xy = conf.GetPositions()[:, :2]
    z = d.depth
    for i in range(m.GetNumAtoms()):
        conf.SetAtomPosition(i, (float(xy[i, 0]), float(xy[i, 1]), float(z[i])))
    conf.Set3D(True)
    block_w = Chem.MolToMolBlock(m)          # 3D coordinates + the projection's wedges
    for b in m.GetBonds():
        if b.GetBondDir() in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH):
            b.SetBondDir(Chem.BondDir.NONE)
    block = Chem.MolToMolBlock(m)
    back_w = Chem.MolFromMolBlock(block_w)
    if back_w is None or Chem.MolToSmiles(back_w) != d.smiles:
        warnings.append("Wedged 3D molfile check failed.")
    back = Chem.MolFromMolBlock(block)
    if back is None or Chem.MolToSmiles(back) != d.smiles:
        warnings.append("3D molfile stereo check failed.")
    cdx3d = None
    if d.cdx:
        try:
            cdx3d = cdx_with_depth(d.cdx, xy, z, z_sign)
            # self-check: decode our own 3D positions the way the spec defines them
            p3 = cdx_3d_positions(cdx3d)
            chk = Chem.Mol(m)
            cc = chk.GetConformer()
            for i, q in enumerate(p3):
                cc.SetAtomPosition(i, q)
            Chem.AssignStereochemistryFrom3D(chk)
            if z_sign > 0 and stereo_smiles(chk) != d.smiles:
                warnings.append("CDX 3D stereo check failed; copied without 3D positions.")
                cdx3d = d.cdx
            if Chem.MolToSmiles(_combine(Chem.MolsFromCDXML(strip_cdx_3d(cdx3d)))) != d.smiles:
                warnings.append("Wedges may be ambiguous in this orientation - rotate slightly "
                                "if ChemDraw shows a different stereocentre.")
        except Exception as ex:  # noqa: BLE001
            warnings.append(f"Could not add 3D positions to CDX ({ex}).")
            cdx3d = d.cdx
    return Drawing3D(d.smiles, block, cdx3d, warnings, block_w)
