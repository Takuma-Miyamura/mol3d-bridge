"""
Transition-state workflow for pysisyphus (no GUI code).

* atom mapping between reactant and product (MCS + geometry + Hungarian assignment)
* manual relabelling helpers
* YAML input generation, running `pysis`, finding the TS and the two IRC minima

pysisyphus (GPL-3.0, run as a separate program): J. Steinmetzer, S. Kupfer, S. Gräfe,
"pysisyphus: Exploring potential energy surfaces in ground and excited states",
Int. J. Quantum Chem. 2021, 121, e26390, doi:10.1002/qua.26390.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFMCS

from . import core


# ======================================================================================
# assignment (Hungarian algorithm; no SciPy needed)
# ======================================================================================
def hungarian(cost: np.ndarray):
    """Minimum-cost assignment for a square matrix. Returns col index for every row."""
    C = np.asarray(cost, float)
    n = C.shape[0]
    INF = float("inf")
    u = np.zeros(n + 1)
    v = np.zeros(n + 1)
    p = np.zeros(n + 1, dtype=int)
    way = np.zeros(n + 1, dtype=int)
    for i in range(1, n + 1):
        p[0] = i
        j0 = 0
        minv = np.full(n + 1, INF)
        used = np.zeros(n + 1, dtype=bool)
        while True:
            used[j0] = True
            i0 = p[j0]
            delta, j1 = INF, 0
            for j in range(1, n + 1):
                if not used[j]:
                    cur = C[i0 - 1, j - 1] - u[i0] - v[j]
                    if cur < minv[j]:
                        minv[j] = cur
                        way[j] = j0
                    if minv[j] < delta:
                        delta, j1 = minv[j], j
            for j in range(n + 1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                else:
                    minv[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while True:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
            if j0 == 0:
                break
    out = [0] * n
    for j in range(1, n + 1):
        if p[j]:
            out[p[j] - 1] = j - 1
    return out


def kabsch(P: np.ndarray, Q: np.ndarray):
    """Rotation R and translation t so that (Q @ R.T + t) best fits P (proper rotation)."""
    pc, qc = P.mean(axis=0), Q.mean(axis=0)
    H = (Q - qc).T @ (P - pc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T)) or 1.0
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    return R, pc - qc @ R.T


def rmsd_after_fit(P, Q):
    R, t = kabsch(P, Q)
    return float(np.sqrt(((Q @ R.T + t - P) ** 2).sum(axis=1).mean()))


# ======================================================================================
# atom mapping
# ======================================================================================
@dataclass
class Mapping:
    perm: list            # perm[i] = product atom that becomes atom i (matches reactant atom i)
    core_size: int        # heavy atoms mapped by the common substructure
    rmsd: float           # all-atom RMSD after mapping + alignment (Å)
    broken: list = field(default_factory=list)   # [(i, j)] bonds only in the reactant
    formed: list = field(default_factory=list)   # [(i, j)] bonds only in the product
    notes: list = field(default_factory=list)
    alternatives: list = field(default_factory=list)  # other equally simple perms, best first


def _formula(m):
    out = {}
    for a in m.GetAtoms():
        out[a.GetSymbol()] = out.get(a.GetSymbol(), 0) + 1
    return out


def _bonds(m, order=None):
    inv = None if order is None else {old: new for new, old in enumerate(order)}
    s = set()
    for b in m.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if inv is not None:
            i, j = inv[i], inv[j]
        s.add((min(i, j), max(i, j)))
    return s


def _any_bonds(smarts: str) -> str:
    """Make every bond of an MCS SMARTS match any bond (the MCS keeps the bond types it saw)."""
    out, i, depth = [], 0, 0
    bond_chars = set("-=#:~@,;!/\\")
    while i < len(smarts):
        c = smarts[i]
        if c == "[":
            depth += 1
        elif c == "]":
            depth -= 1
        if depth == 0 and c in bond_chars:
            j = i
            while j < len(smarts) and smarts[j] in bond_chars:
                j += 1
            out.append("~")
            i = j
            continue
        out.append(c)
        i += 1
    return "".join(out)


def map_atoms(react: Chem.Mol, prod: Chem.Mol, timeout: int = 10, max_candidates: int = 300) -> Mapping:
    """
    Find which product atom corresponds to each reactant atom, following the principle of
    minimal structural change:
      1) maximum common substructure of the heavy atoms (bond orders ignored, so bonds that only
         change order still match); every symmetry-equivalent placement is a candidate,
      2) the remaining heavy atoms and all hydrogens are assigned optimally (Hungarian method):
         keep neighbours / the parent atom wherever possible, distance after superposition breaks ties,
      3) the candidate with the fewest broken + formed bonds wins; the best 3D fit decides ties.
    """
    n = react.GetNumAtoms()
    if prod.GetNumAtoms() != n or _formula(react) != _formula(prod):
        raise ValueError("Reactant and product must contain the same atoms "
                         f"(reactant {core.formula_and_mass(react)[0]}, "
                         f"product {core.formula_and_mass(prod)[0]}). Add or remove hydrogens / "
                         "fragments so both sides have the same formula.")
    P = react.GetConformer().GetPositions()
    Q = prod.GetConformer().GetPositions()
    notes = []
    rnb = [[x.GetIdx() for x in a.GetNeighbors()] for a in react.GetAtoms()]
    pnb = [set(x.GetIdx() for x in a.GetNeighbors()) for a in prod.GetAtoms()]
    rbonds = _bonds(react)
    pbonds_raw = _bonds(prod)

    def heavy_copy(m):
        rw = Chem.RWMol(m)
        for a in rw.GetAtoms():
            a.SetIntProp("orig", a.GetIdx())
        for i in sorted([a.GetIdx() for a in rw.GetAtoms() if a.GetAtomicNum() == 1], reverse=True):
            rw.RemoveAtom(i)
        h = rw.GetMol()
        h.UpdatePropertyCache(strict=False)
        return h

    # ---- candidate placements of the common core
    rh, ph = heavy_copy(react), heavy_copy(prod)
    candidates = []
    if rh.GetNumAtoms() >= 2:
        prm = rdFMCS.MCSParameters()
        prm.AtomTyper = rdFMCS.AtomCompare.CompareElements
        prm.BondTyper = rdFMCS.BondCompare.CompareAny
        prm.BondCompareParameters.RingMatchesRingOnly = False
        prm.BondCompareParameters.CompleteRingsOnly = False
        prm.Timeout = timeout
        res = rdFMCS.FindMCS([rh, ph], prm)
        if res.canceled:
            notes.append("Common-substructure search hit its time limit; result may be partial.")
        q = Chem.MolFromSmarts(_any_bonds(res.smartsString)) if res.smartsString else None
        if q is not None and q.GetNumAtoms() > 0:
            mr = rh.GetSubstructMatch(q)
            ri = [rh.GetAtomWithIdx(k).GetIntProp("orig") for k in mr]
            for mp in ph.GetSubstructMatches(q, uniquify=False, maxMatches=max_candidates, useChirality=False):
                pi = [ph.GetAtomWithIdx(k).GetIntProp("orig") for k in mp]
                candidates.append(dict(zip(ri, pi)))
    if not candidates:
        candidates = [{}]

    def superposed(mapping):
        keys = sorted(mapping)
        if len(keys) >= 3:
            R, t = kabsch(P[keys], Q[[mapping[k] for k in keys]])
            return Q @ R.T + t
        return Q - Q.mean(axis=0) + P.mean(axis=0)

    def complete(partial):
        mapping = dict(partial)
        Qs = superposed(mapping)
        for sym in sorted({a.GetSymbol() for a in react.GetAtoms() if a.GetAtomicNum() > 1}):
            rs = [a.GetIdx() for a in react.GetAtoms() if a.GetSymbol() == sym and a.GetIdx() not in mapping]
            if not rs:
                continue
            used = set(mapping.values())
            ps = [a.GetIdx() for a in prod.GetAtoms() if a.GetSymbol() == sym and a.GetIdx() not in used]
            C = np.zeros((len(rs), len(ps)))
            for a, i in enumerate(rs):
                want = {mapping[k] for k in rnb[i] if k in mapping}
                for b, j in enumerate(ps):
                    C[a, b] = 10.0 * len(want ^ (pnb[j] & set(mapping.values()))) + np.linalg.norm(P[i] - Qs[j])
            for a, b in enumerate(hungarian(C)):
                mapping[rs[a]] = ps[b]
        # hydrogens: stay on the same (mapped) parent whenever the counts allow it
        rH = [a.GetIdx() for a in react.GetAtoms() if a.GetAtomicNum() == 1]
        pH = [a.GetIdx() for a in prod.GetAtoms() if a.GetAtomicNum() == 1]
        if rH:
            Qs = superposed(mapping)
            C = np.zeros((len(rH), len(pH)))
            for a, i in enumerate(rH):
                rpar = {mapping.get(k, -1) for k in rnb[i]}
                for b, j in enumerate(pH):
                    C[a, b] = (0.0 if rpar & pnb[j] else 100.0) + np.linalg.norm(P[i] - Qs[j])
            for a, b in enumerate(hungarian(C)):
                mapping[rH[a]] = pH[b]
        return mapping

    # ---- heavy-atom mapping that minimises the number of bond changes (local search)
    rheavy = [a.GetIdx() for a in react.GetAtoms() if a.GetAtomicNum() > 1]
    elem = {i: react.GetAtomWithIdx(i).GetSymbol() for i in range(n)}
    pelem = {j: prod.GetAtomWithIdx(j).GetSymbol() for j in range(n)}
    ER = {(i, j) for i, j in rbonds if i in set(rheavy) and j in set(rheavy)}
    pheavy = {a.GetIdx() for a in prod.GetAtoms() if a.GetAtomicNum() > 1}
    EP = {(i, j) for i, j in pbonds_raw if i in pheavy and j in pheavy}
    hr = {i: sum(1 for k in rnb[i] if elem[k] == "H") for i in rheavy}
    hp = {j: sum(1 for k in pnb[j] if pelem[k] == "H") for j in pheavy}

    def cost(hm):
        common = sum(1 for i, j in ER if (min(hm[i], hm[j]), max(hm[i], hm[j])) in EP)
        return len(ER) + len(EP) - 2 * common + sum(abs(hr[i] - hp[hm[i]]) for i in rheavy)

    def improve(hm):
        hm = dict(hm)
        c = cost(hm)
        groups = {}
        for i in rheavy:
            groups.setdefault(elem[i], []).append(i)
        better = True
        while better:
            better = False
            for g in groups.values():
                for a in range(len(g)):
                    for b in range(a + 1, len(g)):
                        i, j = g[a], g[b]
                        hm[i], hm[j] = hm[j], hm[i]
                        c2 = cost(hm)
                        if c2 < c:
                            c, better = c2, True
                        else:
                            hm[i], hm[j] = hm[j], hm[i]
        return hm, c

    rng = np.random.default_rng(7)
    starts = []
    for cand in candidates[:60]:
        full = complete(cand)
        starts.append({i: full[i] for i in rheavy})
    for _ in range(25):                                   # random restarts
        hm = {}
        for sym in {elem[i] for i in rheavy}:
            rs = [i for i in rheavy if elem[i] == sym]
            ps = [j for j in sorted(pheavy) if pelem[j] == sym]
            rng.shuffle(ps)
            hm.update(dict(zip(rs, ps)))
        starts.append(hm)
    pool = {}
    for hm in starts:
        hm2, c = improve(hm) if rheavy else (hm, 0)
        pool.setdefault(c, {})[tuple(sorted(hm2.items()))] = hm2
    cmin = min(pool)
    scored = []
    for hm in list(pool[cmin].values())[:60]:             # ties: least motion, then best 3D fit
        mapping = complete(hm)
        perm = [mapping[i] for i in range(n)]
        inv = {p: i for i, p in enumerate(perm)}
        pb = {(min(inv[i], inv[j]), max(inv[i], inv[j])) for i, j in pbonds_raw}
        changes = len(rbonds ^ pb)
        Qm = Q[perm]
        R, t = kabsch(P, Qm)
        rmsd = float(np.sqrt(((Qm @ R.T + t - P) ** 2).sum(axis=1).mean()))
        # formed bonds should already be short in the reactant, broken ones short in the product
        motion = sum(np.linalg.norm(P[i] - P[j]) for i, j in pb - rbonds) + \
            sum(np.linalg.norm(Qm[i] - Qm[j]) for i, j in rbonds - pb)
        scored.append(((changes, round(motion, 2), round(rmsd, 3)), perm, pb, rmsd))
    scored.sort(key=lambda x: x[0])
    _, perm, pb, rmsd = scored[0]
    alternatives = []
    seen_changes = {(tuple(sorted(rbonds - pb)), tuple(sorted(pb - rbonds)))}
    for sc, p2, pb2, _ in scored[1:]:
        if sc[0] != scored[0][0][0]:
            continue
        key = (tuple(sorted(rbonds - pb2)), tuple(sorted(pb2 - rbonds)))
        if key not in seen_changes:                      # only chemically different numberings
            seen_changes.add(key)
            alternatives.append(p2)
    core_size = max((len(c) for c in candidates), default=0)
    if len(candidates) >= max_candidates:
        notes.append(f"Checked the first {max_candidates} placements of the common core.")
    return Mapping(perm, core_size, rmsd, sorted(rbonds - pb), sorted(pb - rbonds), notes, alternatives)


def renumber(mol: Chem.Mol, perm) -> Chem.Mol:
    """New molecule whose atom i is the old atom perm[i] (coordinates kept)."""
    return Chem.RenumberAtoms(mol, [int(k) for k in perm])


def swap_atoms(mol: Chem.Mol, i: int, j: int) -> Chem.Mol:
    perm = list(range(mol.GetNumAtoms()))
    perm[i], perm[j] = perm[j], perm[i]
    return renumber(mol, perm)


def move_label(mol: Chem.Mol, atom: int, new_index: int) -> Chem.Mol:
    """Give `atom` the number `new_index` (the atom that had it takes `atom`'s old number)."""
    return swap_atoms(mol, atom, new_index)


def align_onto(mobile: Chem.Mol, ref: Chem.Mol, conf_m=0, conf_r=0) -> float:
    """Rigidly superpose mobile onto ref (same atom order). Returns RMSD (Å)."""
    P = ref.GetConformer(conf_r).GetPositions()
    Q = mobile.GetConformer(conf_m).GetPositions()
    R, t = kabsch(P, Q)
    X = Q @ R.T + t
    core.apply_positions(mobile, conf_m, X)
    return float(np.sqrt(((X - P) ** 2).sum(axis=1).mean()))


def describe_changes(m: Mapping, react: Chem.Mol) -> str:
    lab = lambda i: f"{react.GetAtomWithIdx(i).GetSymbol()}{i}"  # noqa: E731
    parts = []
    if m.broken:
        parts.append("breaks " + ", ".join(f"{lab(i)}–{lab(j)}" for i, j in m.broken))
    if m.formed:
        parts.append("forms " + ", ".join(f"{lab(i)}–{lab(j)}" for i, j in m.formed))
    return "; ".join(parts) or "no bonds change"


# ======================================================================================
# pysisyphus input
# ======================================================================================
@dataclass
class TSSettings:
    coord_type: str = "dlc"           # dlc | redund | cart
    cos_type: str = "gs"              # gs (growing string) | neb | szts
    max_nodes: int = 18
    climb: bool = True
    opt_type: str = "string"
    align: bool = False
    max_cycles: int = 40
    tsopt_type: str = "rsirfo"        # rsirfo | rsprfo | trim
    do_hess: bool = True
    thresh: str = "gau"
    hessian_recalc: int = 1
    irc_type: str = "eulerpc"
    rms_grad_thresh: float = 0.0005
    preopt: bool = True
    endopt: bool = True
    calc_type: str = "xtb"
    gfn: int = 2
    pal: int = 6
    etemp: float = 1500.0
    charge: int = 0
    mult: int = 1
    solvent: str = "none"


def yaml_input(s: TSSettings, fn_r="reactant.xyz", fn_p="product.xyz") -> str:
    L = ["geom:", f" type: {s.coord_type}", f" fn: [{fn_r}, {fn_p}]"]
    L.append("preopt:" if s.preopt else "# preopt: (off)")
    L += ["cos:", f" type: {s.cos_type}", f" max_nodes: {s.max_nodes}", f" climb: {s.climb}"]
    L += ["opt:", f" type: {s.opt_type}", f" align: {s.align}", f" max_cycles: {s.max_cycles}"]
    L += ["tsopt:", f" type: {s.tsopt_type}", f" do_hess: {s.do_hess}", f" thresh: {s.thresh}",
          f" hessian_recalc: {s.hessian_recalc}"]
    L += ["irc:", f" type: {s.irc_type}", f" rms_grad_thresh: {s.rms_grad_thresh}"]
    L.append("endopt:" if s.endopt else "# endopt: (off)")
    L += ["calc:", f" type: {s.calc_type}", f" pal: {s.pal}"]
    if s.calc_type == "xtb":
        if s.gfn != 2:
            L.append(f" gfn: {s.gfn}")
        L.append(f" etemp: {s.etemp:g}")
        if s.solvent and s.solvent != "none":
            L.append(f" alpb: {s.solvent}")
    L += [f" charge: {s.charge}", f" mult: {s.mult}"]
    return "\n".join(L) + "\n"


# ======================================================================================
# running pysis
# ======================================================================================
PYSIS_VERSION = "1.0.0"

# Portable CPython used when Mol3D Bridge installs pysisyphus by itself (python-build-standalone,
# release 20260929; checksums pinned here so a tampered download is refused).
_PBS = "https://github.com/astral-sh/python-build-standalone/releases/download/20260929/"
PORTABLE_PYTHON = {
    "win": ("cpython-3.12.14+20260929-x86_64-pc-windows-msvc-install_only.tar.gz",
            "28728baf30b65e263f0b25c5a85be8226e7ab0d212fbadd6a8f0f796139fa804"),
    "linux": ("cpython-3.12.14+20260929-x86_64-unknown-linux-gnu-install_only.tar.gz",
              "06c90b93f419b63371c18f20fed0558a1a901f6518c3c24f755077e048447e7f"),
    "mac-arm": ("cpython-3.12.14+20260929-aarch64-apple-darwin-install_only.tar.gz",
                "de6b8f94fa765639b423ea353ab340669704c7186f96ee3cab389dcfde770c3c"),
    "mac-x86": ("cpython-3.12.14+20260929-x86_64-apple-darwin-install_only.tar.gz",
                "f51ec8a7fa0ede129a5e2e942a2a453ba4830adaca2289da4449390e42443292"),
}

WIN = sys.platform.startswith("win")
_NOWIN = {"creationflags": 0x08000000} if WIN else {}


def _scripts_dirs():
    """Folders where pip puts commands, for the usual Python installations."""
    import glob
    home = os.path.expanduser("~")
    pats = []
    if WIN:
        la, ad = os.environ.get("LOCALAPPDATA", ""), os.environ.get("APPDATA", "")
        pats += [os.path.join(la, "Python", "pythoncore-*", "Scripts"),          # Python install manager
                 os.path.join(la, "Programs", "Python", "Python3*", "Scripts"),   # python.org installer
                 os.path.join(ad, "Python", "Python3*", "Scripts"),               # pip install --user
                 r"C:\Program Files\Python3*\Scripts",
                 os.path.join(home, "miniconda3", "Scripts"), os.path.join(home, "anaconda3", "Scripts"),
                 os.path.join(home, "miniforge3", "Scripts"),
                 os.path.join(la, "Packages", "PythonSoftwareFoundation.Python.3*", "LocalCache",
                              "local-packages", "Python3*", "Scripts")]
    else:
        pats += [os.path.join(home, ".local", "bin"), "/opt/homebrew/bin", "/usr/local/bin",
                 os.path.join(home, "miniconda3", "bin"), os.path.join(home, "anaconda3", "bin"),
                 os.path.join(home, "miniforge3", "bin"), os.path.join(home, "Library", "Python", "3.*", "bin")]
    out = []
    for p in pats:
        out += sorted(glob.glob(p), reverse=True)
    return out


def managed_pysis(base_dir: str) -> str:
    """Where Mol3D Bridge's own pysisyphus installation keeps the pysis command."""
    env = os.path.join(base_dir, "env")
    return os.path.join(env, "Scripts", "pysis.exe") if WIN else os.path.join(env, "bin", "pysis")


def find_pysis(hint: str = "", managed_dir: str = "") -> str | None:
    if hint and os.path.isfile(hint):
        return hint
    if managed_dir and os.path.isfile(managed_pysis(managed_dir)):
        return managed_pysis(managed_dir)
    w = shutil.which("pysis")
    if w:
        return w
    dirs = [os.path.dirname(sys.executable), os.path.join(os.path.dirname(sys.executable), "Scripts")]
    for d in dirs + _scripts_dirs():
        for name in ("pysis.exe", "pysis") if WIN else ("pysis",):
            p = os.path.join(d, name)
            if os.path.isfile(p):
                return p
    return None


def _clean_env():
    env = dict(os.environ)
    for k in ("PYTHONHOME", "PYTHONPATH", "PYTHONSTARTUP", "VIRTUAL_ENV", "__PYVENV_LAUNCHER__"):
        env.pop(k, None)
    env.update(PYTHONIOENCODING="utf-8", PYTHONUTF8="1", PIP_PROGRESS_BAR="off",
               PIP_NO_INPUT="1", PIP_DISABLE_PIP_VERSION_CHECK="1", PYTHONUNBUFFERED="1")
    return env


def _stream(cmd, progress, cancel, cwd=None):
    proc = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            env=_clean_env(), text=True, encoding="utf-8", errors="replace",
                            bufsize=1, **_NOWIN)
    for line in proc.stdout:
        if progress:
            progress(line.rstrip())
        if cancel is not None and cancel.is_set():
            proc.kill()
            proc.wait()
            raise core.Cancelled()
    rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"'{os.path.basename(cmd[0])} {' '.join(cmd[1:3])}…' failed (exit code {rc}).")


def _portable_key():
    import platform
    if WIN:
        return "win"
    if sys.platform == "darwin":
        return "mac-arm" if platform.machine() in ("arm64", "aarch64") else "mac-x86"
    return "linux"


def download_portable_python(dest: str, progress=None, cancel=None) -> str:
    """Download + verify + unpack a private Python 3.12. Returns the python executable."""
    import hashlib
    import tarfile
    import urllib.request
    name, sha = PORTABLE_PYTHON[_portable_key()]
    os.makedirs(dest, exist_ok=True)
    arc = os.path.join(dest, name)
    if progress:
        progress(f"Downloading Python 3.12 ({name.split('-install')[0]}) …")
    with urllib.request.urlopen(_PBS + name.replace("+", "%2B"), timeout=60) as r, open(arc, "wb") as fh:
        total = int(r.headers.get("Content-Length") or 0)
        got, last = 0, -10
        while True:
            if cancel is not None and cancel.is_set():
                raise core.Cancelled()
            chunk = r.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
            got += len(chunk)
            pct = int(got * 100 / total) if total else 0
            if progress and pct >= last + 10:
                last = pct
                progress(f"  {pct} %")
    h = hashlib.sha256()
    with open(arc, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    if h.hexdigest() != sha:
        os.remove(arc)
        raise RuntimeError("The Python download is corrupted (checksum mismatch). Please try again.")
    if progress:
        progress("Checksum OK – unpacking …")
    with tarfile.open(arc) as t:
        try:
            t.extractall(dest, filter="data")
        except TypeError:
            t.extractall(dest)
    os.remove(arc)
    exe = os.path.join(dest, "python", "python.exe") if WIN else os.path.join(dest, "python", "bin", "python3")
    if not os.path.isfile(exe):
        raise RuntimeError("Unpacked Python, but python executable is missing.")
    return exe


def install_pysis(base_dir: str, progress=None, cancel=None) -> str:
    """
    One-click installation: private Python 3.12 -> virtual environment -> pip install pysisyphus.
    Everything goes into base_dir (Mol3D Bridge's data folder); nothing else on the PC is changed.
    """
    os.makedirs(base_dir, exist_ok=True)
    env = os.path.join(base_dir, "env")
    py_env = os.path.join(env, "Scripts", "python.exe") if WIN else os.path.join(env, "bin", "python")
    if not os.path.isfile(py_env):
        base = os.path.join(base_dir, "python", "python.exe") if WIN else \
            os.path.join(base_dir, "python", "bin", "python3")
        if not os.path.isfile(base):
            base = download_portable_python(base_dir, progress, cancel)
        if progress:
            progress("Creating a private environment …")
        _stream([base, "-m", "venv", env], progress, cancel)
    if progress:
        progress("Installing pysisyphus and its packages (this takes a few minutes) …")
    _stream([py_env, "-m", "pip", "install", "--upgrade", "pip"], progress, cancel)
    _stream([py_env, "-m", "pip", "install", f"pysisyphus=={PYSIS_VERSION}"], progress, cancel)
    exe = managed_pysis(base_dir)
    if not os.path.isfile(exe):
        raise RuntimeError("pysisyphus was installed, but the pysis command was not found.")
    _stream([py_env, "-c", "import pysisyphus, sys; print('pysisyphus', pysisyphus.__version__, "
             "'ready (Python', sys.version.split()[0] + ')')"], progress, cancel)
    return exe


def write_inputs(workdir: str, react: Chem.Mol, prod: Chem.Mol, settings: TSSettings,
                 yaml_text: str | None, xtb_path: str = "") -> str:
    os.makedirs(workdir, exist_ok=True)
    with open(os.path.join(workdir, "reactant.xyz"), "w", newline="\n") as fh:
        fh.write(core.xyz_block(react, react.GetConformer().GetId(), "reactant"))
    with open(os.path.join(workdir, "product.xyz"), "w", newline="\n") as fh:
        fh.write(core.xyz_block(prod, prod.GetConformer().GetId(), "product"))
    yml = os.path.join(workdir, "ts.yaml")
    with open(yml, "w", newline="\n") as fh:
        fh.write(yaml_text if yaml_text is not None else yaml_input(settings))
    rc = os.path.join(workdir, "pysisyphusrc")
    with open(rc, "w", newline="\n") as fh:
        if xtb_path:
            fh.write(f"[xtb]\ncmd={xtb_path}\n")
    return yml


LAUNCHER = """# Written by Mol3D Bridge: starts pysisyphus like the `pysis` command does.
import os
import sys

_symlink = os.symlink


def _symlink_or_skip(src, dst, *args, **kwargs):
    # pysisyphus links qm_calcs/cur_out to the running calculation's output - a convenience only.
    # Windows allows symbolic links only with administrator rights or Developer Mode
    # ([WinError 1314]); without them the link is skipped instead of aborting the run.
    try:
        return _symlink(src, dst, *args, **kwargs)
    except OSError as err:
        if os.name == "nt" or getattr(err, "winerror", None) == 1314:
            return None
        raise


os.symlink = _symlink_or_skip
from pysisyphus.run import run  # noqa: E402

sys.argv[0] = "pysis"
sys.exit(run())
"""


def pysis_python(pysis: str) -> str | None:
    """The Python interpreter that belongs to a `pysis` command (None if it cannot be found)."""
    base = os.path.basename(pysis).lower()
    if base.startswith("python"):
        return pysis
    d = os.path.dirname(os.path.abspath(pysis))
    if WIN:
        for c in (os.path.join(d, "python.exe"), os.path.join(os.path.dirname(d), "python.exe")):
            if os.path.isfile(c):
                return c
        return None
    try:
        with open(pysis, "rb") as fh:
            first = fh.readline().decode("utf-8", "replace").strip()
        if first.startswith("#!"):
            exe = first[2:].strip().split()[0]
            if os.path.isfile(exe):
                return exe
    except OSError:
        pass
    for c in (os.path.join(d, "python3"), os.path.join(d, "python")):
        if os.path.isfile(c):
            return c
    return None


def run_pysis(pysis: str, workdir: str, xtb_path: str = "", progress=None, cancel=None) -> int:
    """Run pysisyphus on ts.yaml in workdir, streaming each output line through progress()."""
    env = dict(os.environ)
    env["PYSISRC"] = os.path.join(workdir, "pysisyphusrc")
    env["PYTHONUNBUFFERED"] = "1"
    # pysisyphus prints Unicode (e.g. C₅₄H₃₆, kJ mol⁻¹). On Windows a pipe would otherwise use the
    # ANSI code page (cp1252) and pysis crashes with UnicodeEncodeError - force UTF-8 everywhere.
    env["PYTHONIOENCODING"] = "utf-8"
    env["PYTHONUTF8"] = "1"
    env.setdefault("OMP_STACKSIZE", "1G")
    if xtb_path:
        env["PATH"] = os.path.dirname(xtb_path) + os.pathsep + env.get("PATH", "")
    py = pysis_python(pysis)
    if py:
        with open(os.path.join(workdir, "run_pysis.py"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(LAUNCHER)
        cmd = [py, "run_pysis.py", "ts.yaml"]
    else:
        cmd = [pysis, "ts.yaml"]
    kw = {}
    if sys.platform.startswith("win"):
        kw["creationflags"] = 0x08000000
    proc = subprocess.Popen(cmd, cwd=workdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            env=env, bufsize=1, text=True, encoding="utf-8", errors="replace", **kw)
    log = open(os.path.join(workdir, "pysis.log"), "w", encoding="utf-8")
    try:
        for line in proc.stdout:
            log.write(line)
            if progress:
                progress(line.rstrip("\n"))
            if cancel is not None and cancel.is_set():
                proc.kill()
                break
    finally:
        log.close()
    rc = proc.wait()
    if cancel is not None and cancel.is_set():
        raise core.Cancelled()
    return rc


@dataclass
class TSResult:
    ts: Chem.Mol | None = None
    minimum1: Chem.Mol | None = None      # "Left"  = optimised forward IRC end
    minimum2: Chem.Mol | None = None      # "Right" = optimised backward IRC end
    files: dict = field(default_factory=dict)
    energies: dict = field(default_factory=dict)    # Hartree
    relative_kj: dict = field(default_factory=dict) # from pysis' BARRIERS table
    sides: dict = field(default_factory=dict)       # minimum1/2 -> "reactant" | "product" | "other"
    imag_freq: float | None = None
    irc_frames: list = field(default_factory=list)  # [(positions Nx3, energy Eh)] along the IRC
    notes: list = field(default_factory=list)


def read_xyz_frames(path):
    """[(xyz_text, comment)] of a .xyz / .trj file."""
    with open(path, errors="replace") as fh:
        lines = fh.read().splitlines()
    frames, i = [], 0
    while i < len(lines):
        try:
            n = int(lines[i].split()[0])
        except (ValueError, IndexError):
            i += 1
            continue
        block = lines[i:i + n + 2]
        if len(block) < n + 2:
            break
        frames.append(("\n".join(block) + "\n", block[1]))
        i += n + 2
    return frames


def positions_of(xyz_text):
    rows = xyz_text.splitlines()[2:]
    return np.array([[float(v) for v in r.split()[1:4]] for r in rows if r.strip()])


def energy_of(comment):
    m = re.search(r"(-?\d+\.\d{5,})", comment or "")
    return float(m.group(1)) if m else None


def mol_with_positions(template: Chem.Mol, pos, perceive="connectivity") -> Chem.Mol:
    """Same atoms as template, new coordinates; bonds re-perceived from the geometry."""
    from rdkit.Chem import rdDetermineBonds
    rw = Chem.RWMol()
    for a in template.GetAtoms():
        na = Chem.Atom(a.GetAtomicNum())
        na.SetFormalCharge(a.GetFormalCharge())
        na.SetNoImplicit(True)
        rw.AddAtom(na)
    conf = Chem.Conformer(template.GetNumAtoms())
    for i, p in enumerate(pos):
        conf.SetAtomPosition(i, [float(x) for x in p])
    conf.Set3D(True)
    rw.AddConformer(conf, assignId=True)
    m = rw.GetMol()
    try:
        if perceive == "bonds":
            rdDetermineBonds.DetermineBonds(m, charge=Chem.GetFormalCharge(template))
        else:
            rdDetermineBonds.DetermineConnectivity(m)
    except Exception:  # noqa: BLE001
        try:
            rdDetermineBonds.DetermineConnectivity(m)
        except Exception:  # noqa: BLE001
            pass
    m.UpdatePropertyCache(strict=False)
    return m


def _side(m: Chem.Mol, react: Chem.Mol, prod: Chem.Mol) -> str:
    b = _bonds(m)
    dr, dp = len(b ^ _bonds(react)), len(b ^ _bonds(prod))
    if dr == 0 and dp > 0:
        return "reactant"
    if dp == 0 and dr > 0:
        return "product"
    if dr == dp:
        return "other"
    return "reactant" if dr < dp else "product"


def collect_results(workdir: str, react: Chem.Mol, prod: Chem.Mol | None = None) -> TSResult:
    r = TSResult()
    trio = os.path.join(workdir, "left_ts_right_geoms.trj")
    names = {"minimum1": ["forward_end_opt.xyz", "forward_end.xyz", "forward_last.xyz"],
             "ts": ["ts_opt.xyz", "ts_final_geometry.xyz", "splined_hei.xyz"],
             "minimum2": ["backward_end_opt.xyz", "backward_end.xyz", "backward_last.xyz"]}
    found = {}
    if os.path.isfile(trio):
        fr = read_xyz_frames(trio)
        if len(fr) == 3:
            for key, (txt, com) in zip(("minimum1", "ts", "minimum2"), fr):
                found[key] = (txt, com, "left_ts_right_geoms.trj")
    for key, cands in names.items():
        for nm in cands:
            p = os.path.join(workdir, nm)
            if key not in found and os.path.isfile(p):
                fr = read_xyz_frames(p)
                if fr:
                    found[key] = (fr[-1][0], fr[-1][1], nm)
            if os.path.isfile(p) and key in found and energy_of(found[key][1]) is None:
                fr = read_xyz_frames(p)
                if fr and energy_of(fr[-1][1]) is not None:
                    found[key] = (found[key][0], fr[-1][1], found[key][2])
    for key, (txt, com, nm) in found.items():
        try:
            pos = positions_of(txt)
            m = mol_with_positions(react, pos, "connectivity" if key == "ts" else "bonds")
            setattr(r, key, m)
            r.files[key] = nm
            e = energy_of(com)
            if e is not None:
                r.energies[key] = e
        except Exception as ex:  # noqa: BLE001
            r.notes.append(f"{nm}: {ex}")
    irc = os.path.join(workdir, "finished_irc.trj")
    if os.path.isfile(irc):
        for txt, com in read_xyz_frames(irc):
            try:
                r.irc_frames.append((positions_of(txt), energy_of(com)))
            except ValueError:
                pass
    log = os.path.join(workdir, "pysis.log")
    if os.path.isfile(log):
        with open(log, encoding="utf-8", errors="replace") as fh:
            text = fh.read()
        im = re.findall(r"Imaginary frequencies:\s*\[([^\]]*)\]", text)
        if im:
            vals = [float(v) for v in re.findall(r"-?\d+\.?\d*", im[-1])]
            r.imag_freq = vals[0] if vals else None
            if len(vals) > 1:
                r.notes.append(f"{len(vals)} imaginary frequencies at the final TS: {vals}")
        if "BARRIERS" in text:
            tail = text[text.rindex("BARRIERS"):]
            for key, lab in (("minimum1", "Left"), ("ts", "TS"), ("minimum2", "Right")):
                m = re.search(rf"^\s*{lab}:\s+(-?\d+\.\d+)\s*kJ", tail, re.M)
                if m:
                    r.relative_kj[key] = float(m.group(1))
    if prod is not None:
        for key in ("minimum1", "minimum2"):
            m = getattr(r, key)
            if m is not None:
                r.sides[key] = _side(m, react, prod)
    if r.ts is None:
        r.notes.append("No TS geometry was written – see the output log.")
    return r


STAGES = [("preopt", "Pre-optimizing the end points"), ("cos", "Growing string / chain of states"),
          ("tsopt", "TS optimization (+ Hessian)"), ("irc", "IRC"), ("endopt", "Optimizing the IRC ends"),
          ("done", "Finished")]


def current_stage(workdir: str) -> tuple:
    """(stage key, detail) from the files pysis has written so far."""
    import glob
    has = lambda pat: bool(glob.glob(os.path.join(workdir, pat)))  # noqa: E731
    if has("left_ts_right_geoms.trj"):
        return "done", ""
    if has("*_end_*"):
        return "endopt", ""
    if has("*_irc.trj") or has("irc_data.h5"):
        return "irc", ""
    if has("ts_*") and not has("ts_guess*") or has("ts_opt*"):
        return "tsopt", ""
    cyc = sorted(glob.glob(os.path.join(workdir, "cycle_*.trj")))
    if cyc:
        return "cos", os.path.basename(cyc[-1])
    if has("*pre*"):
        return "preopt", ""
    return "start", ""


def latest_path(workdir: str):
    """Images (positions) and energies (Eh) of the newest chain-of-states cycle, or None."""
    import glob
    cyc = sorted(glob.glob(os.path.join(workdir, "cycle_*.trj")))
    if not cyc:
        return None
    try:
        fr = read_xyz_frames(cyc[-1])
    except OSError:
        return None
    if not fr:
        return None
    return (os.path.basename(cyc[-1]), [positions_of(t) for t, _ in fr], [energy_of(c) for _, c in fr])


# ======================================================================================
# permanent results folder
# ======================================================================================
def _safe_name(text: str) -> str:
    text = re.sub(r"[^\w\-. ]+", "_", text.strip(), flags=re.UNICODE)
    return re.sub(r"_+", "_", re.sub(r"\s+", "_", text)).strip("._")[:60]


def unique_folder(root: str, name: str) -> str:
    os.makedirs(root, exist_ok=True)
    path = os.path.join(root, name)
    k = 2
    while os.path.exists(path):
        path = os.path.join(root, f"{name}_{k}")
        k += 1
    os.makedirs(path)
    return path


def save_results(workdir: str, results_root: str, react: Chem.Mol, prod: Chem.Mol,
                 run_name: str = "", status: str = "finished", keep_all_files: bool = True,
                 started: float | None = None) -> str:
    """
    Store a finished (or stopped / failed) run in its own, clearly named folder:

      <root>/2026-10-20_14-15_<name>/
        summary.txt               energies, barriers, imaginary frequency, bond changes, settings
        TS.xyz  minimum_reactant-side.xyz  minimum_product-side.xyz
        results.sdf               the three structures with energies (opens in ChemDraw, PyMOL, …)
        IRC_path.xyz              whole reaction path (animates in Mol3D Bridge, Avogadro, VMD, …)
        reactant_input.xyz  product_input.xyz  ts.yaml  pysis.log
        pysisyphus/               every file pysisyphus wrote (optional)
    """
    import datetime
    stamp = datetime.datetime.fromtimestamp(started or time.time()).strftime("%Y-%m-%d_%H-%M")
    formula = core.formula_and_mass(react)[0]
    name = f"{stamp}_{_safe_name(run_name) or formula}"
    if status != "finished":
        name += f"_{status}"
    dest = unique_folder(results_root, name)

    res = collect_results(workdir, react, prod)
    side = {"reactant": "reactant-side", "product": "product-side", "other": "other"}
    names = {"ts": "TS"}
    for key, default in (("minimum1", "minimum_A"), ("minimum2", "minimum_B")):
        s_ = res.sides.get(key)
        names[key] = f"minimum_{side[s_]}" if s_ in ("reactant", "product") else default
    if names["minimum1"] == names["minimum2"]:
        names["minimum1"] += "_A"
        names["minimum2"] += "_B"

    def write(fn, text):
        with open(os.path.join(dest, fn), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)

    sdf = Chem.SDWriter(os.path.join(dest, "results.sdf"))
    for key in ("ts", "minimum1", "minimum2"):
        m = getattr(res, key)
        if m is None:
            continue
        e = res.energies.get(key)
        comment = f"{names[key]}  {formula}" + (f"  E={e:.8f} Eh" if e is not None else "")
        write(f"{names[key]}.xyz", core.xyz_block(m, 0, comment))
        mm = Chem.Mol(m)
        mm.SetProp("_Name", names[key])
        if e is not None:
            mm.SetProp("Energy_Eh", f"{e:.8f}")
        if key in res.relative_kj:
            mm.SetProp("Relative_kJ_mol", f"{res.relative_kj[key]:.2f}")
        if key == "ts" and res.imag_freq is not None:
            mm.SetProp("Imaginary_frequency_cm-1", f"{res.imag_freq:.2f}")
        sdf.write(mm)
    sdf.close()
    if res.irc_frames:
        syms = [a.GetSymbol() for a in react.GetAtoms()]
        blocks = []
        for k, (pos, e) in enumerate(res.irc_frames):
            lines = [f"{len(syms)}", f"IRC point {k}" + (f"  E={e:.8f} Eh" if e is not None else "")]
            lines += [f"{s:<2s} {x:14.8f} {y:14.8f} {z:14.8f}" for s, (x, y, z) in zip(syms, pos)]
            blocks.append("\n".join(lines))
        write("IRC_path.xyz", "\n".join(blocks) + "\n")
    for src, dst in (("reactant.xyz", "reactant_input.xyz"), ("product.xyz", "product_input.xyz"),
                     ("ts.yaml", "ts.yaml"), ("pysis.log", "pysis.log")):
        p = os.path.join(workdir, src)
        if os.path.isfile(p):
            shutil.copy2(p, os.path.join(dest, dst))

    # ---- human-readable summary
    rb, pb = _bonds(react), _bonds(prod)
    chg = Mapping(list(range(react.GetNumAtoms())), 0, 0.0, sorted(rb - pb), sorted(pb - rb))
    L = ["Mol3D Bridge – transition-state search (pysisyphus)", "=" * 52, "",
         f"Status:      {status}", f"Formula:     {formula}  ({react.GetNumAtoms()} atoms)",
         f"Reaction:    {describe_changes(chg, react)}", f"Saved:       {dest}", ""]
    if started:
        L.insert(5, f"Run time:    {time.time() - started:.0f} s")
    L.append("Structure                       E (Eh)            rel. (kJ/mol)")
    for key in ("minimum1", "ts", "minimum2"):
        if getattr(res, key) is None:
            continue
        e = res.energies.get(key)
        r_ = res.relative_kj.get(key)
        L.append(f"  {names[key]:<28s}{(f'{e:.8f}' if e is not None else '–'):>16s}"
                 f"{(f'{r_:+.2f}' if r_ is not None else '–'):>16s}")
    L.append("")
    if "ts" in res.relative_kj:
        for key in ("minimum1", "minimum2"):
            if key in res.relative_kj:
                L.append(f"Barrier from {names[key]}: {res.relative_kj['ts'] - res.relative_kj[key]:.2f} kJ/mol"
                         f"  ({(res.relative_kj['ts'] - res.relative_kj[key]) / 4.184:.2f} kcal/mol)")
    if res.imag_freq is not None:
        L.append(f"Imaginary frequency of the TS: {res.imag_freq:.2f} cm-1")
    if res.notes:
        L += [""] + [f"Note: {n}" for n in res.notes]
    yml = os.path.join(workdir, "ts.yaml")
    if os.path.isfile(yml):
        with open(yml, encoding="utf-8") as fh:
            L += ["", "pysisyphus input (ts.yaml):", "-" * 27, fh.read().rstrip()]
    L += ["", "Please cite:",
          "  pysisyphus – J. Steinmetzer, S. Kupfer, S. Gräfe, \"pysisyphus: Exploring potential energy surfaces",
          "    in ground and excited states\", Int. J. Quantum Chem. 2021, 121, e26390. doi:10.1002/qua.26390",
          "  GFN2-xTB – C. Bannwarth, S. Ehlert, S. Grimme, J. Chem. Theory Comput. 2019, 15, 1652–1671.",
          "    doi:10.1021/acs.jctc.8b01176",
          "  xtb program – C. Bannwarth, E. Caldeweyher, S. Ehlert, A. Hansen, P. Pracht, J. Seibert,",
          "    S. Spicher, S. Grimme, WIREs Comput. Mol. Sci. 2021, 11, e1493. doi:10.1002/wcms.1493"]
    write("summary.txt", "\n".join(L) + "\n")

    if keep_all_files:
        shutil.copytree(workdir, os.path.join(dest, "pysisyphus"),
                        ignore=shutil.ignore_patterns("qm_calcs", "*.tmp"), dirs_exist_ok=True)
    return dest


def find_run_folder(path: str) -> str | None:
    """Accept a saved results folder or a pysisyphus working folder."""
    for d in (path, os.path.join(path, "pysisyphus")):
        if os.path.isfile(os.path.join(d, "reactant.xyz")) and os.path.isfile(os.path.join(d, "product.xyz")):
            return d
    return None
