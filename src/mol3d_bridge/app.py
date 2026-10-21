#!/usr/bin/env python3
"""
Mol3D Bridge - a Chem3D-style companion for ChemDraw and ReShiki.

  Ctrl+V          paste a structure from ChemDraw / ReShiki -> 3D model (not minimized)
  Ctrl+M          minimize the current geometry again
  Ctrl+C          copy the current 3D view back as a 2D drawing (wedges from real depth)
  Ctrl+Shift+C    copy as a clean 2D drawing (stereo taken from the 3D model)
  Ctrl+D          copy for ChemDraw keeping 3D coordinates (rotatable in 3D, like Chem3D)
  Ctrl+Shift+D    copy for ReShiki keeping 3D coordinates (3D molfile -> ReShiki depth)
  Ctrl+E          ORCA input (preview, copy, save)

Rotate mode (R): left-drag rotate, Shift+drag spin, right/middle-drag (or Ctrl+drag) pan,
wheel zoom, click atoms to measure (2 = distance, 3 = angle, 4 = dihedral), Esc clears.

Edit-atoms mode (E) - adjust the starting geometry before minimizing:
  click atom / drag a box        select (Shift adds/toggles), double-click = whole molecule
  drag an atom                   move the selection in the screen plane (attached H follow)
  Alt+drag                       move toward / away from you (depth)
  Shift+drag a selected atom     rotate the selected group about its centre
  arrow keys / PgUp PgDn         nudge 0.05 Å (Shift: 0.5 Å)
  right-drag                     rotate the view while editing
  Ctrl+Shift+S                   select substituent (atom that stays, then atom whose group moves)
  Ctrl+G                         set distance / angle / dihedral of 2-4 selected atoms
  Ctrl+M                         minimize from the edited geometry (optionally holding selection fixed)

Methods: MMFF94 / MMFF94s / UFF (built in) or GFN2-xTB / GFN-FF (needs the free xtb program;
Calculate ▸ xtb settings… downloads it). Optimizations are animated step by step with a live
energy curve (Settings ▸ animation length, default 5 s). Esc / Stop skips or stops.
"""
from __future__ import annotations

import math
import os
import random
import sys
import threading
import time
import traceback

import numpy as np
from PySide6.QtCore import QPointF, QRectF, QSettings, QStandardPaths, Qt, QThread, QTimer, Signal
from PySide6.QtCore import QUrl
from PySide6.QtGui import (QDesktopServices, QAction, QActionGroup, QBrush, QColor, QFont, QKeySequence,
                           QPainter, QPen, QRadialGradient, QIcon, QPixmap)
from PySide6.QtWidgets import (QApplication, QComboBox, QMenu, QToolButton, QDialog, QDialogButtonBox,
                               QFileDialog, QFormLayout, QHBoxLayout, QInputDialog, QLabel,
                               QLineEdit, QMainWindow, QMessageBox, QPlainTextEdit,
                               QPushButton, QSpinBox, QToolBar, QVBoxLayout, QWidget,
                               QCheckBox, QDoubleSpinBox)
from rdkit import Chem

from . import APP_ID, APP_NAME, __version__
from . import clipboard as clip
from . import core

APP = APP_NAME
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources")


def app_base_dir():
    """Folder of the installed program (PyInstaller bundle) or of the package."""
    if getattr(sys, "frozen", False):
        return getattr(sys, "_MEIPASS", os.path.dirname(sys.executable))
    return os.path.dirname(os.path.abspath(__file__))


def user_data_dir():
    d = QStandardPaths.writableLocation(QStandardPaths.AppDataLocation) or os.path.expanduser("~/.mol3d_bridge")
    os.makedirs(d, exist_ok=True)
    return d

CPK = {
    "H": (238, 238, 238), "C": (125, 125, 125), "N": (48, 80, 248), "O": (240, 20, 20),
    "F": (144, 224, 80), "Cl": (31, 200, 31), "Br": (166, 41, 41), "I": (148, 0, 148),
    "S": (230, 200, 40), "P": (255, 128, 0), "B": (255, 170, 170), "Si": (220, 180, 140),
    "Li": (204, 128, 255), "Na": (171, 92, 242), "K": (143, 64, 212), "Mg": (138, 220, 0),
    "Ca": (61, 255, 0), "Fe": (224, 102, 51), "Co": (240, 144, 160), "Ni": (80, 208, 80),
    "Cu": (200, 128, 51), "Zn": (125, 128, 176), "Ru": (36, 143, 143), "Rh": (10, 125, 140),
    "Pd": (0, 105, 133), "Ag": (192, 192, 192), "Ir": (23, 84, 135), "Pt": (208, 208, 224),
    "Au": (255, 209, 35), "Sn": (102, 128, 128), "Se": (255, 161, 0), "Ti": (191, 194, 199),
    "Mn": (156, 122, 199), "Cr": (138, 153, 199), "Al": (191, 166, 166),
}
PT = Chem.GetPeriodicTable()


def color_of(sym):
    return QColor(*CPK.get(sym, (255, 20, 147)))


def rot(axis, ang):
    c, s = math.cos(ang), math.sin(ang)
    if axis == "x":
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    if axis == "y":
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


# ======================================================================================
# 3D viewer (pure QPainter: no OpenGL needed, works everywhere incl. remote desktop)
# ======================================================================================
class MolView(QWidget):
    selectionChanged = Signal(str)
    editStarted = Signal()          # before atoms are moved (for undo)
    editFinished = Signal(list)     # atoms that were moved
    modeChanged = Signal(str)

    STYLES = ("Ball & Stick", "Sticks", "Wireframe", "Space filling")

    def __init__(self):
        super().__init__()
        self.setMinimumSize(500, 400)
        self.setFocusPolicy(Qt.StrongFocus)
        self.mol = None
        self.disp = None
        self.conf_id = 0
        self.R = np.eye(3)
        self.zoom = 1.0
        self.pan = QPointF(0, 0)
        self.style = "Ball & Stick"
        self.show_h = True
        self.labels = "hetero"  # none | hetero | all | index
        self.dark = False
        self.sel = []
        self._press = None
        self._last = None
        self._moved = False
        self._screen = None
        self.mode = "rotate"          # rotate | edit
        self.move_h = True            # attached H follow moved heavy atoms
        self.center = np.zeros(3)
        self.extent = 5.0
        self.hover = None
        self._drag = None             # None | rotate | pan | move | depth | turn | box
        self._box = None
        self._moving = []
        self.locked = False           # True while an optimisation is animating (view-only)
        self.empty_text = None        # message shown while there is no molecule
        self.extra_bonds = []         # (i, j) drawn dashed, e.g. bonds that form/break at a TS
        self.highlight = None         # QColor: frame drawn around the view (e.g. the active panel)
        self.trace = None             # energy curve of the running / last optimisation
        self.setMouseTracking(True)

    # ---------------------------------------------------------------- model
    def set_molecule(self, mol, conf_id=0, keep_view=False):
        self.mol = mol
        self.conf_id = conf_id
        n = mol.GetNumAtoms() if mol is not None else 0
        self.sel = [i for i in self.sel if i < n] if keep_view else []
        self.hover = None
        self.disp = None
        if mol is not None:
            d = Chem.Mol(mol)
            try:
                Chem.Kekulize(d, clearAromaticFlags=True)
            except Exception:
                pass
            self.disp = d
        if not keep_view:
            self.R = core.best_view_rotation(mol, conf_id) if mol is not None else np.eye(3)
            self.zoom, self.pan = 1.0, QPointF(0, 0)
            self.recenter()
        self.selectionChanged.emit(self.measurement())
        self.update()

    def recenter(self):
        """Fix the rotation centre / scale (kept constant while atoms are being moved)."""
        if self.mol is None:
            return
        pos = self.positions()
        heavy = [a.GetIdx() for a in self.mol.GetAtoms() if a.GetAtomicNum() > 1] or list(range(len(pos)))
        self.center = pos[heavy].mean(axis=0)
        self.extent = max(1.5, float(np.max(np.linalg.norm(pos - self.center, axis=1))) + 1.0)

    def set_mode(self, mode):
        self.mode = mode
        self._drag = None
        self.setCursor(Qt.ArrowCursor)
        self.modeChanged.emit(mode)
        self.update()

    def positions(self):
        return self.mol.GetConformer(self.conf_id).GetPositions()

    def visible(self, i):
        return self.show_h or self.mol.GetAtomWithIdx(i).GetAtomicNum() != 1 or i in self.sel

    def atom_radius(self, sym):
        if self.style == "Space filling":
            return PT.GetRvdw(sym)
        if self.style == "Ball & Stick":
            return 0.22 * PT.GetRvdw(sym)
        if self.style == "Sticks":
            return 0.14
        return 0.0

    # ---------------------------------------------------------------- paint
    def paintEvent(self, _):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        bg = QColor(28, 30, 36) if self.dark else QColor(246, 247, 249)
        fg = QColor(230, 230, 230) if self.dark else QColor(40, 40, 46)
        p.fillRect(self.rect(), bg)
        if self.mol is None:
            p.setPen(fg)
            f = QFont(); f.setPointSize(12); p.setFont(f)
            p.drawText(self.rect(), Qt.AlignCenter, self.empty_text or
                       "Copy a structure in ChemDraw or ReShiki,\nthen press Ctrl+V here.\n\n"
                       "You can also drop a .cdx / .cdxml / .mol / .sdf / .xyz file,\n"
                       "or type a SMILES in the toolbar.")
            self._paint_highlight(p)
            return

        pos = self.positions()
        v = (pos - self.center) @ self.R.T
        s = 0.52 * min(self.width(), self.height()) / self.extent * self.zoom
        cx, cy = self.width() / 2 + self.pan.x(), self.height() / 2 + self.pan.y()
        sx = cx + v[:, 0] * s
        sy = cy - v[:, 1] * s
        self._screen = (sx, sy, v[:, 2], s)
        zmin, zmax = float(v[:, 2].min()), float(v[:, 2].max())
        zr = max(zmax - zmin, 1e-3)

        def shade(col, z):
            f = 0.55 * (1 - (z - zmin) / zr)  # far = fade towards background
            return QColor(int(col.red() + (bg.red() - col.red()) * f),
                          int(col.green() + (bg.green() - col.green()) * f),
                          int(col.blue() + (bg.blue() - col.blue()) * f))

        syms = [a.GetSymbol() for a in self.mol.GetAtoms()]
        items = []
        if self.style != "Space filling":
            for b in self.disp.GetBonds():
                i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
                if not (self.visible(i) and self.visible(j)):
                    continue
                order = {Chem.BondType.DOUBLE: 2, Chem.BondType.TRIPLE: 3}.get(b.GetBondType(), 1)
                mx, my, mz = (sx[i] + sx[j]) / 2, (sy[i] + sy[j]) / 2, (v[i, 2] + v[j, 2]) / 2
                for a, o in ((i, j), (j, i)):
                    items.append(((v[a, 2] + mz) / 2 - 0.01, 0, (a, o, mx, my, order)))
        for a in range(len(syms)):
            if self.visible(a):
                items.append((v[a, 2], 1, a))
        items.sort(key=lambda t: (t[0], t[1]))

        stick = (0.13 if self.style == "Ball & Stick" else 0.14) * s
        for z, typ, data in items:
            if typ == 0:
                a, o, mx, my, order = data
                col = shade(color_of(syms[a]), z)
                x0, y0 = sx[a], sy[a]
                if self.style == "Wireframe":
                    p.setPen(QPen(col, max(1.5, 0.06 * s), Qt.SolidLine, Qt.RoundCap))
                    p.drawLine(QPointF(x0, y0), QPointF(mx, my))
                    continue
                dx, dy = mx - x0, my - y0
                ln = math.hypot(dx, dy) or 1.0
                px, py = -dy / ln, dx / ln
                if order == 1:
                    offs, w = [0.0], stick * 2
                else:
                    gap = 0.11 * s
                    offs = [-gap / 2, gap / 2] if order == 2 else [-gap, 0, gap]
                    w = stick * (1.0 if order == 2 else 0.8)
                for off in offs:
                    a0 = QPointF(x0 + px * off, y0 + py * off)
                    a1 = QPointF(mx + px * off, my + py * off)
                    p.setPen(QPen(col.darker(150), w, Qt.SolidLine, Qt.FlatCap))
                    p.drawLine(a0, a1)
                    p.setPen(QPen(col.lighter(112), w * 0.62, Qt.SolidLine, Qt.FlatCap))
                    p.drawLine(QPointF(a0.x() - px * w * 0.06, a0.y() - py * w * 0.06),
                               QPointF(a1.x() - px * w * 0.06, a1.y() - py * w * 0.06))
            else:
                a = data
                r = self.atom_radius(syms[a]) * s
                if self.style == "Wireframe":
                    if self.mol.GetAtomWithIdx(a).GetDegree() == 0:
                        r = 0.15 * s
                    else:
                        continue
                if self.style == "Sticks":
                    r = stick
                col = shade(color_of(syms[a]), z)
                g = QRadialGradient(QPointF(sx[a] - 0.35 * r, sy[a] - 0.35 * r), 1.35 * r)
                g.setColorAt(0.0, col.lighter(175))
                g.setColorAt(0.35, col)
                g.setColorAt(1.0, col.darker(230))
                p.setBrush(QBrush(g))
                p.setPen(QPen(col.darker(250), 0.6) if self.style != "Sticks" else Qt.NoPen)
                p.drawEllipse(QPointF(sx[a], sy[a]), r, r)

        # bonds being formed / broken (dashed)
        if self.extra_bonds:
            p.setPen(QPen(QColor(230, 130, 0, 220), max(1.5, 0.05 * s), Qt.DashLine, Qt.RoundCap))
            for i, j in self.extra_bonds:
                if i < len(syms) and j < len(syms) and self.visible(i) and self.visible(j):
                    p.drawLine(QPointF(sx[i], sy[i]), QPointF(sx[j], sy[j]))

        # labels on top
        if self.labels != "none":
            f = QFont(); f.setPointSizeF(max(7.5, min(13.0, 0.28 * s))); f.setBold(True)
            p.setFont(f)
            for a, sym in enumerate(syms):
                if not self.visible(a):
                    continue
                atom = self.mol.GetAtomWithIdx(a)
                if self.labels == "index":
                    txt = str(a)
                elif self.labels == "all" or sym not in ("C", "H"):
                    txt = sym
                    ch = atom.GetFormalCharge()
                    if ch:
                        txt += (str(abs(ch)) if abs(ch) > 1 else "") + ("+" if ch > 0 else "−")
                else:
                    continue
                rect = QRectF(sx[a] - 30, sy[a] - 12, 60, 24)
                p.setPen(QColor(0, 0, 0, 160) if not self.dark else QColor(0, 0, 0, 200))
                p.drawText(rect.translated(1, 1), Qt.AlignCenter, txt)
                p.setPen(QColor(255, 255, 255))
                p.drawText(rect, Qt.AlignCenter, txt)

        # hover
        if self.hover is not None and self.hover < len(syms) and self.visible(self.hover):
            a = self.hover
            r = max(self.atom_radius(syms[a]) * s, 0.2 * s) + 3
            p.setBrush(Qt.NoBrush)
            p.setPen(QPen(QColor(0, 160, 255, 200), 2))
            p.drawEllipse(QPointF(sx[a], sy[a]), r, r)

        # selection + measurement
        if self.sel:
            p.setBrush(QColor(255, 190, 0, 70))
            numbered = len(self.sel) <= 4
            for k, a in enumerate(self.sel):
                r = max(self.atom_radius(syms[a]) * s, 0.2 * s) + 4
                p.setPen(QPen(QColor(255, 170, 0), 2.2))
                p.drawEllipse(QPointF(sx[a], sy[a]), r, r)
                if numbered:
                    p.setPen(QColor(255, 140, 0))
                    p.drawText(QPointF(sx[a] + r, sy[a] - r), str(k + 1))
            if numbered:
                p.setPen(QPen(QColor(255, 150, 0), 1.5, Qt.DashLine))
                for a, b in zip(self.sel, self.sel[1:]):
                    p.drawLine(QPointF(sx[a], sy[a]), QPointF(sx[b], sy[b]))
            m = self.measurement()
            f = QFont(); f.setPointSize(11); f.setBold(True); p.setFont(f)
            p.setPen(fg)
            p.drawText(QRectF(12, 10, self.width() - 24, 24), Qt.AlignLeft,
                       m or f"{len(self.sel)} atoms selected")

        if self._box is not None:
            p.setPen(QPen(QColor(0, 140, 255), 1, Qt.DashLine))
            p.setBrush(QColor(0, 140, 255, 35))
            p.drawRect(self._box.normalized())

        if self.trace and self.trace.get("energies"):
            self._paint_trace(p, fg, bg)

        # mode badge + hint
        f = QFont(); f.setPointSize(8); p.setFont(f)
        if self.mode == "edit":
            badge = QRectF(self.width() - 128, 10, 116, 22)
            p.setPen(Qt.NoPen); p.setBrush(QColor(255, 150, 0))
            p.drawRoundedRect(badge, 6, 6)
            p.setPen(QColor(255, 255, 255))
            fb = QFont(); fb.setPointSize(9); fb.setBold(True); p.setFont(fb)
            p.drawText(badge, Qt.AlignCenter, "EDIT ATOMS (E)")
            p.setFont(f)
            hint = ("click/box: select (Shift adds)   drag atom: move   Alt-drag: depth   "
                    "Shift-drag: rotate selection   right-drag: rotate view   dbl-click: molecule")
        else:
            hint = (f"{self.style}  ·  drag: rotate   Shift-drag: spin   right-drag: pan   wheel: zoom   "
                    "click atoms: measure   E: edit atoms")
        p.setPen(QColor(150, 150, 150))
        p.drawText(QRectF(10, self.height() - 22, self.width() - 20, 18), Qt.AlignLeft, hint)
        self._paint_highlight(p)

    def _paint_highlight(self, p):
        if self.highlight is not None:
            p.setBrush(Qt.NoBrush)
            p.setPen(QPen(self.highlight, 3))
            p.drawRect(QRectF(1.5, 1.5, self.width() - 3, self.height() - 3))

    def _paint_trace(self, p, fg, bg):
        """Energy-vs-step curve (log ΔE) with the current position of the animation."""
        tr = self.trace
        E = np.asarray(tr["energies"], float)
        if len(E) < 2:
            return
        W, H = 300, 130
        box = QRectF(self.width() - W - 12, self.height() - H - 34, W, H)
        p.setPen(QPen(QColor(150, 150, 150, 120), 1))
        p.setBrush(QColor(bg.red(), bg.green(), bg.blue(), 225))
        p.drawRoundedRect(box, 8, 8)
        plot = box.adjusted(38, 26, -10, -18)
        emin = float(np.nanmin(E))
        y = np.log10(1.0 + np.clip(E - emin, 0, None))
        ymax = max(float(np.nanmax(y)), 1e-6)
        n = len(E)
        stride = max(1, n // 400)
        idx = float(tr.get("idx", 0.0))

        def pt(i):
            return QPointF(plot.left() + plot.width() * i / (n - 1),
                           plot.bottom() - plot.height() * y[i] / ymax)
        pts = [pt(i) for i in range(0, n, stride)] + [pt(n - 1)]
        cur = int(min(n - 1, idx))
        p.setPen(QPen(QColor(140, 140, 150, 110), 1.4))
        for a, b in zip(pts, pts[1:]):
            p.drawLine(a, b)
        p.setPen(QPen(QColor(0, 150, 255), 2.2))
        done_pts = [pt(i) for i in range(0, cur + 1, stride)] + [pt(cur)]
        for a, b in zip(done_pts, done_pts[1:]):
            p.drawLine(a, b)
        c = pt(cur)
        p.setBrush(QColor(255, 120, 0)); p.setPen(Qt.NoPen)
        p.drawEllipse(c, 4.5, 4.5)
        f = QFont(); f.setPointSize(8); p.setFont(f)
        p.setPen(fg)
        e = tr.get("e", E[cur])
        method = tr.get("method", "")
        p.drawText(QRectF(box.left() + 10, box.top() + 4, W - 20, 18), Qt.AlignLeft,
                   f"{method}   step {cur}/{n - 1 if tr.get('done') else '…'}")
        p.drawText(QRectF(box.left() + 10, box.top() + 4, W - 20, 18), Qt.AlignRight,
                   core.format_energy(e, method))
        p.setPen(QColor(140, 140, 140))
        p.drawText(QRectF(box.left() + 4, plot.top() - 2, 34, 14), Qt.AlignLeft, f"{E[0] - emin:.0f}")
        p.drawText(QRectF(box.left() + 4, plot.bottom() - 10, 34, 14), Qt.AlignLeft, "0")
        p.drawText(QRectF(plot.left(), plot.bottom() + 1, plot.width(), 14), Qt.AlignCenter,
                   "ΔE to lowest (kcal/mol, log scale)  vs  step")

    # ---------------------------------------------------------------- measuring
    def measurement(self):
        if self.mol is None or not 2 <= len(self.sel) <= 4:
            return ""
        pos = self.positions()
        P = [pos[i] for i in self.sel]
        lab = "–".join(f"{self.mol.GetAtomWithIdx(i).GetSymbol()}{i}" for i in self.sel)
        if len(P) == 2:
            return f"{lab}   distance = {np.linalg.norm(P[0] - P[1]):.4f} Å"
        if len(P) == 3:
            a, b = P[0] - P[1], P[2] - P[1]
            ang = math.degrees(math.acos(np.clip(a @ b / np.linalg.norm(a) / np.linalg.norm(b), -1, 1)))
            return f"{lab}   angle = {ang:.2f}°"
        b0, b1, b2 = P[0] - P[1], P[2] - P[1], P[3] - P[2]
        b1n = b1 / np.linalg.norm(b1)
        vv = b0 - (b0 @ b1n) * b1n
        ww = b2 - (b2 @ b1n) * b1n
        x = vv @ ww
        y = np.cross(b1n, vv) @ ww
        return f"{lab}   dihedral = {math.degrees(math.atan2(y, x)):.2f}°"

    def pick(self, pt):
        if self._screen is None or self.mol is None:
            return None
        sx, sy, z, s = self._screen
        best = None
        for a in range(len(sx)):
            if not self.visible(a):
                continue
            r = max(self.atom_radius(self.mol.GetAtomWithIdx(a).GetSymbol()) * s, 0.2 * s, 6)
            if (sx[a] - pt.x()) ** 2 + (sy[a] - pt.y()) ** 2 <= r * r:
                if best is None or z[a] > z[best]:
                    best = a
        return best

    # ---------------------------------------------------------------- editing helpers
    def _to_model(self, dv):
        """Screen-frame vector (x right, y up, z to viewer, in Å) -> model frame."""
        return self.R.T @ np.asarray(dv, float)

    def moving_set(self, atoms):
        atoms = set(atoms)
        if self.move_h:
            atoms |= core.attached_hydrogens(self.mol, [a for a in atoms
                                                         if self.mol.GetAtomWithIdx(a).GetAtomicNum() > 1])
        return sorted(atoms)

    def translate_atoms(self, atoms, dv_view):
        conf = self.mol.GetConformer(self.conf_id)
        d = self._to_model(dv_view)
        for i in atoms:
            q = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (q.x + d[0], q.y + d[1], q.z + d[2]))

    def rotate_atoms(self, atoms, rot_view):
        conf = self.mol.GetConformer(self.conf_id)
        pos = conf.GetPositions()
        pivot_atoms = [a for a in self.sel if a in atoms] or atoms
        c = pos[pivot_atoms].mean(axis=0)
        M = self.R.T @ rot_view @ self.R
        for i in atoms:
            conf.SetAtomPosition(i, (c + M @ (pos[i] - c)).tolist())

    def select_fragment(self, a):
        for frag in Chem.GetMolFrags(self.mol):
            if a in frag:
                return list(frag)
        return [a]

    def _emit_sel(self):
        self.selectionChanged.emit(self.measurement())
        self.update()

    # ---------------------------------------------------------------- input
    def mousePressEvent(self, e):
        self.setFocus()
        self._press = self._last = e.position()
        self._moved = False
        btn, mods = e.button(), e.modifiers()
        hit = self.pick(e.position())
        if self.locked:  # optimisation running: only look around
            self._drag = "pan" if btn != Qt.LeftButton or mods & Qt.ControlModifier else "rotate"
            return
        if btn == Qt.MiddleButton or (btn == Qt.LeftButton and mods & Qt.ControlModifier):
            self._drag = "pan"
        elif self.mode == "edit" and btn == Qt.RightButton:
            self._drag = "rotate"
        elif btn == Qt.RightButton:
            self._drag = "pan"
        elif self.mode == "edit" and btn == Qt.LeftButton and self.mol is not None:
            if hit is None:
                self._drag = "box"
                self._box_add = bool(mods & Qt.ShiftModifier)
                self._box = QRectF(e.position(), e.position())
            else:
                if hit not in self.sel and not mods & Qt.ShiftModifier:
                    self._pending = ("replace", hit)
                else:
                    self._pending = ("toggle", hit) if mods & Qt.ShiftModifier else None
                self._hit = hit
                if mods & Qt.AltModifier:
                    self._drag = "depth"
                elif mods & Qt.ShiftModifier and hit in self.sel:
                    self._drag = "turn"
                else:
                    self._drag = "move"
        else:
            self._drag = "rotate"

    def _begin_edit(self):
        if getattr(self, "_pending", None) and self._pending[0] == "replace":
            self.sel = [self._pending[1]]
        elif self._hit not in self.sel:
            self.sel = self.sel + [self._hit]
        self._pending = None
        self._moving = self.moving_set(self.sel)
        self.editStarted.emit()

    def mouseMoveEvent(self, e):
        if self._last is None:  # hover only
            if self.mol is not None:
                h = self.pick(e.position()) if self.mode == "edit" else None
                if h != self.hover:
                    self.hover = h
                    self.setCursor(Qt.OpenHandCursor if h is not None else Qt.ArrowCursor)
                    self.update()
            return
        d = e.position() - self._last
        self._last = e.position()
        first_move = False
        if not self._moved and (e.position() - self._press).manhattanLength() > 3:
            self._moved = True
            first_move = True
        if not self._moved:
            return
        s = self._screen[3] if self._screen else 50.0
        mods = e.modifiers()
        if self._drag == "pan":
            self.pan += d
        elif self._drag == "rotate":
            if mods & Qt.ShiftModifier and self.mode != "edit":
                self.R = rot("z", -d.x() * 0.01) @ self.R
            else:
                self.R = rot("x", d.y() * 0.01) @ rot("y", d.x() * 0.01) @ self.R
        elif self._drag == "box":
            self._box = QRectF(self._press, e.position())
        elif self._drag in ("move", "depth", "turn"):
            if first_move:
                self._begin_edit()
                self.setCursor(Qt.ClosedHandCursor)
            if self._drag == "move":
                self.translate_atoms(self._moving, (d.x() / s, -d.y() / s, 0.0))
            elif self._drag == "depth":
                self.translate_atoms(self._moving, (0.0, 0.0, -d.y() / s))
            else:
                self.rotate_atoms(self._moving, rot("x", d.y() * 0.01) @ rot("y", d.x() * 0.01))
            self.selectionChanged.emit(self.measurement())
        self.update()

    def mouseReleaseEvent(self, e):
        drag, self._drag = self._drag, None
        if drag == "box":
            if self._moved and self._screen is not None:
                r = self._box.normalized()
                sx, sy, _, _ = self._screen
                inside = [a for a in range(len(sx)) if self.visible(a) and r.contains(QPointF(sx[a], sy[a]))]
                self.sel = (self.sel + [a for a in inside if a not in self.sel]) if self._box_add else inside
            elif not self._box_add:
                self.sel = []
            self._box = None
        elif drag in ("move", "depth", "turn"):
            if self._moved:
                self.setCursor(Qt.OpenHandCursor)
                self.editFinished.emit(list(self._moving))
            else:  # plain click on an atom in edit mode
                kind, a = (self._pending or ("noop", self._hit))
                if kind == "replace":
                    self.sel = [a]
                elif kind == "toggle":
                    self.sel = [x for x in self.sel if x != a] if a in self.sel else self.sel + [a]
                else:  # click on an already selected atom without Shift -> select only it
                    self.sel = [a]
            self._pending = None
        elif not self._moved and e.button() == Qt.LeftButton and self.mode == "rotate":
            a = self.pick(e.position())
            if a is None:
                self.sel = []
            elif a in self.sel:
                self.sel.remove(a)
            else:
                self.sel = (self.sel + [a])[-4:]
        self._last = None
        self._emit_sel()

    def mouseDoubleClickEvent(self, e):
        if self.mode == "edit" and self.mol is not None:
            a = self.pick(e.position())
            if a is not None:
                self.sel = self.select_fragment(a)
                self._emit_sel()

    def wheelEvent(self, e):
        self.zoom *= 1.1 ** (e.angleDelta().y() / 120)
        self.zoom = min(max(self.zoom, 0.1), 30)
        self.update()

    def keyPressEvent(self, e):
        k, a = e.key(), math.radians(5)
        mods = e.modifiers()
        if k == Qt.Key_E and not mods:
            return self.set_mode("rotate" if self.mode == "edit" else "edit")
        if k == Qt.Key_R and not mods:
            return self.set_mode("rotate")
        nudge = {Qt.Key_Left: (-1, 0, 0), Qt.Key_Right: (1, 0, 0),
                 Qt.Key_Up: (0, 1, 0), Qt.Key_Down: (0, -1, 0),
                 Qt.Key_PageUp: (0, 0, 1), Qt.Key_PageDown: (0, 0, -1)}
        if self.mode == "edit" and self.sel and k in nudge and self.mol is not None and not self.locked:
            step = 0.5 if mods & Qt.ShiftModifier else 0.05
            atoms = self.moving_set(self.sel)
            self.editStarted.emit()
            self.translate_atoms(atoms, [c * step for c in nudge[k]])
            self.editFinished.emit(atoms)
            return self._emit_sel()
        if k == Qt.Key_Left:
            self.R = rot("y", -a) @ self.R
        elif k == Qt.Key_Right:
            self.R = rot("y", a) @ self.R
        elif k == Qt.Key_Up:
            self.R = rot("x", -a) @ self.R
        elif k == Qt.Key_Down:
            self.R = rot("x", a) @ self.R
        elif k in (Qt.Key_Plus, Qt.Key_Equal):
            self.zoom *= 1.1
        elif k == Qt.Key_Minus:
            self.zoom /= 1.1
        elif k == Qt.Key_Escape:
            self.sel = []
            self.selectionChanged.emit("")
        else:
            return super().keyPressEvent(e)
        self.update()


# ======================================================================================
# background worker so the UI never freezes during embedding / conformer search
# ======================================================================================
class Worker(QThread):
    done = Signal(object)
    failed = Signal(str)

    def __init__(self, fn, *a, **kw):
        super().__init__()
        self.fn, self.a, self.kw = fn, a, kw

    def run(self):
        try:
            self.done.emit(self.fn(*self.a, **self.kw))
        except Exception as ex:  # noqa: BLE001
            traceback.print_exc()
            self.failed.emit(str(ex))


class StreamWorker(QThread):
    """Runs fn(progress, cancel) in the background; progress() items arrive as `frame` signals."""
    frame = Signal(object)
    done = Signal(object)
    failed = Signal(str)

    def __init__(self, fn):
        super().__init__()
        self.fn = fn
        self.cancel = threading.Event()

    def run(self):
        try:
            self.done.emit(self.fn(self.frame.emit, self.cancel))
        except core.Cancelled:
            self.failed.emit("")
        except Exception as ex:  # noqa: BLE001
            traceback.print_exc()
            self.failed.emit(str(ex) or type(ex).__name__)


# ======================================================================================
# xtb settings dialog
# ======================================================================================
class XtbDialog(QDialog):
    def __init__(self, main):
        super().__init__(main)
        self.main, self.st = main, main.st
        self.setWindowTitle("xtb settings (GFN2-xTB / GFN-FF)")
        self.resize(620, 330)
        cur = main.xtb_settings()
        self.path = QLineEdit(cur.path or self.st.value("xtb/path", ""))
        self.path.setPlaceholderText("not found – download it or browse to xtb / xtb.exe")
        browse = QPushButton("Browse…"); browse.clicked.connect(self.browse)
        self.dl = QPushButton(f"Download xtb {core.XTB_VERSION}")
        self.dl.clicked.connect(self.download)
        if sys.platform == "darwin":
            self.dl.setEnabled(False)
            self.dl.setToolTip("No official macOS binary: conda install -c conda-forge xtb")
        self.state = QLabel("")
        self.state.setWordWrap(True)
        self.threads = QSpinBox(); self.threads.setRange(1, max(1, os.cpu_count() or 1))
        self.threads.setValue(cur.threads)
        self.solvent = QComboBox(); self.solvent.addItems(core.ALPB_SOLVENTS)
        self.solvent.setCurrentText(cur.solvent)
        self.level = QComboBox(); self.level.addItems(["crude", "sloppy", "loose", "normal", "tight", "vtight"])
        self.level.setCurrentText(cur.level)
        self.charge = QSpinBox(); self.charge.setRange(-21, 20); self.charge.setSpecialValueText("auto")
        self.charge.setValue(-21 if cur.charge is None else cur.charge)
        self.uhf = QSpinBox(); self.uhf.setRange(-1, 20); self.uhf.setSpecialValueText("auto")
        self.uhf.setValue(-1 if cur.uhf is None else cur.uhf)

        form = QFormLayout()
        row = QHBoxLayout(); row.addWidget(self.path, 1); row.addWidget(browse); row.addWidget(self.dl)
        form.addRow("xtb program", row)
        form.addRow("", self.state)
        form.addRow("CPU threads", self.threads)
        form.addRow("Implicit solvent (ALPB)", self.solvent)
        form.addRow("Convergence (--opt)", self.level)
        form.addRow("Charge", self.charge)
        form.addRow("Unpaired electrons (--uhf)", self.uhf)
        note = QLabel("GFN2-xTB (Grimme) is a fast semi-empirical method: a few seconds for typical "
                      "organic molecules (< 100 atoms) on a normal PC. GFN-FF is xtb's force field, "
                      "faster still. 'auto' takes charge / spin from the structure. Please cite: "
                      "Bannwarth, Ehlert, Grimme, JCTC 2019, 15, 1652.")
        note.setWordWrap(True); note.setStyleSheet("color: gray")
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.save); bb.rejected.connect(self.reject)
        lay = QVBoxLayout(self)
        lay.addLayout(form); lay.addWidget(note); lay.addStretch(); lay.addWidget(bb)
        self.path.textChanged.connect(self.check)
        self.check()

    def check(self):
        p = self.path.text().strip()
        if p and os.path.isfile(p):
            self.state.setText("✓ xtb found")
            self.state.setStyleSheet("color: #2a8a2a")
        else:
            self.state.setText("xtb not found. Windows/Linux: click Download (official GitHub release, "
                               "checksum-verified). macOS: conda install -c conda-forge xtb, then Browse.")
            self.state.setStyleSheet("color: #b05000")

    def browse(self):
        p, _ = QFileDialog.getOpenFileName(self, "Locate xtb executable", os.path.expanduser("~"))
        if p:
            self.path.setText(p)

    def download(self):
        self.dl.setEnabled(False)
        dest = self.main.xtb_dir()
        w = StreamWorker(lambda progress, cancel: core.download_xtb(dest, progress))
        w.setParent(self)
        w.frame.connect(lambda m: self.state.setText(str(m)))
        w.done.connect(lambda exe: (self.path.setText(exe), self.dl.setEnabled(True)))
        w.failed.connect(lambda e: (self.state.setText("Download failed: " + e), self.dl.setEnabled(True)))
        self._w = w
        w.start()

    def save(self):
        self.st.setValue("xtb/path", self.path.text().strip())
        self.st.setValue("xtb/threads", self.threads.value())
        self.st.setValue("xtb/solvent", self.solvent.currentText())
        self.st.setValue("xtb/level", self.level.currentText())
        self.st.setValue("xtb/charge", "auto" if self.charge.value() == -21 else self.charge.value())
        self.st.setValue("xtb/uhf", "auto" if self.uhf.value() == -1 else self.uhf.value())
        self.accept()


# ======================================================================================
# Gaussian dialog
# ======================================================================================
class GaussianDialog(QDialog):
    def __init__(self, parent, mol, conf_id, settings: QSettings, title=""):
        super().__init__(parent)
        self.setWindowTitle("Gaussian input")
        self.resize(760, 640)
        self.mol, self.conf_id, self.st = mol, conf_id, settings
        ch, mu = core.default_charge_mult(mol)
        self.preset = QComboBox()
        for name, rt in core.GAUSSIAN_PRESETS:
            self.preset.addItem(name, rt)
        self.route = QLineEdit(settings.value("g16/route", core.GAUSSIAN_PRESETS[0][1]))
        self.charge = QSpinBox(); self.charge.setRange(-20, 20); self.charge.setValue(ch)
        self.mult = QSpinBox(); self.mult.setRange(1, 15); self.mult.setValue(mu)
        self.nprocs = QSpinBox(); self.nprocs.setRange(1, 512)
        self.nprocs.setValue(int(settings.value("g16/nprocs", 4)))
        self.mem = QLineEdit(settings.value("g16/mem", "4GB")); self.mem.setMaximumWidth(90)
        self.chk = QCheckBox("%chk (named after the file)")
        self.chk.setChecked(settings.value("g16/chk", True, type=bool))
        self.title = QLineEdit(title)
        self.extra = QPlainTextEdit(settings.value("g16/extra", ""))
        self.extra.setPlaceholderText("Text after the coordinates, e.g. ModRedundant lines "
                                      "(D 1 2 3 4 S 36 10.0), gen basis sets, …")
        self.extra.setMaximumHeight(90)
        self.preview = QPlainTextEdit()
        mono = QFont("Consolas" if sys.platform.startswith("win") else "Menlo" if sys.platform == "darwin" else "Monospace")
        mono.setStyleHint(QFont.Monospace)
        self.preview.setFont(mono)
        form = QFormLayout()
        form.addRow("Preset", self.preset)
        form.addRow("Route (#)", self.route)
        row = QHBoxLayout()
        for lab, w in (("Charge", self.charge), ("Multiplicity", self.mult),
                       ("%nprocshared", self.nprocs), ("%mem", self.mem)):
            row.addWidget(QLabel(lab)); row.addWidget(w)
        row.addWidget(self.chk); row.addStretch()
        form.addRow(row)
        form.addRow("Title", self.title)
        form.addRow("After coordinates", self.extra)
        bb = QDialogButtonBox()
        copy_b = bb.addButton("Copy", QDialogButtonBox.ActionRole)
        save_b = bb.addButton("Save…", QDialogButtonBox.AcceptRole)
        bb.addButton(QDialogButtonBox.Close)
        copy_b.clicked.connect(self.copy); save_b.clicked.connect(self.save)
        bb.rejected.connect(self.reject)
        lay = QVBoxLayout(self)
        lay.addLayout(form)
        lay.addWidget(QLabel("Preview (editable):"))
        lay.addWidget(self.preview, 1)
        lay.addWidget(bb)
        self.preset.currentIndexChanged.connect(
            lambda i: self.preset.itemData(i).strip() != "#p" and self.route.setText(self.preset.itemData(i)))
        for w in (self.route, self.title, self.mem):
            w.textChanged.connect(self.refresh)
        for w in (self.charge, self.mult, self.nprocs):
            w.valueChanged.connect(self.refresh)
        self.extra.textChanged.connect(self.refresh)
        self.chk.toggled.connect(self.refresh)
        self.refresh()

    def text(self, path=None):
        chk = ""
        if self.chk.isChecked():
            chk = (os.path.splitext(os.path.basename(path))[0] if path else
                   core.formula_and_mass(self.mol)[0]) + ".chk"
        return core.gaussian_input(self.mol, self.conf_id, self.route.text(), self.charge.value(),
                                   self.mult.value(), self.nprocs.value(), self.mem.text().strip(),
                                   chk, self.title.text(), self.extra.toPlainText())

    def refresh(self):
        self.preview.setPlainText(self.text())
        nelec = sum(a.GetAtomicNum() for a in self.mol.GetAtoms()) - self.charge.value()
        bad = (nelec % 2 == 0) != (self.mult.value() % 2 == 1)
        self.mult.setStyleSheet("background:#ffd6d6" if bad else "")

    def _remember(self):
        self.st.setValue("g16/route", self.route.text())
        self.st.setValue("g16/nprocs", self.nprocs.value())
        self.st.setValue("g16/mem", self.mem.text().strip())
        self.st.setValue("g16/chk", self.chk.isChecked())
        self.st.setValue("g16/extra", self.extra.toPlainText())

    def copy(self):
        clip.write_text(self.preview.toPlainText())
        self._remember()

    def save(self):
        last = self.st.value("lastdir", os.path.expanduser("~"))
        name = core.formula_and_mass(self.mol)[0] + ".gjf"
        path, _ = QFileDialog.getSaveFileName(self, "Save Gaussian input", os.path.join(last, name),
                                              "Gaussian input (*.gjf *.com);;All files (*)")
        if not path:
            return
        self.st.setValue("lastdir", os.path.dirname(path))
        txt = self.preview.toPlainText()
        if self.chk.isChecked() and txt == self.text():
            txt = self.text(path)  # %chk matches the chosen file name
        with open(path, "w", newline="\n") as fh:
            fh.write(txt if txt.endswith("\n\n") else txt.rstrip("\n") + "\n\n")
        self._remember()
        self.accept()


# ======================================================================================
# ORCA dialog
# ======================================================================================
class OrcaDialog(QDialog):
    def __init__(self, parent, mol, conf_id, settings: QSettings, title=""):
        super().__init__(parent)
        self.setWindowTitle("ORCA input")
        self.resize(760, 640)
        self.mol, self.conf_id, self.st = mol, conf_id, settings
        ch, mu = core.default_charge_mult(mol)

        self.preset = QComboBox()
        for name, kw in core.ORCA_PRESETS:
            self.preset.addItem(name, kw)
        self.keywords = QLineEdit(settings.value("orca/keywords", core.ORCA_PRESETS[0][1]))
        self.charge = QSpinBox(); self.charge.setRange(-20, 20); self.charge.setValue(ch)
        self.mult = QSpinBox(); self.mult.setRange(1, 15); self.mult.setValue(mu)
        self.nprocs = QSpinBox(); self.nprocs.setRange(1, 256)
        self.nprocs.setValue(int(settings.value("orca/nprocs", 4)))
        self.maxcore = QSpinBox(); self.maxcore.setRange(0, 1_000_000); self.maxcore.setSingleStep(500)
        self.maxcore.setValue(int(settings.value("orca/maxcore", 2000)))
        self.title = QLineEdit(title)
        self.extra = QPlainTextEdit(settings.value("orca/extra", ""))
        self.extra.setPlaceholderText("%geom ... end, %scf ... end, %cpcm smd true SMDsolvent \"water\" end ...")
        self.extra.setMaximumHeight(90)
        self.xyzfile = QCheckBox("Write coordinates to a separate .xyz (* xyzfile)")
        self.preview = QPlainTextEdit()
        mono = QFont("Consolas" if sys.platform.startswith("win") else "Menlo" if sys.platform == "darwin" else "Monospace")
        mono.setStyleHint(QFont.Monospace)
        self.preview.setFont(mono)

        form = QFormLayout()
        form.addRow("Preset", self.preset)
        form.addRow("Keywords", self.keywords)
        row = QHBoxLayout()
        for lab, w in (("Charge", self.charge), ("Multiplicity", self.mult),
                       ("nprocs", self.nprocs), ("maxcore (MB)", self.maxcore)):
            row.addWidget(QLabel(lab)); row.addWidget(w)
        row.addStretch()
        form.addRow(row)
        form.addRow("Title", self.title)
        form.addRow("Extra blocks", self.extra)
        form.addRow("", self.xyzfile)

        bb = QDialogButtonBox()
        copy_b = bb.addButton("Copy", QDialogButtonBox.ActionRole)
        save_b = bb.addButton("Save…", QDialogButtonBox.AcceptRole)
        bb.addButton(QDialogButtonBox.Close)
        copy_b.clicked.connect(self.copy)
        save_b.clicked.connect(self.save)
        bb.rejected.connect(self.reject)

        lay = QVBoxLayout(self)
        lay.addLayout(form)
        lay.addWidget(QLabel("Preview (editable):"))
        lay.addWidget(self.preview, 1)
        lay.addWidget(bb)

        self.preset.currentIndexChanged.connect(self._preset)
        for w in (self.keywords, self.title):
            w.textChanged.connect(self.refresh)
        for w in (self.charge, self.mult, self.nprocs, self.maxcore):
            w.valueChanged.connect(self.refresh)
        self.extra.textChanged.connect(self.refresh)
        self.xyzfile.toggled.connect(self.refresh)
        self.refresh()

    def _preset(self, i):
        kw = self.preset.itemData(i)
        if kw.strip() != "!":
            self.keywords.setText(kw)

    def _xyzname(self, inp_path=None):
        base = os.path.splitext(os.path.basename(inp_path))[0] if inp_path else "structure"
        return base + ".xyz"

    def text(self, inp_path=None):
        return core.orca_input(self.mol, self.conf_id, self.keywords.text(), self.charge.value(),
                               self.mult.value(), self.nprocs.value(), self.maxcore.value(),
                               self.extra.toPlainText(), self.title.text(),
                               self._xyzname(inp_path) if self.xyzfile.isChecked() else None)

    def refresh(self):
        self.preview.setPlainText(self.text())
        nelec = sum(a.GetAtomicNum() for a in self.mol.GetAtoms()) - self.charge.value()
        bad = (nelec % 2 == 0) != (self.mult.value() % 2 == 1)
        self.mult.setStyleSheet("background:#ffd6d6" if bad else "")
        self.mult.setToolTip(f"{nelec} electrons" + (" – impossible multiplicity!" if bad else ""))

    def _remember(self):
        self.st.setValue("orca/keywords", self.keywords.text())
        self.st.setValue("orca/nprocs", self.nprocs.value())
        self.st.setValue("orca/maxcore", self.maxcore.value())
        self.st.setValue("orca/extra", self.extra.toPlainText())

    def copy(self):
        clip.write_text(self.preview.toPlainText())
        self._remember()

    def save(self):
        last = self.st.value("lastdir", os.path.expanduser("~"))
        name = core.formula_and_mass(self.mol)[0] + ".inp"
        path, _ = QFileDialog.getSaveFileName(self, "Save ORCA input", os.path.join(last, name),
                                              "ORCA input (*.inp);;All files (*)")
        if not path:
            return
        self.st.setValue("lastdir", os.path.dirname(path))
        txt = self.preview.toPlainText()
        if self.xyzfile.isChecked():
            txt = self.text(path)  # make the xyzfile name match the chosen file
            xyz = os.path.join(os.path.dirname(path), self._xyzname(path))
            with open(xyz, "w", newline="\n") as fh:
                fh.write(core.xyz_block(self.mol, self.conf_id, self.title.text() or "Mol3D Bridge"))
        with open(path, "w", newline="\n") as fh:
            fh.write(txt)
        self._remember()
        self.accept()


# ======================================================================================
# main window
# ======================================================================================
class Main(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"{APP} {__version__}")
        self.resize(1180, 780)
        self.setAcceptDrops(True)
        self.st = QSettings("Mol3DBridge", "Mol3DBridge")
        self.view = MolView()
        self.setCentralWidget(self.view)
        self.mol = None          # 3D model with H
        self.source = None       # what was pasted / opened (the drawing)
        self.energies = []       # per conformer, kcal/mol
        self.emeth = []          # per conformer, method that produced the energy
        self.undo = []
        self.anim = None         # running optimisation animation
        self.timer = QTimer(self)
        self.timer.setInterval(33)
        self.timer.timeout.connect(self._tick)
        self.seed = 0xF00D
        self.worker = None
        self._build_ui()
        self.view.selectionChanged.connect(lambda m: self.statusBar().showMessage(m) if m else None)
        self.view.editStarted.connect(self.push_undo)
        self.view.editFinished.connect(self.on_edit)
        self.view.modeChanged.connect(self._mode_changed)
        self.statusBar().showMessage(
            "Ready.  RDKit ChemDraw (CDX) support: " + ("yes" if core.HAS_CDX else "no – MOL/SMILES only"))

    # ------------------------------------------------------------------ UI
    def act(self, text, slot, shortcut=None, menu=None, tip=None, checkable=False, checked=False):
        a = QAction(text, self)
        if shortcut:
            a.setShortcut(QKeySequence(shortcut))
        if tip:
            a.setStatusTip(tip); a.setToolTip(tip)
        if checkable:
            a.setCheckable(True); a.setChecked(checked)
            a.toggled.connect(slot)
        else:
            a.triggered.connect(slot)
        if menu is not None:
            menu.addAction(a)
        return a

    def _build_ui(self):
        mb = self.menuBar()
        fm = mb.addMenu("&File")
        self.act("&Open…", self.open_file, QKeySequence.Open, fm)
        self.act("Save 3D &As…  (xyz / mol / sdf / pdb)", self.save_3d, QKeySequence.Save, fm)
        self.act("&ORCA input…", self.orca, "Ctrl+E", fm)
        self.act("&Gaussian input…", self.gaussian, "Ctrl+Shift+G", fm)
        fm.addSeparator()
        self.act("&Quit", self.close, QKeySequence.Quit, fm)

        em = mb.addMenu("&Edit")
        self.a_paste = self.act("&Paste structure (ChemDraw / ReShiki)", self.paste, QKeySequence.Paste, em)
        em.addSeparator()
        self.a_copy = self.act("&Copy as 2D drawing (current view)", lambda: self.copy_2d("view"),
                               QKeySequence.Copy, em,
                               "Projection of the current orientation; wedges from real depth")
        self.act("Copy as &clean 2D drawing", lambda: self.copy_2d("clean"), "Ctrl+Shift+C", em)
        self.a_copy3d_cd = self.act("Copy 3D structure for Chem&Draw (rotatable in 3D)",
                                    lambda: self.copy_3d("chemdraw"), "Ctrl+D", em,
                                    "Like Chem3D: CDX with 3D atom positions, so ChemDraw can rotate it in 3D")
        self.a_copy3d_rs = self.act("Copy 3D structure for Re&Shiki (keeps depth)",
                                    lambda: self.copy_3d("reshiki"), "Ctrl+Shift+D", em,
                                    "3D molfile: current view as the drawing, z kept as ReShiki depth")
        self.act("Copy &3D molfile", self.copy_3d_mol, "Ctrl+Alt+C", em)
        self.act("Copy &XYZ", lambda: self.mol and clip.write_text(core.xyz_block(self.mol, self.view.conf_id)),
                 None, em)
        self.act("Copy ORCA &input (last settings)", self.copy_orca_quick, "Ctrl+Shift+E", em)
        self.act("Copy &image of 3D view", self.copy_image, "Ctrl+Shift+I", em)
        em.addSeparator()
        self.act("&Undo", self.do_undo, QKeySequence.Undo, em)
        em.addSeparator()
        mg = QActionGroup(self)
        self.a_rot = self.act("&Rotate / measure mode", lambda: self.view.set_mode("rotate"), "R", em,
                              checkable=False)
        self.a_edit = self.act("&Edit atoms mode (select & move)", lambda: self.view.set_mode("edit"), "E", em)
        for a in (self.a_rot, self.a_edit):
            a.setCheckable(True); mg.addAction(a)
        self.a_rot.setChecked(True)
        self.act("Select &all atoms", self.select_all, QKeySequence.SelectAll, em)
        self.act("Select &none", self.select_none, None, em)
        self.act("Select &molecule (of selected atoms)", self.select_molecule, None, em)
        self.act("Select &substituent (2 bonded atoms: side of the 2nd)", self.select_substituent,
                 "Ctrl+Shift+S", em, "Select atom 1 (stays) then atom 2 (moves): selects the group on atom 2")
        self.act("Set &distance / angle / dihedral…", self.set_internal, "Ctrl+G", em,
                 "Select 2, 3 or 4 atoms in order; the side of the last atom moves")
        em.addSeparator()
        self.a_moveh = self.act("Move attached hydrogens with atoms", self._moveh, None, em,
                                checkable=True, checked=self.st.value("moveh", True, type=bool))
        self.view.move_h = self.a_moveh.isChecked()
        self.a_fix = self.act("Hold selected atoms fixed when minimizing", lambda v: self.st.setValue("fix", v),
                              None, em, checkable=True, checked=self.st.value("fix", False, type=bool))
        self.a_relax = self.act("Relax the rest after each move (moved atoms held)",
                                lambda v: self.st.setValue("relax", v), None, em, checkable=True,
                                checked=self.st.value("relax", False, type=bool))
        em.addSeparator()
        self.act("Show what is on the clipboard", lambda: QMessageBox.information(
            self, APP, "Readable clipboard formats:\n\n" + clip.available_summary()), None, em)

        cm = mb.addMenu("&Calculate")
        self.a_min = self.act("&Minimize from current geometry (incl. your edits)", self.minimize, "Ctrl+M", cm,
                              "Minimization (selected method) starting from the structure as it is now")
        self.a_stop = self.act("&Stop / skip animation", self.stop_anim, "Esc", cm)
        self.a_stop.setEnabled(False)
        self.a_opt = self.act("Re&build 3D from the drawing (not minimized)", self.preoptimize, "Ctrl+B", cm,
                              "Discard 3D edits: build the pasted drawing in 3D again (no minimization)")
        self.act("&New random 3D start (not minimized)", self.reembed, "Ctrl+R", cm)
        self.a_conf = self.act("&Conformer search…", self.conformer_search, "Ctrl+K", cm)
        self.a_tsw = self.act("&Transition-state search (pysisyphus)…", self.open_ts, "Ctrl+T", cm,
                              "Reactant + product → growing string → TS → IRC → minima")
        cm.addSeparator()
        self.a_refine = self.act("Refine all conformers with xtb", self.refine_all, None, cm,
                                 "Re-optimize every conformer with GFN2-xTB / GFN-FF and re-rank")
        cm.addSeparator()
        ffm = cm.addMenu("&Method")
        self.ffgroup = QActionGroup(self)
        cur = self.st.value("ff", "MMFF94")
        if cur not in core.METHODS:
            cur = "MMFF94"
        for ff in core.METHODS:
            a = QAction(ff, self, checkable=True)
            a.setChecked(ff == cur)
            self.ffgroup.addAction(a); ffm.addAction(a)
            if ff == core.FORCE_FIELDS[-1]:
                ffm.addSeparator()
        self.ffgroup.triggered.connect(lambda a: self.ffcombo.setCurrentText(a.text()))
        self.act("&xtb settings…", self.xtb_dialog, None, cm)

        vm = mb.addMenu("&View")
        sg = QActionGroup(self)
        for s_ in MolView.STYLES:
            a = QAction(s_, self, checkable=True)
            a.setChecked(s_ == self.view.style)
            a.triggered.connect(lambda _=False, s_=s_: self.set_style(s_))
            sg.addAction(a); vm.addAction(a)
        vm.addSeparator()
        self.act("Show &hydrogens", self.toggle_h, "Ctrl+H", vm, checkable=True, checked=True)
        lg = QActionGroup(self)
        lm = vm.addMenu("&Labels")
        for key, name in (("none", "None"), ("hetero", "Heteroatoms"), ("all", "All elements"),
                          ("index", "Atom index (0-based, as ORCA)")):
            a = QAction(name, self, checkable=True)
            a.setChecked(key == self.view.labels)
            a.triggered.connect(lambda _=False, k=key: self.set_labels(k))
            lg.addAction(a); lm.addAction(a)
        vm.addSeparator()
        self.act("&Best 2D orientation", self.best_view, "Ctrl+0", vm,
                 "Lay the molecule flat in the screen plane (good before copying back)")
        self.act("&Reset zoom / pan", self.reset_zoom, None, vm)
        self.act("&Dark background", self.toggle_dark, None, vm, checkable=True, checked=False)

        sm = mb.addMenu("&Settings")
        self.a_auto = self.act("Minimize automatically after pasting", lambda v: self.st.setValue("auto_min", v),
                               None, sm, checkable=True, checked=self.st.value("auto_min", False, type=bool))
        self.a_keeph = self.act("Keep hydrogens in copied 2D drawing", lambda v: self.st.setValue("keeph", v),
                                None, sm, checkable=True, checked=self.st.value("keeph", False, type=bool))
        self.a_cdx = self.act("Put ChemDraw CDX on the clipboard when copying", lambda v: self.st.setValue("cdx", v),
                              None, sm, checkable=True,
                              checked=self.st.value("cdx", True, type=bool) and core.HAS_CDX)
        self.a_cdx.setEnabled(core.HAS_CDX)
        self.a_textmol = self.act("Plain-text copy = molfile (instead of SMILES)",
                                  lambda v: self.st.setValue("textmol", v), None, sm, checkable=True,
                                  checked=self.st.value("textmol", False, type=bool))

        sm.addSeparator()
        self.a_anim = self.act("Animate optimizations (watch it fall into the minimum)",
                               lambda v: self.st.setValue("anim", v), None, sm, checkable=True,
                               checked=self.st.value("anim", True, type=bool))
        self.act("Animation length…", self.set_anim_len, None, sm)
        self.act("Force-field step size…", self.set_step, None, sm)

        sm.addSeparator()
        self.act("ChemDraw bond length for copies…", self.set_cd_bond, None, sm,
                 "Match your ChemDraw document style (ACS 1996 = 14.4 pt, ChemDraw default = 30 pt)")
        self.a_cdx3d = self.act("ChemDraw 3D copy: use CDX 3D positions (off = 3D molfile)",
                                lambda v: self.st.setValue("cdx3d", v), None, sm, checkable=True,
                                checked=self.st.value("cdx3d", True, type=bool))
        self.a_zflip = self.act("ChemDraw 3D copy: invert depth (only if ChemDraw shows the mirror image after a 3D rotation)",
                                lambda v: self.st.setValue("zflip", v), None, sm, checkable=True,
                                checked=self.st.value("zflip", False, type=bool))

        hm = mb.addMenu("&Help")
        self.act("&How to use", self.help, QKeySequence.HelpContents, hm)
        self.act("Open &log / data folder", lambda: QDesktopServices.openUrl(
            QUrl.fromLocalFile(user_data_dir())), None, hm)
        hm.addSeparator()
        self.act(f"&About {APP}", self.about, None, hm)

        tb = QToolBar("Main")
        tb.setMovable(False)
        tb.setToolButtonStyle(Qt.ToolButtonTextOnly)
        self.addToolBar(tb)
        tb.addAction(self.a_paste); self.a_paste.setIconText("Paste ▸ 3D")
        tb.addSeparator()
        tb.addAction(self.a_rot); self.a_rot.setIconText("⟳ Rotate")
        tb.addAction(self.a_edit); self.a_edit.setIconText("✥ Edit atoms")
        tb.addSeparator()
        tb.addAction(self.a_min); self.a_min.setIconText("Minimize")
        tb.addAction(self.a_stop); self.a_stop.setIconText("■ Stop")
        tb.addAction(self.a_opt); self.a_opt.setIconText("Rebuild")
        tb.addAction(self.a_conf); self.a_conf.setIconText("Conformers…")
        tb.addAction(self.a_tsw); self.a_tsw.setIconText("TS search…")
        tb.addSeparator()
        tb.addWidget(QLabel(" Method "))
        self.ffcombo = QComboBox(); self.ffcombo.addItems(core.METHODS)
        self.ffcombo.setCurrentText(cur)
        self.ffcombo.currentTextChanged.connect(self._ff_changed)
        tb.addWidget(self.ffcombo)
        tb.addWidget(QLabel("  Conformer "))
        self.confcombo = QComboBox(); self.confcombo.setMinimumWidth(200)
        self.confcombo.currentIndexChanged.connect(self.pick_conf)
        tb.addWidget(self.confcombo)
        tb.addSeparator()
        tb.addAction(self.a_copy); self.a_copy.setIconText("Copy ▸ 2D")
        b3 = QToolButton()
        b3.setText("Copy ▸ 3D")
        b3.setToolTip("Copy keeping the 3D coordinates (Ctrl+D ChemDraw, Ctrl+Shift+D ReShiki)")
        m3 = QMenu(b3)
        m3.addAction(self.a_copy3d_cd)
        m3.addAction(self.a_copy3d_rs)
        b3.setMenu(m3)
        b3.setPopupMode(QToolButton.InstantPopup)
        tb.addWidget(b3)
        b = QAction("ORCA…", self); b.triggered.connect(self.orca); tb.addAction(b)
        b = QAction("Gaussian…", self); b.triggered.connect(self.gaussian); tb.addAction(b)
        tb.addSeparator()
        self.smiles = QLineEdit(); self.smiles.setPlaceholderText("SMILES / InChI  ↵")
        self.smiles.setMaximumWidth(260)
        self.smiles.returnPressed.connect(self.from_smiles)
        tb.addWidget(self.smiles)

        self.info = QLabel("")
        self.statusBar().addPermanentWidget(self.info)

    # ------------------------------------------------------------------ helpers
    def ff(self):
        """Selected method (MMFF94 / MMFF94s / UFF / GFN2-xTB / GFN-FF)."""
        return self.ffcombo.currentText()

    def search_ff(self):
        """Force field used for embedding-based conformer search."""
        return self.ff() if not core.is_xtb(self.ff()) else "MMFF94"

    def _ensure_e(self):
        n = self.mol.GetNumConformers() if self.mol is not None else 0
        self.energies = (list(self.energies) + [float("nan")] * n)[:n]
        self.emeth = (list(self.emeth) + [""] * n)[:n]

    def set_energy(self, cid, e, method):
        self._ensure_e()
        self.energies[cid], self.emeth[cid] = e, method

    def _ff_changed(self, t):
        self.st.setValue("ff", t)
        for a in self.ffgroup.actions():
            a.setChecked(a.text() == t)

    def busy(self, on, msg=""):
        for a in (self.a_paste, self.a_opt, self.a_min, self.a_conf, self.a_refine):
            a.setEnabled(not on)
        if on:
            QApplication.setOverrideCursor(Qt.WaitCursor)
            self.statusBar().showMessage(msg)
        else:
            QApplication.restoreOverrideCursor()

    def run(self, fn, done, msg, *a, **kw):
        if self.worker and self.worker.isRunning():
            return
        self.busy(True, msg)
        self.worker = Worker(fn, *a, **kw)
        self.worker.setParent(self)
        self.worker.done.connect(lambda r: (self.busy(False), done(r)))
        self.worker.failed.connect(lambda e: (self.busy(False), self.error(e)))
        self.worker.start()

    def closeEvent(self, e):
        if self.anim and self.anim.get("cancel"):
            self.anim["cancel"].set()
        if self.worker is not None and self.worker.isRunning():
            if getattr(self.worker, "cancel", None):
                self.worker.cancel.set()
            self.worker.wait(5000)
        super().closeEvent(e)

    def error(self, msg):
        self.statusBar().showMessage("Error: " + msg)
        QMessageBox.warning(self, APP, msg)

    def push_undo(self):
        if self.mol is not None:
            self.undo.append((Chem.Mol(self.mol), list(self.energies), list(self.emeth), self.view.conf_id))
            self.undo = self.undo[-30:]

    def do_undo(self):
        if not self.undo or self.anim:
            return
        self.mol, self.energies, self.emeth, cid = self.undo.pop()
        self.view.trace = None
        self.show_model(keep_view=True, conf=cid)
        self.statusBar().showMessage("Undone.")

    def need_mol(self):
        if self.mol is None:
            self.statusBar().showMessage("No structure yet – paste one from ChemDraw / ReShiki (Ctrl+V).")
            return False
        return True

    def show_model(self, keep_view=False, conf=0):
        self.confcombo.blockSignals(True)
        self.confcombo.clear()
        self._ensure_e()
        n = self.mol.GetNumConformers()
        for i in range(n):
            e, m = self.energies[i], self.emeth[i]
            if math.isnan(e) or not m:
                self.confcombo.addItem(f"#{i + 1}")
                continue
            same = [x for x, mm in zip(self.energies, self.emeth) if mm == m and not math.isnan(x)]
            de = e - min(same)
            self.confcombo.addItem(f"#{i + 1}  {m}  ΔE = {de:.2f} kcal/mol")
        self.confcombo.setCurrentIndex(min(conf, n - 1))
        self.confcombo.blockSignals(False)
        self.view.set_molecule(self.mol, min(conf, n - 1), keep_view=keep_view)
        self.update_info()

    def update_info(self):
        if self.mol is None:
            self.info.setText("")
            return
        f, mw = core.formula_and_mass(self.mol)
        cid = self.view.conf_id
        self._ensure_e()
        e, m = self.energies[cid], self.emeth[cid]
        es = f"   {m}: {core.format_energy(e, m)}" if m and not math.isnan(e) else ""
        self.info.setText(f"{f}   MW {mw:.2f}   {self.mol.GetNumAtoms()} atoms{es}")

    def pick_conf(self, i):
        if self.mol is not None and i >= 0:
            self.view.conf_id = i
            self.view.sel = []
            self.view.recenter()
            self.view.update()
            self.update_info()

    # ------------------------------------------------------------------ input
    def accept_structure(self, mol, origin):
        if mol is None or mol.GetNumAtoms() == 0:
            self.error(f"Could not read a structure from {origin}.")
            return
        if self.anim is not None or (self.worker and self.worker.isRunning()):
            self.statusBar().showMessage("Busy – wait for the current calculation (or press Esc) and paste again.")
            return
        mol = Chem.Mol(mol)
        try:
            Chem.SanitizeMol(mol)
        except Exception as ex:  # noqa: BLE001
            self.error(f"The structure from {origin} is not chemically valid:\n{ex}\n\n"
                       "Check valences / charges in the drawing.")
            return
        if core.normalize_scale(mol):
            pass  # coordinates were in other units (points, 1/65536 pt, ...)
        self.push_undo()
        self.source = Chem.Mol(mol)
        if core.is_3d(mol):
            # already 3D (Chem3D / ReShiki 3D / xyz / sdf): keep the coordinates as they are
            m = Chem.AddHs(mol, addCoords=True)
            self.mol, self.energies, self.emeth = m, [], []
            self.view.trace = None
            self.show_model()
            self.statusBar().showMessage(f"Loaded 3D structure from {origin}. Ctrl+M to minimize.")
            return
        self.run(core.build_3d, lambda r: self._embedded(r, origin),
                 "Building 3D model (this can take a few seconds for strained rings)…", mol, None, 1, self.seed)

    def _embedded(self, r: core.BuildResult, origin=""):
        """3D embedding done -> show it, then (optionally) the animated pre-optimisation."""
        self.mol, self.energies, self.emeth = r.mol, [], []
        self.view.trace = None
        self.show_model()
        self.view.R = np.eye(3)   # the model was aligned to the drawing: show it the same way
        self.view.update()
        msg = f"3D model from {origin}, not minimized." if origin else "New 3D start, not minimized."
        msg += "  Press Ctrl+M (Minimize) when ready."
        if r.messages:
            msg += "  " + "  ".join(r.messages)
        self.statusBar().showMessage(msg)
        if self.a_auto.isChecked():
            self.minimize(push=False)

    def _built(self, r: core.BuildResult, origin="", keep_view=False):
        self.mol = r.mol
        self.energies, self.emeth = list(r.energies), [r.forcefield] * len(r.energies)
        self.view.trace = None
        self.show_model(keep_view=keep_view)
        e = r.energies[0] if r.energies else float("nan")
        msg = f"3D model from {origin}: " if origin else ""
        msg += f"{r.forcefield} E = {e:.3f} kcal/mol" if r.forcefield else "embedded (no FF)"
        if r.mol.GetNumConformers() > 1:
            msg += f"   ·   {r.mol.GetNumConformers()} conformers"
        if r.messages:
            msg += "   ·   " + "  ".join(r.messages)
        self.statusBar().showMessage(msg)

    def paste(self):
        items = clip.read_items()
        if not items:
            self.error("The clipboard has no structure.\n\nIn ChemDraw/ReShiki select the structure "
                       "and Copy (Ctrl+C), then paste here.")
            return
        errors, parsed = [], []
        for kind, data, src in items:
            try:
                m = core.parse_clipboard_item(kind, data)
            except Exception as ex:  # noqa: BLE001
                errors.append(f"{src}: {ex}")
                continue
            if m is not None and m.GetNumAtoms():
                parsed.append((kind, m))
            else:
                errors.append(f"{src}: not a structure")
        if parsed:
            # real 3D coordinates (e.g. a 3D molfile next to the CDX) win over a 2D drawing
            kind, m = next(((k, x) for k, x in parsed if core.is_3d(x)), parsed[0])
            self.accept_structure(m, f"clipboard ({kind})")
            return
        self.error("Could not read the clipboard:\n" + "\n".join(errors) +
                   "\n\nTip: in ChemDraw use Edit ▸ Copy As ▸ MOL Text or SMILES.")

    def from_smiles(self):
        t = self.smiles.text().strip()
        if not t:
            return
        m = core.parse_text(t)
        if m is None:
            self.error("Could not parse that SMILES/InChI.")
            return
        self.accept_structure(m, "SMILES")

    def open_file(self):
        last = self.st.value("lastdir", os.path.expanduser("~"))
        path, _ = QFileDialog.getOpenFileName(
            self, "Open structure", last,
            "Structures (*.cdx *.cdxml *.mol *.sdf *.xyz *.pdb *.smi *.txt);;All files (*)")
        if path:
            self.load_path(path)

    def load_path(self, path):
        self.st.setValue("lastdir", os.path.dirname(path))
        try:
            m = core.load_file(path)
        except Exception as ex:  # noqa: BLE001
            self.error(str(ex))
            return
        self.accept_structure(m, os.path.basename(path))

    def dragEnterEvent(self, e):
        if e.mimeData().hasUrls():
            e.acceptProposedAction()

    def dropEvent(self, e):
        for u in e.mimeData().urls():
            if u.isLocalFile():
                self.load_path(u.toLocalFile())
                break

    # ------------------------------------------------------------------ calculations
    def preoptimize(self):
        if self.source is None:
            return self.need_mol()
        self.push_undo()
        self.run(core.build_3d, self._embedded, "Embedding…", self.source, None, 1, self.seed, force_embed=True)

    def reembed(self):
        if self.source is None:
            return self.need_mol()
        self.push_undo()
        self.seed = random.randint(1, 2**31 - 1)
        self.run(core.build_3d, self._embedded, "Re-embedding…", self.source, None, 1, self.seed,
                 force_embed=True)

    def minimize(self, _checked=False, fixed_override=None, push=True, animate=None):
        """Minimize from the current geometry with the selected method, step by step."""
        if not self.need_mol() or self.anim is not None:
            return
        if self.worker and self.worker.isRunning():
            return
        method = self.ff()
        cid = self.view.conf_id
        mol = Chem.Mol(self.mol)
        fixed = list(self.view.sel) if (self.a_fix.isChecked() and fixed_override is None) else fixed_override
        animate = self.a_anim.isChecked() if animate is None else animate
        if core.is_xtb(method):
            st = self.xtb_settings()
            if not st.path:
                if QMessageBox.question(
                        self, APP, f"{method} needs the xtb program, which was not found.\n\n"
                        "Set it up now? (download ~60 MB on Windows/Linux, or choose an installed xtb)"
                        ) == QMessageBox.Yes:
                    self.xtb_dialog()
                return

            def fn(progress, cancel):
                return core.xtb_optimize(mol, cid, method, st, fixed, progress, cancel)
        else:
            step = float(self.st.value("step", 0.02))

            def fn(progress, cancel):
                return core.ff_trajectory(mol, cid, method, fixed, max_step=step, cancel=cancel)
        if push:
            self.push_undo()
        self.anim = dict(frames=[], idx=0.0, done=False, result=None, cid=cid, method=method,
                         t0=time.monotonic(), animate=animate, fixed=fixed,
                         T=float(self.st.value("animlen", 5.0)))
        self.view.locked = True
        self.view.trace = dict(energies=[], idx=0.0, method=method, done=False)
        self.busy(True, f"Optimizing with {method}…")
        QApplication.restoreOverrideCursor()  # keep the normal cursor: the user can rotate meanwhile
        self.a_stop.setEnabled(True)
        w = StreamWorker(fn)
        w.setParent(self)
        self.anim["cancel"] = w.cancel
        w.frame.connect(self._frame)
        w.done.connect(self._opt_done)
        w.failed.connect(self._opt_failed)
        self.worker = w
        w.start()
        if animate:
            self.timer.start()

    def _frame(self, fr):
        if self.anim is not None and isinstance(fr, tuple):
            self.anim["frames"].append(fr)
            self.view.trace["energies"].append(fr[1])

    def _opt_done(self, r: core.OptResult):
        if self.anim is None:
            return
        self.anim.update(frames=r.frames, done=True, result=r)
        self.view.trace["energies"] = [f[1] for f in r.frames]
        self.view.trace["done"] = True
        if not self.anim["animate"] or r.cancelled:
            self._finish_anim()

    def _opt_failed(self, msg):
        self.timer.stop()
        a, self.anim = self.anim, None
        self.view.locked = False
        self.view.trace = None
        self.busy(False)
        self.a_stop.setEnabled(False)
        if self.undo and a is not None:
            self.mol, self.energies, self.emeth, cid = self.undo.pop()
            self.show_model(keep_view=True, conf=cid)
        if "Cancelled" in msg or not msg:
            self.statusBar().showMessage("Stopped.")
        else:
            self.error(msg)

    def _show_frame_at(self, idx):
        fr = self.anim["frames"]
        i0 = int(math.floor(idx))
        i1 = min(i0 + 1, len(fr) - 1)
        t = idx - i0
        pos = fr[i0][0] * (1 - t) + fr[i1][0] * t
        core.apply_positions(self.mol, self.anim["cid"], pos)
        self.view.trace["idx"] = idx
        e = fr[i0][1] * (1 - t) + fr[i1][1] * t
        self.view.trace["e"] = e
        self.view.update()

    def _tick(self):
        """Advance the animation by wall-clock time, so it lasts ~T seconds on any PC."""
        a = self.anim
        if a is None:
            self.timer.stop()
            return
        n = len(a["frames"])
        if n == 0:
            return
        now = time.monotonic()
        if a["done"]:
            if "t_done" not in a:
                a["t_done"], a["idx_done"] = now, a["idx"]
                fr = a["frames"]
                rmsd = float(np.sqrt(((fr[-1][0] - fr[0][0]) ** 2).sum(axis=1).mean()))
                T = a["T"] if rmsd > 0.05 else min(a["T"], 1.0)  # nothing much to show
                a["t_end"] = max(a["t0"] + T, now + 0.4)
            span = max(a["t_end"] - a["t_done"], 1e-3)
            u = min(1.0, (now - a["t_done"]) / span)
            a["idx"] = a["idx_done"] + (n - 1 - a["idx_done"]) * u
            self._show_frame_at(min(a["idx"], n - 1.0))
            if u >= 1.0:
                self._finish_anim()
        else:  # still computing (xtb): follow the incoming steps smoothly
            target = n - 1.0
            a["idx"] = min(target, a["idx"] + max(0.2, (target - a["idx"]) * 0.12))
            self._show_frame_at(a["idx"])

    def stop_anim(self):
        a = self.anim
        if a is None:
            return
        if not a["done"] or a.get("result") is None:
            a["cancel"].set()           # worker returns what it has
            self.statusBar().showMessage("Stopping…")
        else:
            self._finish_anim()         # skip the rest of the animation

    def _finish_anim(self):
        self.timer.stop()
        a, self.anim = self.anim, None
        r: core.OptResult = a["result"]
        self.view.locked = False
        self.busy(False)
        self.a_stop.setEnabled(False)
        cid = a["cid"]
        core.apply_positions(self.mol, cid, r.final)
        self.set_energy(cid, r.energy, r.method)
        if self.view.trace is not None:
            self.view.trace.update(idx=len(r.frames) - 1.0, e=r.energy, done=True)
        self.show_model(keep_view=True, conf=cid)
        warn = core.stereo_changes(self.source, self.mol, cid)
        steps = len(r.frames) - 1
        msg = (f"{r.method}: E = {core.format_energy(r.energy, r.method)}  ·  {steps} steps  " +
               " ".join(r.messages) + (f"   ⚠ {warn}" if warn else ""))
        self.statusBar().showMessage(msg)

    def cd_bond(self) -> float:
        return float(self.st.value("cdbond", 14.4))

    def set_cd_bond(self):
        items = ["14.4 pt – ACS Document 1996", "30 pt – ChemDraw New Document", "Custom…"]
        cur = self.cd_bond()
        idx = 0 if abs(cur - 14.4) < 1e-6 else 1 if abs(cur - 30) < 1e-6 else 2
        it, ok = QInputDialog.getItem(self, "ChemDraw bond length",
                                      "Bond length of copied drawings (use your ChemDraw document "
                                      "style so ChemDraw doesn't ask about settings):", items, idx, False)
        if not ok:
            return
        if it.startswith("14.4"):
            v = 14.4
        elif it.startswith("30"):
            v = 30.0
        else:
            v, ok = QInputDialog.getDouble(self, "ChemDraw bond length", "Bond length (pt):", cur, 2.0, 200.0, 2)
            if not ok:
                return
        self.st.setValue("cdbond", v)
        self.statusBar().showMessage(f"ChemDraw copies now use a {v:g} pt bond length.")

    def set_anim_len(self):
        v, ok = QInputDialog.getDouble(self, "Animation length",
                                       "Seconds an optimization animation should take (3–8 recommended):",
                                       float(self.st.value("animlen", 5.0)), 1.0, 20.0, 1)
        if ok:
            self.st.setValue("animlen", v)

    def set_step(self):
        v, ok = QInputDialog.getDouble(self, "Force-field step size",
                                       "Largest move of any atom per optimization step (Å).\n"
                                       "Smaller = more, smoother steps (the result is the same):",
                                       float(self.st.value("step", 0.02)), 0.002, 0.2, 3)
        if ok:
            self.st.setValue("step", v)

    # ------------------------------------------------------------------ xtb
    def xtb_dir(self):
        return os.path.join(user_data_dir(), "xtb")

    def xtb_settings(self) -> core.XtbSettings:
        path = core.find_xtb(self.st.value("xtb/path", ""),
                             [self.xtb_dir(), os.path.join(app_base_dir(), "xtb")]) or ""
        ch = self.st.value("xtb/charge", "auto")
        uhf = self.st.value("xtb/uhf", "auto")
        return core.XtbSettings(
            path=path, solvent=self.st.value("xtb/solvent", "none"),
            level=self.st.value("xtb/level", "normal"),
            threads=int(self.st.value("xtb/threads", min(4, os.cpu_count() or 1))),
            charge=None if ch in ("auto", None, "") else int(ch),
            uhf=None if uhf in ("auto", None, "") else int(uhf))

    def xtb_dialog(self):
        XtbDialog(self).exec()

    def refine_all(self):
        if not self.need_mol():
            return
        method = self.ff() if core.is_xtb(self.ff()) else "GFN2-xTB"
        st = self.xtb_settings()
        if not st.path:
            return self.xtb_dialog()
        self.push_undo()
        mol = Chem.Mol(self.mol)

        def fn(progress, cancel):
            n = mol.GetNumConformers()
            ens, msgs = [], []
            for k, c in enumerate(list(mol.GetConformers())):
                if cancel.is_set():
                    raise core.Cancelled()
                progress(f"{method}: conformer {k + 1}/{n}…")
                r = core.xtb_optimize(mol, c.GetId(), method, st, cancel=cancel)
                core.apply_positions(mol, c.GetId(), r.final)
                ens.append(r.energy)
                if not r.converged:
                    msgs.append(f"#{k + 1} not converged")
            order = sorted(range(n), key=lambda i: ens[i])
            confs = [Chem.Conformer(mol.GetConformer(list(mol.GetConformers())[i].GetId())) for i in order]
            mol.RemoveAllConformers()
            for i, c in enumerate(confs):
                c.SetId(i)
                mol.AddConformer(c, assignId=False)
            return mol, [ens[i] for i in order], msgs

        def done(res):
            self.timer.stop()
            self.busy(False)
            self.a_stop.setEnabled(False)
            self.anim = None
            m, ens, msgs = res
            self.mol, self.energies, self.emeth = m, ens, [method] * len(ens)
            self.view.trace = None
            self.show_model(keep_view=True, conf=0)
            self.statusBar().showMessage(f"Refined {len(ens)} conformers with {method}; sorted by energy. "
                                         + " ".join(msgs))

        w = StreamWorker(fn)
        w.setParent(self)
        w.frame.connect(lambda m: isinstance(m, str) and self.statusBar().showMessage(m))
        w.done.connect(done)
        w.failed.connect(lambda e: (self.busy(False), self.a_stop.setEnabled(False),
                                    setattr(self, "anim", None),
                                    self.statusBar().showMessage("Stopped.") if "Cancelled" in e or not e
                                    else self.error(e)))
        self.worker = w
        self.anim = dict(done=False, cancel=w.cancel, frames=[])  # lets Stop cancel the batch
        self.busy(True, f"Refining conformers with {method}…")
        self.a_stop.setEnabled(True)
        w.start()

    # ------------------------------------------------------------------ 3D editing
    def _mode_changed(self, mode):
        (self.a_edit if mode == "edit" else self.a_rot).setChecked(True)
        if mode == "edit":
            self.statusBar().showMessage(
                "Edit atoms: click / box-select atoms, drag to move (Alt = depth, Shift = rotate group), "
                "right-drag rotates the view. Ctrl+M minimizes from the edited geometry.")

    def _moveh(self, v):
        self.st.setValue("moveh", v)
        self.view.move_h = v

    def on_edit(self, moved):
        cid = self.view.conf_id
        self.view.trace = None
        ff = self.ff() if not core.is_xtb(self.ff()) else "MMFF94"
        try:
            e = core.single_point(self.mol, cid, ff)
        except Exception:  # noqa: BLE001
            e = float("nan")
        self.set_energy(cid, e, ff if not math.isnan(e) else "")
        self.update_info()
        warn = core.stereo_changes(self.source, self.mol, cid)
        if self.a_relax.isChecked():
            self.minimize(fixed_override=list(moved), push=False, animate=False)
            return
        self.statusBar().showMessage(
            f"Moved {len(moved)} atom(s)" + (f"  ·  E = {e:.2f} kcal/mol (not minimized)" if not math.isnan(e) else "")
            + "  ·  Ctrl+M: minimize from here" + (f"   ⚠ {warn}" if warn else ""))

    def select_all(self):
        if self.need_mol():
            self.view.sel = list(range(self.mol.GetNumAtoms()))
            self.view._emit_sel()

    def select_none(self):
        self.view.sel = []
        self.view._emit_sel()

    def select_molecule(self):
        if self.need_mol() and self.view.sel:
            out = []
            for a in self.view.sel:
                out += [x for x in self.view.select_fragment(a) if x not in out]
            self.view.sel = out
            self.view._emit_sel()

    def select_substituent(self):
        if not self.need_mol():
            return
        if len(self.view.sel) != 2:
            return self.error("Select exactly two bonded atoms: first the atom that stays, "
                              "then the atom whose group should be selected.")
        a, b = self.view.sel
        side = core.side_of_bond(self.mol, a, b)
        if not side:
            return self.error("Those atoms are not bonded, or the bond is in a ring.")
        self.view.sel = side
        self.view._emit_sel()
        self.statusBar().showMessage(f"Selected {len(side)} atoms - drag them (edit mode) or Shift-drag to rotate.")

    def set_internal(self):
        if not self.need_mol():
            return
        sel = list(self.view.sel)
        if not 2 <= len(sel) <= 4:
            return self.error("Select 2 (distance), 3 (angle) or 4 (dihedral) atoms in order.\n"
                              "The part of the molecule on the side of the last atom is moved.")
        from rdkit.Chem import rdMolTransforms as T
        conf = self.mol.GetConformer(self.view.conf_id)
        if len(sel) == 2:
            cur, name, lo, hi, unit = T.GetBondLength(conf, *sel), "Distance", 0.5, 20.0, "Å"
        elif len(sel) == 3:
            cur, name, lo, hi, unit = T.GetAngleDeg(conf, *sel), "Angle", 1.0, 180.0, "°"
        else:
            cur, name, lo, hi, unit = T.GetDihedralDeg(conf, *sel), "Dihedral", -180.0, 180.0, "°"
        val, ok = QInputDialog.getDouble(self, f"Set {name.lower()}", f"{name} ({unit}):", cur, lo, hi, 2)
        if not ok:
            return
        self.push_undo()
        try:
            core.set_internal(self.mol, self.view.conf_id, sel, val)
        except Exception as ex:  # noqa: BLE001
            self.undo.pop()
            return self.error(f"Could not set {name.lower()}: {ex}\n"
                              "(For angles both bonds must exist; ring bonds can't be rotated.)")
        self.view.update()
        self.on_edit(sel)
        self.view._emit_sel()

    def conformer_search(self):
        if self.source is None:
            return self.need_mol()
        n, ok = QInputDialog.getInt(self, "Conformer search",
                                    "Number of ETKDG starting geometries\n(each is minimized; "
                                    "duplicates within 0.25 Å RMSD are pruned):", 30, 2, 2000)
        if not ok:
            return
        self.push_undo()
        def done(r):
            self._built(r)
            if core.is_xtb(self.ff()) and r.mol.GetNumConformers() > 1:
                if QMessageBox.question(self, APP, f"Conformer search used {r.forcefield}.\n"
                                        f"Refine all {r.mol.GetNumConformers()} conformers with "
                                        f"{self.ff()} and re-rank?") == QMessageBox.Yes:
                    self.refine_all()
        self.run(core.build_3d, done, f"Conformer search ({n} starts, {self.search_ff()})…",
                 self.source, self.search_ff(), n, self.seed, force_embed=True)

    # ------------------------------------------------------------------ output
    def copy_2d(self, mode):
        if not self.need_mol():
            return
        try:
            d = core.drawing_from_3d(self.mol, self.view.conf_id, self.view.R, mode,
                                     keep_h=self.a_keeph.isChecked(), cdx_bond_pt=self.cd_bond())
        except Exception as ex:  # noqa: BLE001
            traceback.print_exc()
            return self.error(str(ex))
        text = d.molblock if self.a_textmol.isChecked() else d.smiles
        where = clip.write_items(cdx=d.cdx if self.a_cdx.isChecked() else None,
                                 molfile=d.molblock, smiles=d.smiles, text=text)
        fmts = ("CDX + " if (d.cdx and self.a_cdx.isChecked()) else "") + "MOL + SMILES"
        msg = f"Copied {'current view' if mode == 'view' else 'clean 2D'} drawing to the {where} ({fmts})." \
              " Paste in ChemDraw / ReShiki."
        if d.warnings:
            msg += "   ⚠ " + " ".join(d.warnings)
        self.statusBar().showMessage(msg)

    def copy_3d(self, target):
        """Copy the current view as a drawing that keeps the 3D coordinates."""
        if not self.need_mol():
            return
        try:
            d = core.drawing_3d_from_view(self.mol, self.view.conf_id, self.view.R,
                                          keep_h=self.a_keeph.isChecked(),
                                          z_sign=-1.0 if self.a_zflip.isChecked() else 1.0,
                                          cdx_bond_pt=self.cd_bond())
        except Exception as ex:  # noqa: BLE001
            traceback.print_exc()
            return self.error(str(ex))
        text = d.molblock3d if self.a_textmol.isChecked() else d.smiles
        if target == "chemdraw" and self.a_cdx3d.isChecked() and d.cdx3d:
            clip.write_items(cdx=d.cdx3d, molfile=d.molblock3d_wedged, smiles=d.smiles, text=text)
            msg = ("Copied 3D structure for ChemDraw (CDX with 3D positions, like Chem3D). Paste in "
                   "ChemDraw - 3D rotation there keeps the geometry.")
        elif target == "chemdraw":
            # A 3D molfile (MDLCT): ChemDraw converts the Å coordinates itself, so no unit guessing.
            # No CDX on purpose - ChemDraw would read that first.
            clip.write_items(cdx=None, molfile=d.molblock3d_wedged, smiles=d.smiles, text=text)
            msg = ("Copied 3D structure for ChemDraw (3D molfile with wedges). Paste in ChemDraw; "
                   "the current view is the drawing and the z coordinates are kept.")
        else:
            # no CDX on purpose: ReShiki prefers CDX but keeps depth only from a 3D molfile
            clip.write_items(cdx=None, molfile=d.molblock3d, smiles=d.smiles, text=text)
            msg = ("Copied 3D structure for ReShiki (3D molfile, depth kept). Paste in ReShiki; its 3D "
                   "tool starts from this geometry.")
        if d.warnings:
            msg += "   ⚠ " + " ".join(d.warnings)
        self.statusBar().showMessage(msg)

    def copy_3d_mol(self):
        if not self.need_mol():
            return
        mb = core.molblock_3d(self.mol, self.view.conf_id)
        clip.write_items(molfile=mb, text=mb)
        self.statusBar().showMessage("Copied 3D molfile (with H).")

    def copy_image(self):
        if self.mol is not None:
            QApplication.clipboard().setImage(self.view.grab().toImage())
            self.statusBar().showMessage("Copied picture of the 3D view.")

    def orca(self):
        if not self.need_mol():
            return
        OrcaDialog(self, self.mol, self.view.conf_id, self.st,
                   core.formula_and_mass(self.mol)[0]).exec()

    def open_ts(self):
        from .ts_window import TSWindow
        if getattr(self, "ts_win", None) is None:
            self.ts_win = TSWindow(self)
        self.ts_win.show()
        self.ts_win.raise_()
        self.ts_win.activateWindow()

    def load_external_3d(self, mol, origin):
        """Show a ready 3D structure (e.g. a TS) as it is – no sanitizing, no rebuilding."""
        self.push_undo()
        self.source = None
        self.mol, self.energies, self.emeth = mol, [], []
        self.view.trace = None
        self.show_model()
        self.statusBar().showMessage(f"Loaded {origin}.")

    def gaussian(self):
        if not self.need_mol():
            return
        GaussianDialog(self, self.mol, self.view.conf_id, self.st,
                       core.formula_and_mass(self.mol)[0]).exec()

    def copy_orca_quick(self):
        if not self.need_mol():
            return
        ch, mu = core.default_charge_mult(self.mol)
        clip.write_text(core.orca_input(
            self.mol, self.view.conf_id, self.st.value("orca/keywords", core.ORCA_PRESETS[0][1]),
            ch, mu, int(self.st.value("orca/nprocs", 4)), int(self.st.value("orca/maxcore", 2000)),
            self.st.value("orca/extra", "")))
        self.statusBar().showMessage(f"Copied ORCA input (charge {ch}, multiplicity {mu}).")

    def save_3d(self):
        if not self.need_mol():
            return
        last = self.st.value("lastdir", os.path.expanduser("~"))
        name = os.path.join(last, core.formula_and_mass(self.mol)[0])
        path, flt = QFileDialog.getSaveFileName(
            self, "Save 3D coordinates", name,
            "XYZ (*.xyz);;MOL (*.mol);;SDF – all conformers (*.sdf);;PDB (*.pdb);;ORCA input (*.inp);;"
            "Gaussian input (*.gjf *.com)")
        if not path:
            return
        self.st.setValue("lastdir", os.path.dirname(path))
        ext = os.path.splitext(path)[1].lower()
        if not ext:
            ext = {"XYZ": ".xyz", "MOL": ".mol", "SDF": ".sdf", "PDB": ".pdb", "ORCA": ".inp",
                   "Gaussian": ".gjf"}[flt.split()[0]]
            path += ext
        cid = self.view.conf_id
        try:
            if ext == ".xyz":
                self._ensure_e()
                e, m = self.energies[cid], self.emeth[cid]
                comment = f"{core.formula_and_mass(self.mol)[0]}" + \
                          (f"  {m} E={core.format_energy(e, m)}" if m and not math.isnan(e) else "")
                txt = core.xyz_block(self.mol, cid, comment)
            elif ext == ".mol":
                txt = core.molblock_3d(self.mol, cid)
            elif ext == ".sdf":
                self._ensure_e()
                core.write_sdf(path, self.mol, self.energies, ",".join(sorted(set(m for m in self.emeth if m))))
                txt = None
            elif ext == ".pdb":
                txt = Chem.MolToPDBBlock(self.mol, confId=cid)
            elif ext == ".inp":
                ch, mu = core.default_charge_mult(self.mol)
                txt = core.orca_input(self.mol, cid, self.st.value("orca/keywords", core.ORCA_PRESETS[0][1]),
                                      ch, mu, int(self.st.value("orca/nprocs", 4)),
                                      int(self.st.value("orca/maxcore", 2000)), self.st.value("orca/extra", ""))
            elif ext in (".gjf", ".com"):
                ch, mu = core.default_charge_mult(self.mol)
                txt = core.gaussian_input(
                    self.mol, cid, self.st.value("g16/route", core.GAUSSIAN_PRESETS[0][1]), ch, mu,
                    int(self.st.value("g16/nprocs", 4)), self.st.value("g16/mem", "4GB"),
                    os.path.splitext(os.path.basename(path))[0] + ".chk"
                    if self.st.value("g16/chk", True, type=bool) else "",
                    "", self.st.value("g16/extra", ""))
            else:
                return self.error("Unknown extension " + ext)
            if txt is not None:
                with open(path, "w", newline="\n") as fh:
                    fh.write(txt)
            self.statusBar().showMessage("Saved " + path)
        except Exception as ex:  # noqa: BLE001
            self.error(str(ex))

    # ------------------------------------------------------------------ view
    def set_style(self, s_):
        self.view.style = s_; self.view.update()

    def set_labels(self, k):
        self.view.labels = k; self.view.update()

    def toggle_h(self, v):
        self.view.show_h = v; self.view.update()

    def toggle_dark(self, v):
        self.view.dark = v; self.view.update()

    def best_view(self):
        if self.mol is not None:
            self.view.R = core.best_view_rotation(self.mol, self.view.conf_id)
            self.view.update()

    def reset_zoom(self):
        self.view.zoom, self.view.pan = 1.0, QPointF(0, 0)
        self.view.update()

    def about(self):
        import rdkit
        from PySide6 import __version__ as pyside_version
        xtb = self.xtb_settings().path or "not installed"
        QMessageBox.about(self, f"About {APP}", (
            f"<h3>{APP} {__version__}</h3><p>© 2026 Takuma Miyamura · MIT licence</p>"
            "<p>A Chem3D-style 3D companion for ChemDraw and ReShiki: paste a drawing, build and "
            "pre-optimise it in 3D (MMFF94 / UFF / GFN2-xTB / GFN-FF), adjust it by hand, and copy it "
            "back as a 2D or rotatable 3D structure, or write ORCA / Gaussian inputs.</p>"
            f"<p>RDKit {rdkit.__version__} · Qt for Python {pyside_version} · "
            f"Python {sys.version.split()[0]}<br>xtb: {xtb}</p>"
            "<p>Please cite RDKit (rdkit.org) and, when using xtb, Bannwarth, Ehlert, Grimme, "
            "<i>J. Chem. Theory Comput.</i> 2019, 15, 1652.</p>"
            "<p>Third-party licences: THIRD_PARTY_NOTICES.md and the licenses folder in the "
            "installation folder. ChemDraw and Chem3D are trademarks of Revvity, Inc.; Mol3D Bridge is "
            "not affiliated with or endorsed by Revvity.</p>"))

    def help(self):
        QMessageBox.information(self, APP, __doc__.strip() + (
            "\n\nChemDraw ▸ here: select the structure in ChemDraw or ReShiki and Copy, then Ctrl+V here."
            "\nHere ▸ ChemDraw: rotate the model to the view you like, Ctrl+C, then paste in ChemDraw/ReShiki."
            "\n\nIf ChemDraw's normal copy can't be read, use ChemDraw ▸ Edit ▸ Copy As ▸ MOL Text."))


def app_icon():
    for name in ("icon.png", "icon_256.png"):
        path = os.path.join(RES, name)
        if os.path.exists(path):
            return QIcon(path)
    return _drawn_icon()


def _drawn_icon():
    pm = QPixmap(64, 64); pm.fill(Qt.transparent)
    p = QPainter(pm); p.setRenderHint(QPainter.Antialiasing)
    for (x, y, c) in ((20, 40, (125, 125, 125)), (44, 40, (240, 20, 20)), (32, 18, (48, 80, 248))):
        g = QRadialGradient(QPointF(x - 4, y - 4), 16)
        g.setColorAt(0, QColor(*c).lighter(170)); g.setColorAt(1, QColor(*c).darker(200))
        p.setBrush(g); p.setPen(Qt.NoPen); p.drawEllipse(QPointF(x, y), 12, 12)
    p.end()
    return QIcon(pm)


def _setup_logging():
    """Windowed builds have no console: send output and crashes to a log file."""
    import logging
    log_path = os.path.join(user_data_dir(), "mol3d_bridge.log")
    if sys.stderr is None or sys.stdout is None or getattr(sys, "frozen", False):
        fh = open(log_path, "a", encoding="utf-8", buffering=1)
        sys.stdout = sys.stdout or fh
        sys.stderr = fh
    logging.basicConfig(level=logging.INFO, stream=sys.stderr,
                        format="%(asctime)s %(levelname)s %(message)s")

    def hook(kind, value, tb):
        text = "".join(traceback.format_exception(kind, value, tb))
        try:
            sys.stderr.write(text)
        except Exception:  # noqa: BLE001
            pass
        app = QApplication.instance()
        if app is not None:
            QMessageBox.critical(None, APP, "An unexpected error occurred. Your structure is still "
                                 f"open.\n\n{value}\n\nDetails were written to:\n{log_path}")
    sys.excepthook = hook
    return log_path


def selftest() -> int:
    """Headless check used by the build pipeline: chemistry, clipboard formats, GUI start-up."""
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    app = QApplication.instance() or QApplication(sys.argv[:1])
    m = core.parse_text("C[C@H](N)C(=O)O")
    r = core.build_3d(m)
    ff = core.ff_trajectory(r.mol, 0, "MMFF94")
    core.apply_positions(r.mol, 0, ff.final)
    d2 = core.drawing_from_3d(r.mol, 0, np.eye(3))
    d3 = core.drawing_3d_from_view(r.mol, 0, np.eye(3))
    assert d2.smiles == Chem.MolToSmiles(m), "2D round trip"
    assert d3.molblock3d and (d3.cdx3d or not core.HAS_CDX), "3D copy"
    from .api import Molecule3D
    api_m = Molecule3D.from_smiles("C[C@H](N)C(=O)O")
    api_m.minimize("MMFF94")
    assert api_m.stereo_ok() and "* xyz" in api_m.orca_input(), "Python API"
    core.orca_input(r.mol, 0)
    core.gaussian_input(r.mol, 0)
    from . import ts as tsmod
    mp = tsmod.map_atoms(core.build_3d(Chem.MolFromSmiles("C=CO"), None).mol,
                         core.build_3d(Chem.MolFromSmiles("O=CC"), None).mol)
    assert len(mp.broken) + len(mp.formed) == 2, "TS atom mapping"
    w = Main()
    w.show()
    w.open_ts()
    app.processEvents()
    w.ts_win.close()
    w.close()
    msg = (f"{APP} {__version__} selftest OK  (RDKit CDX support: {core.HAS_CDX}, "
           f"xtb: {w.xtb_settings().path or 'none'})")
    xtb = w.xtb_settings()
    if xtb.path:
        r2 = core.xtb_optimize(r.mol, 0, "GFN2-xTB", xtb)
        msg += f"\n  GFN2-xTB check: {core.format_energy(r2.energy, 'GFN2-xTB')}"
    out = os.environ.get("MOL3D_SELFTEST_OUT")
    if out:
        with open(out, "w") as fh:
            fh.write(msg + "\n")
    if sys.stdout is not None:
        print(msg)
    return 0


def main():
    if "--selftest" in sys.argv:
        sys.exit(selftest())
    if sys.platform.startswith("win"):
        try:  # own taskbar icon/group instead of python.exe
            import ctypes
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(APP_ID)
        except Exception:  # noqa: BLE001
            pass
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    app.setApplicationName("Mol3D Bridge")
    app.setOrganizationName("Mol3DBridge")
    app.setApplicationVersion(__version__)
    _setup_logging()
    app.setWindowIcon(app_icon())
    w = Main()
    w.show()
    files = [a for a in sys.argv[1:] if os.path.exists(a)]
    if files:
        w.load_path(files[0])
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
