"""Transition-state search window (pysisyphus) for Mol3D Bridge."""
from __future__ import annotations

import datetime
import math
import os
import time

import numpy as np
from PySide6.QtCore import QEvent, QObject, QPointF, QRectF, QStandardPaths, Qt, QTimer, QUrl, Signal
from PySide6.QtGui import QAction, QColor, QDesktopServices, QFont, QKeySequence, QPainter, QPen
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox,
                               QFileDialog, QFrame, QFormLayout, QGridLayout, QGroupBox, QHBoxLayout,
                               QInputDialog, QLabel, QLineEdit, QMainWindow, QMessageBox,
                               QPlainTextEdit, QPushButton, QSpinBox, QSplitter, QTabWidget,
                               QToolBar, QVBoxLayout, QWidget)
from rdkit import Chem

from . import clipboard as clip
from . import core, ts
from .app import GaussianDialog, MolView, OrcaDialog, StreamWorker, user_data_dir

MONO = QFont("Consolas" if os.name == "nt" else "Menlo")
MONO.setStyleHint(QFont.Monospace)


def single_conf(mol: Chem.Mol, conf_id: int = -1) -> Chem.Mol:
    m = Chem.Mol(mol, confId=conf_id) if conf_id >= 0 else Chem.Mol(mol)
    m.GetConformer().SetId(0)
    return m


# ======================================================================================
# one structure (reactant or product)
# ======================================================================================
ACCENT = "#1e88e5"
PANEL_QSS = f"""
QFrame#panel {{ border: 1px solid palette(mid); border-radius: 8px; }}
QFrame#panel[active="true"] {{ border: 3px solid {ACCENT}; }}
QLabel#head {{ border-radius: 5px; padding: 3px 8px; }}
QLabel#head[active="true"] {{ background: {ACCENT}; color: white; }}
"""


class _Activator(QObject):
    """Any click or focus inside a panel makes it the active one."""

    def __init__(self, panel):
        super().__init__(panel)
        self.panel = panel

    def eventFilter(self, obj, ev):
        if ev.type() in (QEvent.MouseButtonPress, QEvent.FocusIn):
            self.panel.win.set_active(self.panel)
        return False


class StructurePanel(QFrame):
    changed = Signal()

    def __init__(self, win, role: str):
        super().__init__()
        self.setObjectName("panel")
        self.win, self.role = win, role
        self.mol = None
        self.energy, self.method = float("nan"), ""
        self.worker = None

        self.title = QLabel(f"<b>{role}</b>")
        self.title.setObjectName("head")
        self.info = QLabel("empty – paste, open or type a SMILES")
        self.info.setStyleSheet("color: gray")
        self.view = MolView()
        self.view.setMinimumSize(300, 260)
        self.view.empty_text = (f"{role}\n\nPaste from ChemDraw / ReShiki, open a file,\n"
                                "type a SMILES, or take the main window’s structure.")
        self.smiles = QLineEdit()
        self.smiles.setPlaceholderText("SMILES ↵")
        self.smiles.returnPressed.connect(lambda: self.from_text(self.smiles.text()))

        def btn(text, slot, tip=""):
            b = QPushButton(text)
            b.clicked.connect(slot)
            b.setToolTip(tip)
            return b

        self.b_edit = QPushButton("Edit atoms")
        self.b_edit.setCheckable(True)
        self.b_edit.setToolTip("Drag atoms to adjust the geometry (right-drag rotates)")
        self.b_edit.toggled.connect(lambda on: self.view.set_mode("edit" if on else "rotate"))
        self.b_opt = btn("Pre-optimize", self.preoptimize, "Minimize with the method chosen in the toolbar")
        other = "Product" if role == "Reactant" else "Reactant"
        arrow = "→" if role == "Reactant" else "←"
        row1 = QHBoxLayout()
        for w in (btn("Copy 3D", self.copy3d, "Copy this 3D structure (with hydrogens and numbering) – "
                      "paste it into the other side, the main window, ChemDraw or ReShiki  (Ctrl+C)"),
                  btn("Paste", self.paste, "3D structure copied from the other side, or a structure from "
                      "ChemDraw / ReShiki / the main window  (Ctrl+V)"),
                  btn(f"{arrow} Copy to {other.lower()}" if role == "Reactant" else f"{arrow} Copy to {other.lower()}",
                      self.copy_to_other, f"Put this exact 3D structure into the {other.lower()} panel"),
                  btn("From main window", self.from_main), btn("Open…", self.open_file)):
            row1.addWidget(w)
        row1.addWidget(self.smiles, 1)
        row2 = QHBoxLayout()
        row2.addWidget(self.b_opt)
        row2.addWidget(self.b_edit)
        row2.addWidget(btn("Reset view", self.reset_view))
        row2.addStretch()
        head = QHBoxLayout()
        head.addWidget(self.title)
        head.addWidget(self.info, 1)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.addLayout(head)
        lay.addLayout(row1)
        lay.addWidget(self.view, 1)
        lay.addLayout(row2)
        self.view.selectionChanged.connect(self._sel)
        self.view.editFinished.connect(lambda _: (self._energy_reset(), self.changed.emit()))
        self._activator = _Activator(self)
        for wdg in [self] + self.findChildren(QWidget):
            wdg.installEventFilter(self._activator)
        self.set_active(False)

    def set_active(self, on: bool):
        self.active = on
        for wdg in (self, self.title):
            wdg.setProperty("active", "true" if on else "false")
            wdg.style().unpolish(wdg)
            wdg.style().polish(wdg)
        self.title.setText(f"<b>{self.role}</b>" + ("  ●" if on else ""))
        self.view.highlight = QColor(ACCENT) if on else None
        self.view.update()

    # ------------------------------------------------------------------ copy between panels
    def copy3d(self):
        if self.mol is None:
            return self.win.error(f"The {self.role.lower()} panel is empty.")
        mb = Chem.MolToMolBlock(self.mol)
        m = Chem.Mol(self.mol)
        Chem.AssignStereochemistryFrom3D(m)
        smi = core.stereo_smiles(m)
        clip.write_items(molfile=mb, smiles=smi, text=mb)
        self.win.copied = (mb.strip(), Chem.Mol(self.mol), self.role)
        self.win.statusBar().showMessage(f"{self.role} copied (3D, with hydrogens and atom numbers). "
                                         "Click the other side and Paste (Ctrl+V).")

    def copy_to_other(self):
        if self.mol is None:
            return self.win.error(f"The {self.role.lower()} panel is empty.")
        other = self.win.prod if self is self.win.react else self.win.react
        other.set_mol(Chem.Mol(self.mol), note=f"copied from the {self.role.lower()}")
        other.energy, other.method = self.energy, self.method
        other.update_info(f"copied from the {self.role.lower()}")
        other.view.R = self.view.R.copy()
        other.view.update()
        self.win.set_active(other)

    # ------------------------------------------------------------------ structure in
    def set_mol(self, mol: Chem.Mol, keep_view=False, note=""):
        self.mol = single_conf(mol)
        self.view.set_molecule(self.mol, 0, keep_view=keep_view)
        self.view.labels = "index" if self.win.a_numbers.isChecked() else "hetero"
        self._energy_reset(keep=bool(keep_view and self.method))
        self.update_info(note)
        self.changed.emit()

    def _energy_reset(self, keep=False):
        if not keep:
            self.energy, self.method = float("nan"), ""
        self.update_info()

    def update_info(self, note=""):
        if self.mol is None:
            return
        f = core.formula_and_mass(self.mol)[0]
        e = f"  ·  {self.method}: {core.format_energy(self.energy, self.method)}" \
            if self.method and not math.isnan(self.energy) else "  ·  not minimized"
        self.info.setText(f"{f}, {self.mol.GetNumAtoms()} atoms{e}" + (f"  ·  {note}" if note else ""))

    def _accept(self, mol, origin):
        if mol is None or mol.GetNumAtoms() == 0:
            return self.win.error(f"{self.role}: could not read a structure from {origin}.")
        mol = Chem.Mol(mol)
        try:
            Chem.SanitizeMol(mol)
        except Exception as ex:  # noqa: BLE001
            return self.win.error(f"{self.role}: invalid structure ({ex}).")
        core.normalize_scale(mol)
        if core.is_3d(mol):
            return self.set_mol(Chem.AddHs(mol, addCoords=True), note=f"3D from {origin}")
        self._run(lambda progress, cancel: core.build_3d(mol, None, 1, 0xF00D),
                  lambda r: self.set_mol(r.mol, note=f"built from {origin}"), "Building 3D…")

    def paste(self):
        cp = getattr(self.win, "copied", None)
        if cp:
            for kind, data, _src in clip.read_items():
                if kind in ("molfile", "text"):
                    txt = core.molfile_text_from_bytes(data) if kind == "molfile" else data.decode("utf-8", "replace")
                    if txt and txt.strip() == cp[0]:
                        self.set_mol(Chem.Mol(cp[1]), note=f"pasted from the {cp[2].lower()}")
                        return
        parsed = []
        for kind, data, _src in clip.read_items():
            try:
                m = core.parse_clipboard_item(kind, data)
            except Exception:  # noqa: BLE001
                continue
            if m is not None and m.GetNumAtoms():
                parsed.append((kind, m))
        if not parsed:
            return self.win.error("The clipboard has no structure. Copy one in ChemDraw / ReShiki first.")
        kind, m = next(((k, x) for k, x in parsed if core.is_3d(x)), parsed[0])
        self._accept(m, f"the clipboard ({kind})")

    def from_text(self, text):
        if text.strip():
            self._accept(core.parse_text(text), "SMILES")

    def from_main(self):
        main = self.win.main
        if main.mol is None:
            return self.win.error("The main window has no structure.")
        self.set_mol(single_conf(main.mol, main.view.conf_id), note="from the main window")

    def open_file(self):
        path, _ = QFileDialog.getOpenFileName(self, f"Open {self.role.lower()}", self.win.lastdir(),
                                              "Structures (*.xyz *.mol *.sdf *.cdx *.cdxml *.pdb);;All files (*)")
        if path:
            self.win.st.setValue("lastdir", os.path.dirname(path))
            try:
                self._accept(core.load_file(path), os.path.basename(path))
            except Exception as ex:  # noqa: BLE001
                self.win.error(str(ex))

    # ------------------------------------------------------------------ calculations
    def _run(self, fn, done, msg):
        if self.worker and self.worker.isRunning():
            return
        self.b_opt.setEnabled(False)
        self.info.setText(msg)
        w = StreamWorker(fn)
        w.setParent(self)
        w.done.connect(lambda r: (self.b_opt.setEnabled(True), done(r)))
        w.failed.connect(lambda e: (self.b_opt.setEnabled(True), self.update_info(),
                                    e and self.win.error(f"{self.role}: {e}")))
        self.worker = w
        w.start()

    def preoptimize(self):
        if self.mol is None:
            return self.win.error(f"Load a {self.role.lower()} first.")
        method = self.win.method.currentText()
        mol = Chem.Mol(self.mol)
        fixed = list(self.view.sel) if self.win.a_fix.isChecked() else None
        if core.is_xtb(method):
            st = self.win.main.xtb_settings()
            if not st.path:
                return self.win.error("xtb was not found – set it up in the main window "
                                      "(Calculate ▸ xtb settings…).")

            def fn(progress, cancel):
                return core.xtb_optimize(mol, 0, method, st, fixed, cancel=cancel)
        else:
            def fn(progress, cancel):
                return core.ff_trajectory(mol, 0, method, fixed, cancel=cancel)

        def done(r: core.OptResult):
            core.apply_positions(self.mol, 0, r.final)
            self.energy, self.method = r.energy, r.method
            self.view.update()
            self.update_info(" ".join(r.messages[:1]))
            self.changed.emit()
        self._run(fn, done, f"Pre-optimizing with {method}…")

    # ------------------------------------------------------------------ numbering
    def _sel(self, _msg):
        if self.mol is None or not self.win.a_relabel.isChecked() or len(self.view.sel) != 1:
            return
        i = self.view.sel[0]
        sym = self.mol.GetAtomWithIdx(i).GetSymbol()
        QTimer.singleShot(0, lambda: self.win.relabel(self, i, sym))

    def renumber(self, perm, note=""):
        self.set_mol(ts.renumber(self.mol, perm), keep_view=True, note=note)

    def reset_view(self):
        if self.mol is not None:
            self.view.R = core.best_view_rotation(self.mol, 0)
            self.view.zoom, self.view.pan = 1.0, QPointF(0, 0)
            self.view.recenter()
            self.view.update()


# ======================================================================================
# energy profile of the current string / IRC
# ======================================================================================
class ProfilePlot(QWidget):
    def __init__(self):
        super().__init__()
        self.setMinimumHeight(110)
        self.values, self.cur, self.title = [], -1, ""

    def set_data(self, values_kj, title=""):
        self.values, self.title = [v for v in values_kj], title
        self.update()

    def set_current(self, i):
        self.cur = i
        self.update()

    def paintEvent(self, _):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        pal = self.palette()
        p.fillRect(self.rect(), pal.base())
        fg = pal.text().color()
        p.setPen(fg)
        p.drawText(QRectF(8, 4, self.width() - 16, 16), Qt.AlignLeft, self.title)
        vals = [v for v in self.values if v is not None]
        if len(vals) < 2:
            p.drawText(self.rect(), Qt.AlignCenter, "energy profile appears here")
            return
        lo, hi = min(vals), max(vals)
        rng = max(hi - lo, 1e-6)
        x0, y0, w, h = 40, 24, self.width() - 52, self.height() - 40
        p.setPen(QColor(150, 150, 150))
        p.drawText(QRectF(2, y0 - 6, 36, 14), Qt.AlignRight, f"{hi - lo:.0f}")
        p.drawText(QRectF(2, y0 + h - 8, 36, 14), Qt.AlignRight, "0")
        p.drawText(QRectF(x0, y0 + h + 2, w, 14), Qt.AlignCenter, "kJ/mol relative to the lowest point")
        n = len(self.values)
        pts = [QPointF(x0 + w * k / (n - 1), y0 + h - h * ((v if v is not None else lo) - lo) / rng)
               for k, v in enumerate(self.values)]
        p.setPen(QPen(QColor(0, 140, 200), 2))
        for a, b in zip(pts, pts[1:]):
            p.drawLine(a, b)
        imax = max(range(n), key=lambda k: self.values[k] if self.values[k] is not None else -1e9)
        for k, q in enumerate(pts):
            p.setBrush(QColor(230, 120, 0) if k == imax else QColor(0, 140, 200))
            p.setPen(Qt.NoPen)
            r = 5 if k == self.cur else 3
            p.drawEllipse(q, r, r)
        if 0 <= self.cur < n:
            p.setPen(QPen(QColor(230, 150, 0), 1, Qt.DashLine))
            p.drawLine(QPointF(pts[self.cur].x(), y0), QPointF(pts[self.cur].x(), y0 + h))


# ======================================================================================
# one result structure
# ======================================================================================
class ResultPane(QWidget):
    def __init__(self, win, name):
        super().__init__()
        self.win, self.name = win, name
        self.mol = None
        self.view = MolView()
        self.view.empty_text = "Run a TS search (or load a previous run)."
        self.info = QLabel("")
        self.info.setWordWrap(True)
        self.play = QPushButton("▶ Play IRC")
        self.play.setCheckable(True)
        self.play.setVisible(False)
        self.play.toggled.connect(self._play)
        self.frames, self.k = [], 0
        self.timer = QTimer(self)
        self.timer.setInterval(70)
        self.timer.timeout.connect(self._tick)
        row = QHBoxLayout()
        for text, slot in (("Copy 2D", lambda: self.copy("2d")),
                           ("Copy 3D for ChemDraw", lambda: self.copy("chemdraw")),
                           ("Copy 3D for ReShiki", lambda: self.copy("reshiki")),
                           ("Save XYZ…", self.save), ("ORCA…", self.orca), ("Gaussian…", self.gaussian),
                           ("Open in main window", self.to_main)):
            b = QPushButton(text)
            b.clicked.connect(slot)
            row.addWidget(b)
        row.addWidget(self.play)
        row.addStretch()
        lay = QVBoxLayout(self)
        lay.addWidget(self.info)
        lay.addWidget(self.view, 1)
        lay.addLayout(row)

    def set_result(self, mol, text, frames=None):
        self.mol = mol
        self.info.setText(text)
        if mol is not None:
            self.view.set_molecule(mol, 0)
            self.view.labels = "index" if self.win.a_numbers.isChecked() else "hetero"
        self.frames = frames or []
        self.play.setVisible(bool(self.frames))

    def _play(self, on):
        if on and self.frames and self.mol is not None:
            self.k = 0
            self.timer.start()
            self.play.setText("■ Stop")
        else:
            self.timer.stop()
            self.play.setText("▶ Play IRC")
            if self.mol is not None and self._ts_pos is not None:
                core.apply_positions(self.mol, 0, self._ts_pos)
                self.view.update()

    _ts_pos = None

    def _tick(self):
        if self._ts_pos is None:
            self._ts_pos = self.mol.GetConformer().GetPositions()
        n = len(self.frames)
        k = self.k % (2 * n - 2) if n > 1 else 0
        k = k if k < n else 2 * n - 2 - k                 # back and forth
        core.apply_positions(self.mol, 0, self.frames[k])
        self.view.update()
        self.k += 1

    def _need(self):
        if self.mol is None:
            self.win.error("No structure here yet.")
            return False
        return True

    def copy(self, what):
        if not self._need():
            return
        main = self.win.main
        try:
            if what == "2d":
                d = core.drawing_from_3d(self.mol, 0, self.view.R, "view", cdx_bond_pt=main.cd_bond())
                clip.write_items(cdx=d.cdx, molfile=d.molblock, smiles=d.smiles, text=d.smiles)
            else:
                d = core.drawing_3d_from_view(self.mol, 0, self.view.R, cdx_bond_pt=main.cd_bond())
                if what == "chemdraw":
                    clip.write_items(cdx=d.cdx3d if main.a_cdx3d.isChecked() else None,
                                     molfile=d.molblock3d_wedged, smiles=d.smiles, text=d.smiles)
                else:
                    clip.write_items(molfile=d.molblock3d, smiles=d.smiles, text=d.smiles)
            self.win.statusBar().showMessage(f"{self.name}: copied.")
        except Exception as ex:  # noqa: BLE001
            self.win.error(str(ex))

    def save(self):
        if not self._need():
            return
        base = os.path.join(self.win.lastdir(), self.name.split()[0].lower() + ".xyz")
        path, _ = QFileDialog.getSaveFileName(self, f"Save {self.name}", base, "XYZ (*.xyz)")
        if path:
            with open(path, "w", newline="\n") as fh:
                fh.write(core.xyz_block(self.mol, 0, self.name))

    def orca(self):
        if self._need():
            OrcaDialog(self, self.mol, 0, self.win.st, self.name).exec()

    def gaussian(self):
        if self._need():
            GaussianDialog(self, self.mol, 0, self.win.st, self.name).exec()

    def to_main(self):
        if self._need():
            self.win.main.load_external_3d(Chem.Mol(self.mol), self.name)
            self.win.main.raise_()
            self.win.main.activateWindow()


# ======================================================================================
# pysisyphus program location
# ======================================================================================
class PysisDialog(QDialog):
    def __init__(self, win):
        super().__init__(win)
        self.win = win
        self.worker = None
        self.setWindowTitle("pysisyphus program")
        self.resize(720, 520)
        self.path = QLineEdit(win.pysis_path() or "")
        self.path.setPlaceholderText("not found – install it below or browse to pysis / pysis.exe")
        browse = QPushButton("Browse…")
        browse.clicked.connect(self.browse)
        row = QHBoxLayout()
        row.addWidget(self.path, 1)
        row.addWidget(browse)
        self.state = QLabel("")
        self.state.setWordWrap(True)
        self.install_btn = QPushButton(f"Install pysisyphus {ts.PYSIS_VERSION} automatically")
        self.install_btn.setStyleSheet("font-weight: 600; padding: 6px 12px")
        self.install_btn.clicked.connect(self.install)
        self.cancel_btn = QPushButton("Cancel installation")
        self.cancel_btn.setEnabled(False)
        self.cancel_btn.clicked.connect(lambda: self.worker and self.worker.cancel.set())
        note = QLabel(
            "The automatic installation sets up a <b>private</b> Python 3.12 with pysisyphus inside "
            "Mol3D Bridge’s data folder – nothing else on the PC is changed, no administrator rights are "
            "needed, and an existing Python is left alone. It downloads about 300 MB once and uses "
            "about 1.5 GB of disk space (pysisyphus needs SciPy, JAX, Numba, …).<br>"
            "Already have pysisyphus (<code>pip install pysisyphus</code>)? It is found automatically, "
            "or choose its <code>pysis</code> program with Browse.<br>"
            "pysisyphus is free software under the GNU GPL v3 by J. Steinmetzer et al. "
            "(github.com/eljost/pysisyphus); it is downloaded from PyPI onto this PC and runs as a "
            "separate program.<br>Please cite: J. Steinmetzer, S. Kupfer, S. Gräfe, "
            "<i>Int. J. Quantum Chem.</i> <b>2021</b>, 121, e26390 (doi:10.1002/qua.26390).")
        note.setWordWrap(True)
        note.setTextFormat(Qt.RichText)
        note.setStyleSheet("color: gray")
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setFont(MONO)
        self.log.setPlaceholderText("Installation messages appear here.")
        irow = QHBoxLayout()
        irow.addWidget(self.install_btn)
        irow.addWidget(self.cancel_btn)
        irow.addStretch()
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.save)
        bb.rejected.connect(self.reject)
        lay = QVBoxLayout(self)
        lay.addWidget(QLabel("pysis program:"))
        lay.addLayout(row)
        lay.addWidget(self.state)
        lay.addLayout(irow)
        lay.addWidget(note)
        lay.addWidget(self.log, 1)
        lay.addWidget(bb)
        self.path.textChanged.connect(self.check)
        self.check()

    def check(self):
        p = self.path.text().strip()
        ok = bool(p) and os.path.isfile(p)
        self.state.setText("✓ pysis found – ready." if ok else
                           "pysis not found. Click “Install pysisyphus automatically”, or Browse.")
        self.state.setStyleSheet("color: #2a8a2a" if ok else "color: #b05000")
        self.install_btn.setText("Reinstall / repair pysisyphus" if ok and p == ts.managed_pysis(
            self.win.pysis_dir()) else f"Install pysisyphus {ts.PYSIS_VERSION} automatically")

    def browse(self):
        p, _ = QFileDialog.getOpenFileName(self, "Locate pysis", os.path.expanduser("~"))
        if p:
            self.path.setText(p)

    def install(self):
        if self.worker and self.worker.isRunning():
            return
        self.install_btn.setEnabled(False)
        self.cancel_btn.setEnabled(True)
        self.log.clear()
        self.state.setText("Installing … you can keep working in the other windows.")
        self.state.setStyleSheet("")
        base = self.win.pysis_dir()
        w = StreamWorker(lambda progress, cancel: ts.install_pysis(base, progress, cancel))
        w.setParent(self)
        w.frame.connect(lambda line: isinstance(line, str) and self.log.appendPlainText(line))
        w.done.connect(self._done)
        w.failed.connect(self._failed)
        self.worker = w
        w.start()

    def _done(self, exe):
        self.install_btn.setEnabled(True)
        self.cancel_btn.setEnabled(False)
        self.path.setText(exe)
        self.win.st.setValue("pysis/path", exe)
        self.log.appendPlainText("\nDone – pysisyphus is installed. Click OK.")

    def _failed(self, msg):
        self.install_btn.setEnabled(True)
        self.cancel_btn.setEnabled(False)
        self.state.setText("Installation stopped." if not msg else f"Installation failed: {msg}")
        self.state.setStyleSheet("color: #b05000")

    def save(self):
        self.win.st.setValue("pysis/path", self.path.text().strip())
        self.accept()

    def reject(self):
        if self.worker and self.worker.isRunning():
            if QMessageBox.question(self, "pysisyphus", "The installation is still running. Cancel it?") \
                    != QMessageBox.Yes:
                return
            self.worker.cancel.set()
            self.worker.wait(10000)
        super().reject()


# ======================================================================================
# the window
# ======================================================================================
class TSWindow(QMainWindow):
    def __init__(self, main):
        super().__init__()
        self.main = main
        self.st = main.st
        self.setWindowTitle("Transition-state search (pysisyphus) – Mol3D Bridge")
        self.resize(1320, 900)
        self.mapping = None
        self.alt_index = 0
        self.worker = None
        self.workdir = ""
        self.t_start = 0.0
        self.pending = []
        self.path_frames, self.path_k = [], 0
        self._build()

    # ------------------------------------------------------------------ layout
    def _build(self):
        tb = QToolBar("TS")
        tb.setMovable(False)
        self.addToolBar(tb)
        tb.addWidget(QLabel(" Pre-optimize with "))
        self.method = QComboBox()
        self.method.addItems(core.METHODS)
        self.method.setCurrentText(self.st.value("ts/premethod", "GFN2-xTB"))
        self.method.currentTextChanged.connect(lambda t: self.st.setValue("ts/premethod", t))
        tb.addWidget(self.method)
        tb.addSeparator()
        self.a_numbers = QAction("Atom numbers", self, checkable=True)
        self.a_numbers.setChecked(True)
        self.a_numbers.setShortcut(QKeySequence("Ctrl+L"))
        self.a_numbers.toggled.connect(self._labels)
        tb.addAction(self.a_numbers)
        a = QAction("Auto-number product", self)
        a.setToolTip("Renumber the product so every atom has the same number as in the reactant "
                     "(fewest bonds broken/formed, then least motion)")
        a.triggered.connect(self.auto_map)
        tb.addAction(a)
        self.a_alt = QAction("Other numbering", self)
        self.a_alt.setToolTip("Cycle through other, equally simple numberings")
        self.a_alt.triggered.connect(self.next_alternative)
        self.a_alt.setEnabled(False)
        tb.addAction(self.a_alt)
        self.a_relabel = QAction("Relabel by click", self, checkable=True)
        self.a_relabel.setToolTip("Click an atom (in either structure) and type its new number; "
                                  "the atom that had that number takes the old one")
        self.a_relabel.toggled.connect(lambda on: self.statusBar().showMessage(
            "Relabel: click an atom, then type its new number." if on else ""))
        tb.addAction(self.a_relabel)
        tb.addSeparator()
        a = QAction("Copy 3D", self)
        a.setShortcut(QKeySequence.Copy)
        a.setToolTip("Copy the active structure (blue frame) – Ctrl+C")
        a.triggered.connect(lambda: self.active.copy3d())
        tb.addAction(a)
        a = QAction("Paste", self)
        a.setShortcut(QKeySequence.Paste)
        a.setToolTip("Paste into the active structure (blue frame) – Ctrl+V")
        a.triggered.connect(lambda: self.active.paste())
        tb.addAction(a)
        tb.addSeparator()
        a = QAction("Swap 2 selected", self)
        a.setToolTip("Select two atoms of the same structure (click them) and swap their numbers")
        a.triggered.connect(self.swap_selected)
        tb.addAction(a)
        a = QAction("Align product onto reactant", self)
        a.triggered.connect(self.align)
        tb.addAction(a)
        tb.addSeparator()
        a = QAction("TS results folder", self)
        a.setToolTip("Open the folder where every TS search is saved")
        a.triggered.connect(lambda: self.open_path(self.results_root()))
        tb.addAction(a)
        tb.addSeparator()
        self.a_fix = QAction("Hold selected fixed", self, checkable=True)
        self.a_fix.setToolTip("Pre-optimize with the selected atoms held in place")
        tb.addAction(self.a_fix)

        self.setStyleSheet(PANEL_QSS)
        self.copied = None
        self.react = StructurePanel(self, "Reactant")
        self.prod = StructurePanel(self, "Product")
        self.active = None
        self.set_active(self.react)
        for p in (self.react, self.prod):
            p.changed.connect(self.check_numbering)
        top = QSplitter(Qt.Horizontal)
        top.addWidget(self.react)
        top.addWidget(self.prod)
        self.mapinfo = QLabel("")
        self.mapinfo.setWordWrap(True)
        topw = QWidget()
        tl = QVBoxLayout(topw)
        tl.setContentsMargins(0, 0, 0, 0)
        tl.addWidget(top, 1)
        tl.addWidget(self.mapinfo)

        self.tabs = QTabWidget()
        self.tabs.addTab(self._setup_tab(), "Calculation")
        self.tabs.addTab(self._monitor_tab(), "Monitor")
        self.results = QTabWidget()
        self.res_ts = ResultPane(self, "Transition state")
        self.res_a = ResultPane(self, "Minimum A")
        self.res_b = ResultPane(self, "Minimum B")
        self.results.addTab(self.res_ts, "Transition state")
        self.results.addTab(self.res_a, "Minimum A")
        self.results.addTab(self.res_b, "Minimum B")
        resw = QWidget()
        rl = QVBoxLayout(resw)
        self.summary = QLabel("Results appear here when pysisyphus has finished.")
        self.summary.setWordWrap(True)
        rl.addWidget(self.summary)
        rl.addWidget(self.results, 1)
        self.saved_lbl = QLabel("")
        self.saved_lbl.setOpenExternalLinks(True)
        self.saved_lbl.setTextFormat(Qt.RichText)
        brow = QHBoxLayout()
        b = QPushButton("Open this run’s folder")
        b.clicked.connect(lambda: self.open_path(getattr(self, "saved_dir", "") or self.workdir or self.results_root()))
        brow.addWidget(b)
        b = QPushButton("Open a saved run…")
        b.clicked.connect(self.load_run)
        brow.addWidget(b)
        b = QPushButton("All results")
        b.clicked.connect(lambda: self.open_path(self.results_root()))
        brow.addWidget(b)
        brow.addWidget(self.saved_lbl, 1)
        rl.addLayout(brow)
        self.tabs.addTab(resw, "Results")

        split = QSplitter(Qt.Vertical)
        split.addWidget(topw)
        split.addWidget(self.tabs)
        split.setSizes([480, 420])
        self.setCentralWidget(split)
        self.statusBar().showMessage("Build reactant and product, check the numbering, then Run.")
        self.poll = QTimer(self)
        self.poll.setInterval(1000)
        self.poll.timeout.connect(self._poll)
        self.flush = QTimer(self)
        self.flush.setInterval(200)
        self.flush.timeout.connect(self._flush)
        self.anim = QTimer(self)
        self.anim.setInterval(140)
        self.anim.timeout.connect(self._anim)

    def _setup_tab(self):
        w = QWidget()
        g = QGridLayout(w)

        def combo(items, key, default):
            c = QComboBox()
            c.addItems(items)
            c.setCurrentText(str(self.st.value(f"ts/{key}", default)))
            c.currentTextChanged.connect(self.refresh_yaml)
            return c

        def spin(lo, hi, key, default, step=1):
            s = QSpinBox()
            s.setRange(lo, hi)
            s.setSingleStep(step)
            s.setValue(int(self.st.value(f"ts/{key}", default)))
            s.valueChanged.connect(self.refresh_yaml)
            return s

        def check(key, default, text):
            c = QCheckBox(text)
            c.setChecked(str(self.st.value(f"ts/{key}", default)).lower() in ("true", "1"))
            c.toggled.connect(self.refresh_yaml)
            return c

        self.f = {}
        geo = QGroupBox("Path and TS")
        fl = QFormLayout(geo)
        self.f["coord_type"] = combo(["dlc", "redund", "cart"], "coord_type", "dlc")
        self.f["cos_type"] = combo(["gs", "neb", "szts"], "cos_type", "gs")
        self.f["max_nodes"] = spin(5, 60, "max_nodes", 18)
        self.f["climb"] = check("climb", True, "climbing image")
        self.f["opt_type"] = combo(["string", "lbfgs", "fire", "qm", "sd"], "opt_type", "string")
        self.f["max_cycles"] = spin(5, 2000, "max_cycles", 40)
        self.f["align"] = check("align", False, "align images")
        self.f["tsopt_type"] = combo(["rsirfo", "rsprfo", "trim"], "tsopt_type", "rsirfo")
        self.f["thresh"] = combo(["gau_loose", "gau", "gau_tight", "gau_vtight"], "thresh", "gau")
        self.f["hessian_recalc"] = spin(0, 100, "hessian_recalc", 1)
        self.f["do_hess"] = check("do_hess", True, "Hessian at the final TS")
        self.f["irc_type"] = combo(["eulerpc", "dvv", "lqa", "euler", "gs"], "irc_type", "eulerpc")
        self.f["rms_grad_thresh"] = QDoubleSpinBox()
        self.f["rms_grad_thresh"].setDecimals(5)
        self.f["rms_grad_thresh"].setRange(1e-5, 0.01)
        self.f["rms_grad_thresh"].setSingleStep(0.0001)
        self.f["rms_grad_thresh"].setValue(float(self.st.value("ts/rms_grad_thresh", 0.0005)))
        self.f["rms_grad_thresh"].valueChanged.connect(self.refresh_yaml)
        self.f["preopt"] = check("preopt", True, "pre-optimize ends (preopt)")
        self.f["endopt"] = check("endopt", True, "optimize IRC ends (endopt)")
        for lab, key in (("Coordinates", "coord_type"), ("Chain of states", "cos_type"),
                         ("Max. nodes", "max_nodes"), ("", "climb"), ("COS optimizer", "opt_type"),
                         ("Max. cycles", "max_cycles"), ("", "align"), ("TS optimizer", "tsopt_type"),
                         ("Convergence", "thresh"), ("Hessian recalc", "hessian_recalc"), ("", "do_hess"),
                         ("IRC", "irc_type"), ("IRC rms(grad)", "rms_grad_thresh"), ("", "preopt"),
                         ("", "endopt")):
            fl.addRow(lab, self.f[key])
        calc = QGroupBox("Calculator (xtb)")
        cl = QFormLayout(calc)
        self.f["gfn"] = combo(["2", "1", "0", "ff"], "gfn", "2")
        self.f["pal"] = spin(1, 256, "pal", min(6, os.cpu_count() or 1))
        self.f["etemp"] = spin(0, 20000, "etemp", 1500, 100)
        self.f["charge"] = QSpinBox()
        self.f["charge"].setRange(-20, 20)
        self.f["charge"].valueChanged.connect(self.refresh_yaml)
        self.f["mult"] = QSpinBox()
        self.f["mult"].setRange(1, 15)
        self.f["mult"].valueChanged.connect(self.refresh_yaml)
        self.f["solvent"] = combo(list(core.ALPB_SOLVENTS), "solvent", "none")
        for lab, key in (("GFN", "gfn"), ("Cores (pal)", "pal"), ("Electronic temp. (K)", "etemp"),
                         ("Charge", "charge"), ("Multiplicity", "mult"), ("Solvent (ALPB)", "solvent")):
            cl.addRow(lab, self.f[key])
        self.pysis_lbl = QLabel("")
        b = QPushButton("pysisyphus program…")
        b.clicked.connect(lambda: (PysisDialog(self).exec(), self._pysis_status()))
        cl.addRow(b, self.pysis_lbl)
        self._pysis_status()

        self.custom = QCheckBox("Edit the YAML by hand (form changes no longer overwrite it)")
        self.yaml = QPlainTextEdit()
        self.yaml.setFont(MONO)
        self.run_btn = QPushButton("▶ Run TS search")
        self.run_btn.setStyleSheet("font-weight: 600; padding: 6px 14px")
        self.run_btn.clicked.connect(self.run)
        self.stop_btn = QPushButton("■ Stop")
        self.stop_btn.setEnabled(False)
        self.stop_btn.clicked.connect(self.stop)
        self.folder_lbl = QLabel("")
        self.folder_lbl.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.folder_lbl.setWordWrap(True)
        runrow = QHBoxLayout()
        for x in (self.run_btn, self.stop_btn):
            runrow.addWidget(x)
        runrow.addWidget(self.folder_lbl, 1)
        # where results are stored
        self.run_name = QLineEdit()
        self.run_name.setPlaceholderText("optional, e.g. Claisen chair – used in the folder name")
        self.autosave = QCheckBox("Save every run automatically")
        self.autosave.setChecked(self.st.value("ts/autosave", True, type=bool))
        self.autosave.toggled.connect(lambda v: self.st.setValue("ts/autosave", v))
        self.keepall = QCheckBox("including all pysisyphus files")
        self.keepall.setChecked(self.st.value("ts/keepall", True, type=bool))
        self.keepall.toggled.connect(lambda v: self.st.setValue("ts/keepall", v))
        self.root_lbl = QLabel("")
        self.root_lbl.setTextInteractionFlags(Qt.TextSelectableByMouse)
        b_root = QPushButton("Change…")
        b_root.clicked.connect(self.choose_root)
        b_open = QPushButton("Open")
        b_open.clicked.connect(lambda: self.open_path(self.results_root()))
        save = QGroupBox("Results")
        sl = QGridLayout(save)
        sl.addWidget(QLabel("Run name"), 0, 0)
        sl.addWidget(self.run_name, 0, 1, 1, 3)
        sl.addWidget(QLabel("Results folder"), 1, 0)
        sl.addWidget(self.root_lbl, 1, 1)
        sl.addWidget(b_root, 1, 2)
        sl.addWidget(b_open, 1, 3)
        sl.addWidget(self.autosave, 2, 1)
        sl.addWidget(self.keepall, 2, 2, 1, 2)
        sl.setColumnStretch(1, 1)
        self._root_label()
        right = QVBoxLayout()
        right.addWidget(QLabel("pysisyphus input (ts.yaml):"))
        right.addWidget(self.yaml, 1)
        right.addWidget(self.custom)
        right.addWidget(save)
        right.addLayout(runrow)
        g.addWidget(geo, 0, 0, 2, 1)
        g.addWidget(calc, 0, 1)
        g.addLayout(right, 0, 2, 2, 1)
        g.setColumnStretch(2, 1)
        self.refresh_yaml()
        return w

    def _monitor_tab(self):
        w = QWidget()
        lay = QHBoxLayout(w)
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setFont(MONO)
        self.log.setMaximumBlockCount(20000)
        self.log.setPlaceholderText("pysisyphus output appears here while the calculation runs.")
        right = QVBoxLayout()
        self.stage_lbl = QLabel("")
        self.stage_lbl.setTextFormat(Qt.RichText)
        self.path_view = MolView()
        self.path_view.setMinimumSize(320, 220)
        self.path_view.empty_text = "The reaction path appears here while the\nchain of states is optimized."
        self.profile = ProfilePlot()
        right.addWidget(self.stage_lbl)
        right.addWidget(QLabel("Current reaction path (latest chain-of-states cycle):"))
        right.addWidget(self.path_view, 1)
        right.addWidget(self.profile)
        sp = QSplitter(Qt.Horizontal)
        sp.addWidget(self.log)
        rw = QWidget()
        rw.setLayout(right)
        sp.addWidget(rw)
        sp.setSizes([700, 500])
        lay.addWidget(sp)
        self._stage("start")
        return w

    # ------------------------------------------------------------------ helpers
    def set_active(self, panel):
        if panel is self.active:
            return
        self.active = panel
        for p in (self.react, self.prod):
            p.set_active(p is panel)
        self.statusBar().showMessage(f"{panel.role} is active – Copy, Paste, Swap and Relabel act on it.")

    def error(self, msg):
        self.statusBar().showMessage(msg.splitlines()[0])
        QMessageBox.warning(self, "TS search", msg)

    def lastdir(self):
        return self.st.value("lastdir", os.path.expanduser("~"))

    def _labels(self, on):
        for v in (self.react.view, self.prod.view, self.res_ts.view, self.res_a.view, self.res_b.view,
                  self.path_view):
            v.labels = "index" if on else "hetero"
            v.update()

    def results_root(self):
        default = os.path.join(QStandardPaths.writableLocation(QStandardPaths.DocumentsLocation)
                               or os.path.expanduser("~"), "Mol3D Bridge", "TS results")
        return self.st.value("ts/results_root", default) or default

    def _root_label(self):
        self.root_lbl.setText(self.results_root())

    def choose_root(self):
        d = QFileDialog.getExistingDirectory(self, "Folder for TS search results", self.results_root())
        if d:
            self.st.setValue("ts/results_root", d)
            self._root_label()

    def open_path(self, path):
        os.makedirs(path, exist_ok=True)
        QDesktopServices.openUrl(QUrl.fromLocalFile(path))

    def store(self, status):
        """Copy the run into the permanent results folder (automatic after every run)."""
        if not self.autosave.isChecked() or not self.workdir or not os.path.isdir(self.workdir):
            return None
        try:
            dest = ts.save_results(self.workdir, self.results_root(), self.react_snapshot, self.prod_snapshot,
                                   self.run_name.text(), status, self.keepall.isChecked(), self.t_wall)
        except Exception as ex:  # noqa: BLE001
            self.error(f"The results could not be saved to {self.results_root()}:\n{ex}\n\n"
                       f"They are still in {self.workdir}.")
            return None
        if self.keepall.isChecked():
            import shutil
            shutil.rmtree(self.workdir, ignore_errors=True)        # the full copy is in the results folder
            self.workdir = ts.find_run_folder(dest) or dest
        self.saved_dir = dest
        self.folder_lbl.setText(f"Saved in {dest}")
        self.saved_lbl.setText(f"Saved in <a href='{QUrl.fromLocalFile(dest).toString()}'>{dest}</a>")
        return dest

    def pysis_dir(self):
        return os.path.join(user_data_dir(), "pysisyphus")

    def pysis_path(self):
        return ts.find_pysis(self.st.value("pysis/path", ""), self.pysis_dir())

    def _pysis_status(self):
        p = self.pysis_path()
        self.pysis_lbl.setText("✓ found" if p else "not installed – click to install")
        self.pysis_lbl.setStyleSheet("color: #2a8a2a" if p else "color: #b05000")
        self.pysis_lbl.setToolTip(p or "pip install pysisyphus")

    def settings(self) -> ts.TSSettings:
        f = self.f
        gfn = f["gfn"].currentText()
        return ts.TSSettings(
            coord_type=f["coord_type"].currentText(), cos_type=f["cos_type"].currentText(),
            max_nodes=f["max_nodes"].value(), climb=f["climb"].isChecked(),
            opt_type=f["opt_type"].currentText(), align=f["align"].isChecked(),
            max_cycles=f["max_cycles"].value(), tsopt_type=f["tsopt_type"].currentText(),
            do_hess=f["do_hess"].isChecked(), thresh=f["thresh"].currentText(),
            hessian_recalc=f["hessian_recalc"].value(), irc_type=f["irc_type"].currentText(),
            rms_grad_thresh=f["rms_grad_thresh"].value(), preopt=f["preopt"].isChecked(),
            endopt=f["endopt"].isChecked(), gfn=gfn if gfn == "ff" else int(gfn),
            pal=f["pal"].value(), etemp=f["etemp"].value(), charge=f["charge"].value(),
            mult=f["mult"].value(), solvent=f["solvent"].currentText())

    def refresh_yaml(self, *_):
        if not hasattr(self, "yaml") or self.custom.isChecked():
            return
        s = self.settings()
        for k, v in s.__dict__.items():
            if k not in ("charge", "mult"):
                self.st.setValue(f"ts/{k}", v)
        self.yaml.setPlainText(ts.yaml_input(s))

    # ------------------------------------------------------------------ numbering
    def check_numbering(self):
        r, p = self.react.mol, self.prod.mol
        if r is not None:
            ch, mu = core.default_charge_mult(r)
            self.f["charge"].blockSignals(True)
            self.f["mult"].blockSignals(True)
            self.f["charge"].setValue(ch)
            self.f["mult"].setValue(mu)
            self.f["charge"].blockSignals(False)
            self.f["mult"].blockSignals(False)
            self.refresh_yaml()
        if r is None or p is None:
            self.mapinfo.setText("")
            for v in (self.react.view, self.prod.view):
                v.extra_bonds = []
            return
        if r.GetNumAtoms() != p.GetNumAtoms() or \
                [a.GetSymbol() for a in r.GetAtoms()] != [a.GetSymbol() for a in p.GetAtoms()]:
            same = ts._formula(r) == ts._formula(p)
            self.mapinfo.setText(
                "⚠ Atom numbering differs (atom i is not the same element on both sides) – click "
                "<b>Auto-number product</b>." if same else
                "⚠ Reactant and product have different formulas – both sides need the same atoms.")
            self.mapinfo.setStyleSheet("color: #b05000; font-weight: 600")
            for v in (self.react.view, self.prod.view):
                v.extra_bonds = []
                v.update()
            return
        rb, pb = ts._bonds(r), ts._bonds(p)
        m = ts.Mapping(list(range(r.GetNumAtoms())), 0, 0.0, sorted(rb - pb), sorted(pb - rb))
        nchg = len(rb ^ pb)
        # show the changes in the structures: on the reactant the bonds that will form, on the
        # product the bonds that were broken (dashed orange)
        self.react.view.extra_bonds = sorted(pb - rb)
        self.prod.view.extra_bonds = sorted(rb - pb)
        self.react.view.update()
        self.prod.view.update()
        warn = nchg > 6
        self.mapinfo.setText(("⚠ " if warn else "") + f"With this numbering the reaction "
                             f"{ts.describe_changes(m, r)} ({nchg} bond changes; dashed orange: bonds "
                             "formed – left, broken – right)."
                             + (" That is a lot – the numbering may be off. Swap atoms (select two, "
                                "“Swap 2 selected”, or “Relabel by click”) or use Auto-number." if warn else ""))
        self.mapinfo.setStyleSheet("color: #b05000; font-weight: 600" if warn else "")

    def auto_map(self):
        r, p = self.react.mol, self.prod.mol
        if r is None or p is None:
            return self.error("Load both a reactant and a product first.")
        try:
            m = ts.map_atoms(r, p)
        except ValueError as ex:
            return self.error(str(ex))
        self.mapping, self.alt_index = m, 0
        self.prod.renumber(m.perm, "renumbered to match the reactant")
        ts.align_onto(self.prod.mol, self.react.mol)
        self.prod.view.recenter()
        self.prod.view.R = self.react.view.R.copy()
        self.prod.view.update()
        self.a_alt.setEnabled(bool(m.alternatives))
        extra = f"  {len(m.alternatives)} other equally simple numbering(s): 'Other numbering'." \
            if m.alternatives else ""
        self.statusBar().showMessage(f"Product renumbered: {ts.describe_changes(m, r)}.{extra}")
        self._orig_prod_for_alts = None

    def next_alternative(self):
        m = self.mapping
        if not m or not m.alternatives:
            return
        # alternatives are permutations of the product as it was when auto-mapping ran
        cur_to_orig = m.perm if self.alt_index == 0 else m.alternatives[self.alt_index - 1]
        self.alt_index = (self.alt_index + 1) % (len(m.alternatives) + 1)
        target = m.perm if self.alt_index == 0 else m.alternatives[self.alt_index - 1]
        inv = {o: i for i, o in enumerate(cur_to_orig)}            # original index -> current index
        self.prod.renumber([inv[o] for o in target], f"numbering option {self.alt_index + 1}")
        ts.align_onto(self.prod.mol, self.react.mol)
        self.prod.view.update()
        self.check_numbering()

    def relabel(self, panel, i, sym):
        n = panel.mol.GetNumAtoms()
        new, ok = QInputDialog.getInt(self, "Relabel atom", f"New number for {sym}{i} in the "
                                      f"{panel.role.lower()} (0–{n - 1}):", i, 0, n - 1)
        panel.view.sel = []
        if ok and new != i:
            panel.renumber([new if k == i else i if k == new else k for k in range(n)],
                           f"{sym}{i} is now number {new}")
            self.mapping = None
            self.a_alt.setEnabled(False)
        panel.view.update()

    def swap_selected(self):
        order = [self.active] + [p for p in (self.react, self.prod) if p is not self.active]
        for panel in order:
            if panel.mol is not None and len(panel.view.sel) == 2:
                i, j = panel.view.sel
                n = panel.mol.GetNumAtoms()
                panel.view.sel = []
                panel.renumber([j if k == i else i if k == j else k for k in range(n)], f"swapped {i} ↔ {j}")
                self.mapping = None
                self.a_alt.setEnabled(False)
                return
        self.error("Click two atoms of the same structure first (they get orange rings), then Swap.")

    def align(self):
        r, p = self.react.mol, self.prod.mol
        if r is None or p is None or r.GetNumAtoms() != p.GetNumAtoms():
            return self.error("Both structures with the same numbering are needed.")
        rmsd = ts.align_onto(p, r)
        self.prod.view.recenter()
        self.prod.view.R = self.react.view.R.copy()
        self.prod.view.update()
        self.statusBar().showMessage(f"Product superposed on the reactant: RMSD {rmsd:.2f} Å.")

    # ------------------------------------------------------------------ running
    def _stage(self, key, detail=""):
        keys = [k for k, _ in ts.STAGES]
        idx = keys.index(key) if key in keys else -1
        rows = []
        if key == "cos" and detail:
            detail = "cycle " + detail.replace("cycle_", "").replace(".trj", "").lstrip("0") or "0"
        for k, (kk, lab) in enumerate(ts.STAGES):
            mark = "✓" if k < idx or key == "done" else ("▶" if k == idx else "○")
            style = "font-weight:600" if k == idx and key != "done" else ("color:gray" if k > idx else "")
            rows.append(f"<span style='{style}'>{mark} {lab}{(' – ' + detail) if k == idx and detail else ''}</span>")
        t = f"  ·  {time.monotonic() - self.t_start:.0f} s" if self.t_start and key != "start" else ""
        self.stage_lbl.setText("<br>".join(rows) + f"<br><span style='color:gray'>{t}</span>")

    def run(self):
        if self.worker and self.worker.isRunning():
            return
        r, p = self.react.mol, self.prod.mol
        if r is None or p is None:
            return self.error("Load both a reactant and a product first.")
        if [a.GetSymbol() for a in r.GetAtoms()] != [a.GetSymbol() for a in p.GetAtoms()]:
            bad = [i for i, (a, b) in enumerate(zip(r.GetAtoms(), p.GetAtoms())) if a.GetSymbol() != b.GetSymbol()]
            return self.error("pysisyphus needs the same element at every atom number. These numbers "
                              f"differ between reactant and product: {', '.join(map(str, bad[:12]))}"
                              f"{' …' if len(bad) > 12 else ''}.\n\nSwap atoms (select two, “Swap 2 "
                              "selected”, or “Relabel by click”) or use Auto-number product.")
        nchg = len(ts._bonds(r) ^ ts._bonds(p))
        if nchg > 6:      # only a reminder – the user decides the numbering
            self.statusBar().showMessage(f"Note: {nchg} bond changes with the current numbering.")
        pysis = self.pysis_path()
        if not pysis:
            PysisDialog(self).exec()
            self._pysis_status()
            pysis = self.pysis_path()
            if not pysis:
                return
        xtb = self.main.xtb_settings().path
        if self.settings().calc_type == "xtb" and not xtb:
            return self.error("xtb was not found – set it up in the main window (Calculate ▸ xtb settings…).")
        stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        self.workdir = os.path.join(user_data_dir(), "ts_runs", f"{stamp}_{core.formula_and_mass(r)[0]}")
        ts.write_inputs(self.workdir, r, p, self.settings(), self.yaml.toPlainText(), xtb)
        self.react_snapshot, self.prod_snapshot = Chem.Mol(r), Chem.Mol(p)
        self.t_wall = time.time()
        self.saved_dir = ""
        self.saved_lbl.setText("")
        self.folder_lbl.setText(self.workdir)
        self.log.clear()
        self.t_start = time.monotonic()
        self._stage("start")
        self.tabs.setCurrentIndex(1)
        self.run_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)
        w = StreamWorker(lambda progress, cancel: ts.run_pysis(pysis, self.workdir, xtb, progress, cancel))
        w.setParent(self)
        w.frame.connect(lambda line: self.pending.append(line))
        w.done.connect(self._finished)
        w.failed.connect(self._failed)
        self.worker = w
        w.start()
        self.poll.start()
        self.flush.start()
        self.statusBar().showMessage(f"pysisyphus running in {self.workdir}")

    def stop(self):
        if self.worker and self.worker.isRunning():
            self.worker.cancel.set()
            self.statusBar().showMessage("Stopping pysisyphus…")

    def _flush(self):
        if self.pending:
            lines, self.pending = self.pending, []
            self.log.appendPlainText("\n".join(lines))

    def _poll(self):
        if not self.workdir:
            return
        key, detail = ts.current_stage(self.workdir)
        self._stage(key, detail)
        lp = ts.latest_path(self.workdir)
        if lp and lp[0] != getattr(self, "_last_cycle", None):
            self._last_cycle = lp[0]
            name, frames, energies = lp
            self.path_frames = frames
            tmpl = self.react_snapshot
            if self.path_view.mol is None or self.path_view.mol.GetNumAtoms() != tmpl.GetNumAtoms():
                self.path_view.set_molecule(ts.mol_with_positions(tmpl, frames[0]), 0)
                self.path_view.labels = "index" if self.a_numbers.isChecked() else "hetero"
            e = [x for x in energies if x is not None]
            if e:
                lo = min(e)
                kj = [None if x is None else (x - lo) * 2625.4996 for x in energies]
                hi = max(v for v in kj if v is not None)
                self.profile.set_data(kj, f"{name}: {len(frames)} images, highest +{hi:.0f} kJ/mol")
            if not self.anim.isActive():
                self.anim.start()

    def _anim(self):
        n = len(self.path_frames)
        if not n or self.path_view.mol is None:
            return
        k = self.path_k % n
        if self.path_view.mol.GetNumAtoms() == len(self.path_frames[k]):
            core.apply_positions(self.path_view.mol, 0, self.path_frames[k])
            self.path_view.update()
        self.profile.set_current(k)
        self.path_k += 1

    def _end(self):
        self.poll.stop()
        self._flush()
        self.flush.stop()
        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)

    def _failed(self, msg):
        self._end()
        self.anim.stop()
        dest = self.store("stopped" if not msg else "failed")
        where = f" Partial results saved in {dest}." if dest else ""
        if not msg:
            self.statusBar().showMessage("Stopped." + where)
            self._stage("start", "stopped")
        else:
            self.error(f"pysisyphus could not be run:\n{msg}" + (f"\n\n{where.strip()}" if where else ""))

    def _finished(self, rc):
        self._end()
        self._poll()
        self.anim.stop()
        res = ts.collect_results(self.workdir, self.react_snapshot, self.prod_snapshot)
        status = "finished" if res.ts is not None else "no-TS"
        dest = self.store(status)
        self.show_results(self.workdir, self.react_snapshot, self.prod_snapshot,
                          f"pysisyphus exited with code {rc}." if rc else "")
        if dest:
            self.statusBar().showMessage(f"TS search {'finished' if res.ts is not None else 'ended'} – "
                                         f"saved in {dest}")

    def load_run(self):
        d = QFileDialog.getExistingDirectory(self, "Saved TS search", self.results_root())
        if not d:
            return
        saved = d
        d = ts.find_run_folder(d) or d
        try:
            r = core.mol_from_xyz(open(os.path.join(d, "reactant.xyz")).read())
            p = core.mol_from_xyz(open(os.path.join(d, "product.xyz")).read())
        except OSError:
            return self.error("That folder has no reactant.xyz / product.xyz from Mol3D Bridge.")
        self.workdir = d
        self.saved_dir = saved
        self.folder_lbl.setText(saved)
        self.saved_lbl.setText(f"Opened <a href='{QUrl.fromLocalFile(saved).toString()}'>{saved}</a>")
        self.show_results(d, r, p)

    def show_results(self, workdir, react, prod, note=""):
        try:
            res = ts.collect_results(workdir, react, prod)
        except Exception as ex:  # noqa: BLE001
            return self.error(f"Could not read the results: {ex}")
        side_name = {"reactant": "reactant side", "product": "product side", "other": "other"}
        names = {}
        for key, pane, default in (("minimum1", self.res_a, "Minimum A"), ("minimum2", self.res_b, "Minimum B")):
            side = res.sides.get(key)
            names[key] = f"{default} ({side_name[side]})" if side else default
            pane.name = names[key]
        rel = res.relative_kj

        def text(key, title):
            parts = [f"<b>{title}</b>"]
            if key in res.energies:
                parts.append(f"E = {res.energies[key]:.6f} Eh")
            if key in rel:
                parts.append(f"{rel[key]:+.1f} kJ/mol vs. the lower minimum")
            if key == "ts" and res.imag_freq is not None:
                parts.append(f"imaginary frequency {res.imag_freq:.0f} cm⁻¹")
            if key in res.files:
                parts.append(f"<span style='color:gray'>{res.files[key]}</span>")
            return " · ".join(parts)

        irc = [f[0] for f in res.irc_frames]
        self.res_ts.set_result(res.ts, text("ts", "Transition state"), irc if res.ts is not None else None)
        rb, pb = ts._bonds(react), ts._bonds(prod)
        self.res_ts.view.extra_bonds = sorted(rb ^ pb)          # forming / breaking bonds, dashed
        self.res_a.set_result(res.minimum1, text("minimum1", names["minimum1"]))
        self.res_b.set_result(res.minimum2, text("minimum2", names["minimum2"]))
        self.results.setTabText(1, names["minimum1"])
        self.results.setTabText(2, names["minimum2"])
        for pane in (self.res_ts, self.res_a, self.res_b):
            if pane.mol is not None and self.react.mol is not None and \
                    pane.mol.GetNumAtoms() == self.react.mol.GetNumAtoms():
                ts.align_onto(pane.mol, self.react.mol)
                pane.view.R = self.react.view.R.copy()
                pane.view.recenter()
                pane.view.update()
        s = []
        if "ts" in rel:
            lo = {k: rel[k] for k in ("minimum1", "minimum2") if k in rel}
            for k, v in lo.items():
                s.append(f"barrier from {names[k]}: <b>{rel['ts'] - v:.1f} kJ/mol</b>")
        if res.imag_freq is not None:
            s.append(f"imaginary frequency {res.imag_freq:.0f} cm⁻¹")
        s += res.notes
        if note:
            s.append(note)
        self.summary.setText("  ·  ".join(s) or "Finished.")
        self.tabs.setCurrentIndex(2)
        self.results.setCurrentIndex(0)
        self._stage("done" if res.ts is not None else "start")
        self.statusBar().showMessage("TS search finished." if res.ts is not None else "No TS found – see Monitor.")

    def closeEvent(self, e):
        if self.worker and self.worker.isRunning():
            if QMessageBox.question(self, "TS search", "pysisyphus is still running. Stop it?") != QMessageBox.Yes:
                e.ignore()
                return
            self.worker.cancel.set()
            self.worker.wait(5000)
        super().closeEvent(e)
