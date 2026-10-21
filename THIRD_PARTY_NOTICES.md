# Third-party software in Mol3D Bridge

Mol3D Bridge itself is © 2026 Takuma Miyamura and released under the MIT licence (`LICENSE`).
The installers also contain the components below; each keeps its own licence. Full licence
texts are in the `licenses/` folder (and, in installers, `third_party_licenses/`).

| Component | Purpose | Licence | Source |
|---|---|---|---|
| RDKit | cheminformatics, 3D embedding, MMFF94/UFF, ChemDraw file I/O | BSD-3-Clause | https://github.com/rdkit/rdkit |
| Qt 6 / Qt for Python (PySide6, Shiboken6) | user interface | LGPL-3.0 | see "Qt" below |
| NumPy | numerics | BSD-3-Clause (bundled libraries: see `third_party_licenses/`) | https://github.com/numpy/numpy |
| Python | runtime | PSF License | https://www.python.org |
| PyObjC (macOS only) | clipboard access | MIT | https://github.com/ronaldoussoren/pyobjc |
| xtb 6.7.1 (Windows/Linux installers) | GFN2-xTB / GFN-FF | LGPL-3.0 | see "xtb" below |
| PyInstaller bootloader | starts the program | GPL-2.0 with bootloader exception | https://github.com/pyinstaller/pyinstaller |
| Inno Setup (Windows) | setup program | Inno Setup licence | https://jrsoftware.org |

## Qt / PySide6 (LGPL-3.0)

Qt and PySide6 are used unmodified, as separate shared libraries (the `PySide6` folder inside
the installation). You may replace them with other compatible builds; the source of Mol3D
Bridge is public, so the program can also be rebuilt against a modified Qt. The complete
corresponding source code of Qt and PySide6 is available from The Qt Company at
https://download.qt.io/official_releases/QtForPython/ (PySide6 / Shiboken6) and
https://download.qt.io/official_releases/qt/ (Qt). The exact versions included in an installer
are listed in its `third_party_licenses/VERSIONS.txt`. On request, the author will provide the
corresponding source for the bundled versions.

## xtb (LGPL-3.0)

xtb 6.7.1 is included unmodified as a separate program in the `xtb` folder, with its licence
files. Source: https://github.com/grimme-lab/xtb/releases/tag/v6.7.1 . When publishing results
obtained with GFN2-xTB please cite C. Bannwarth, S. Ehlert, S. Grimme,
*J. Chem. Theory Comput.* 2019, 15, 1652.

## pysisyphus (not included)

pysisyphus © Johannes Steinmetzer and contributors, licensed under the **GNU General Public
License v3** (checked in the licence file of pysisyphus 1.0.0 on PyPI;
https://github.com/eljost/pysisyphus). It is **not distributed** with Mol3D Bridge: the user
installs it (`pip install pysisyphus`, or the one-click installation, which downloads it from PyPI
together with a portable CPython 3.12 from python-build-standalone, PSF License,
https://github.com/astral-sh/python-build-standalone). Mol3D Bridge runs it as a separate process
– it writes `ts.yaml` and the xyz inputs, starts pysisyphus through a 20-line launcher script
(`run_pysis.py`, which only skips a Windows symbolic link) and reads the output files. The MIT
licence of Mol3D Bridge is compatible with the GPL, so this use is fine; if you redistribute
pysisyphus yourself, the GPL v3 terms apply to it.

Citation requested by the pysisyphus authors:
J. Steinmetzer, S. Kupfer, S. Gräfe, "pysisyphus: Exploring potential energy surfaces in ground
and excited states", *Int. J. Quantum Chem.* **2021**, 121, e26390. doi:10.1002/qua.26390

## Trademarks

ChemDraw and Chem3D are trademarks of Revvity, Inc. Mol3D Bridge is an independent project and
is not affiliated with, sponsored or endorsed by Revvity. ReShiki is an independent open-source
project; Mol3D Bridge contains no code from it. The names are used only to describe compatibility.
