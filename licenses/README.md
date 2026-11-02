# Licence texts

Full texts of the licences of the software that the Mol3D Bridge installers include.
See `../THIRD_PARTY_NOTICES.md` for which licence applies to which component.

| File | Licence | Used by |
|---|---|---|
| `LGPL-3.0.txt` | GNU Lesser General Public License v3 | Qt 6 / PySide6, xtb |
| `GPL-3.0.txt` | GNU General Public License v3 (the LGPL builds on it) | Qt 6 / PySide6, xtb |
| `RDKit-BSD-3-Clause.txt` | BSD 3-Clause | RDKit |
| `NumPy-BSD-3-Clause.txt` | BSD 3-Clause | NumPy |
| `Python-PSF.txt` | Python Software Foundation License | Python runtime |
| `PyObjC-MIT.txt` | MIT | PyObjC (macOS installer only) |
| `PyInstaller-bootloader.txt` | GPL-2.0 with bootloader exception | PyInstaller start-up program |
| `InnoSetup.txt` | Inno Setup licence | Windows setup program |

Installers built with `packaging/build.py` additionally contain, in
`third_party_licenses/`, the licence files shipped inside the exact versions of the Python
packages that were bundled (including the libraries NumPy and RDKit wheels carry, such as
OpenBLAS, Boost, FreeType and the GCC runtime).
