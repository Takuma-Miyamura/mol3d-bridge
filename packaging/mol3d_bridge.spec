# PyInstaller spec for Mol3D Bridge  (built by packaging/build.py)
import os, sys
from PyInstaller.utils.hooks import collect_data_files, collect_dynamic_libs

ROOT = os.path.abspath(os.path.join(SPECPATH, ".."))
SRC = os.path.join(ROOT, "src")
sys.path.insert(0, SRC)
from mol3d_bridge import __version__  # noqa: E402

NAME = "Mol3D Bridge"
XTB_DIR = os.environ.get("M3B_XTB_DIR", "")

datas = [(os.path.join(SRC, "mol3d_bridge", "resources"), "mol3d_bridge/resources"),
         (os.path.join(ROOT, "LICENSE"), "."),
         (os.path.join(ROOT, "THIRD_PARTY_NOTICES.md"), "."),
         (os.path.join(ROOT, "licenses"), "licenses")]
_tpl = os.path.join(ROOT, "build", "third_party_licenses")
if os.path.isdir(_tpl):
    datas.append((_tpl, "third_party_licenses"))
datas += collect_data_files("rdkit", excludes=["**/Contrib/**", "**/test_data/**", "**/Docs/**"])
if XTB_DIR and os.path.isdir(XTB_DIR):
    datas.append((XTB_DIR, "xtb"))

binaries = collect_dynamic_libs("rdkit")

hidden = ["rdkit.Chem.rdDetermineBonds", "rdkit.Chem.rdMolTransforms", "rdkit.Chem.rdChemDraw",
          "rdkit.Chem.rdDepictor", "rdkit.Chem.rdMolDescriptors", "rdkit.Chem.Descriptors",
          "rdkit.Chem.rdForceFieldHelpers", "rdkit.ForceField.rdForceField", "rdkit.Chem.rdDistGeom"]
if sys.platform == "darwin":
    hidden += ["AppKit", "Foundation", "objc"]

excludes = ["tkinter", "matplotlib", "pandas", "IPython", "PIL", "scipy", "pytest",
            "PySide6.QtWebEngineCore", "PySide6.QtWebEngineWidgets", "PySide6.QtQml",
            "PySide6.QtQuick", "PySide6.Qt3DCore", "PySide6.QtMultimedia", "PySide6.QtPdf",
            "PySide6.QtCharts", "PySide6.QtDataVisualization", "PySide6.QtSql", "PySide6.QtTest"]

a = Analysis([os.path.join(ROOT, "packaging", "launcher.py")], pathex=[SRC], binaries=binaries,
             datas=datas, hiddenimports=hidden, excludes=excludes, noarchive=False)
# Qt parts the app never uses (QML/Quick, PDF, virtual keyboard) - pulled in only by plugins
_drop = ("Qt6Quick", "Qt6Qml", "Qt6Pdf", "Qt6VirtualKeyboard", "qpdf", "qtvirtualkeyboard",
         "Qt6OpenGLWidgets", "Qt6QuickWidgets")
a.binaries = [b for b in a.binaries if not any(k in os.path.basename(b[0]) for k in _drop)]
a.datas = [d for d in a.datas if "/translations/qtwebengine" not in d[0].replace("\\", "/")]
pyz = PYZ(a.pure)
icon = os.path.join(SRC, "mol3d_bridge", "resources",
                    "icon.icns" if sys.platform == "darwin" else "icon.ico")
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name=NAME, console=False,
          icon=icon, upx=False,
          version=os.path.join(ROOT, "build", "version_info.txt") if sys.platform.startswith("win") else None)
coll = COLLECT(exe, a.binaries, a.datas, name=NAME, upx=False)

if sys.platform == "darwin":
    exts = ["cdx", "cdxml", "mol", "sdf", "xyz", "pdb"]
    app = BUNDLE(coll, name=f"{NAME}.app", icon=icon, bundle_identifier="io.github.mol3dbridge",
                 info_plist={
                     "CFBundleShortVersionString": __version__,
                     "CFBundleVersion": __version__,
                     "NSHighResolutionCapable": True,
                     "LSMinimumSystemVersion": "11.0",
                     "NSHumanReadableCopyright": "Mol3D Bridge contributors, MIT licence",
                     "CFBundleDocumentTypes": [{
                         "CFBundleTypeName": "Chemical structure",
                         "CFBundleTypeRole": "Viewer",
                         "LSHandlerRank": "Alternate",
                         "CFBundleTypeExtensions": exts,
                     }],
                 })
