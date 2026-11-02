#!/usr/bin/env python3
"""
Build the Mol3D Bridge desktop application and its installer for the current platform.

    python packaging/build.py              # app + installer
    python packaging/build.py --with-xtb   # also bundle the official xtb 6.7.1 (Windows/Linux)

Outputs in dist/:
  Windows : "Mol3D Bridge/" folder + Mol3DBridge-<ver>-Windows-x64-Setup.exe   (needs Inno Setup 6)
  macOS   : "Mol3D Bridge.app"   + Mol3DBridge-<ver>-macOS-<arch>.dmg
  Linux   : "Mol3D Bridge/" folder + Mol3DBridge-<ver>-Linux-x86_64.tar.gz (with install.sh)
"""
import argparse
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC = os.path.join(ROOT, "src")
BUILD = os.path.join(ROOT, "build")
DIST = os.path.join(ROOT, "dist")
NAME = "Mol3D Bridge"
sys.path.insert(0, SRC)
from mol3d_bridge import __version__ as VERSION  # noqa: E402

WIN = sys.platform.startswith("win")
MAC = sys.platform == "darwin"


def run(cmd, **kw):
    print("+", " ".join(f'"{c}"' if " " in str(c) else str(c) for c in cmd), flush=True)
    subprocess.run(cmd, check=True, **kw)


def version_info():
    v = tuple(int(x) for x in (VERSION.split(".") + ["0", "0", "0"])[:4])
    text = f"""VSVersionInfo(
  ffi=FixedFileInfo(filevers={v}, prodvers={v}),
  kids=[StringFileInfo([StringTable('040904B0', [
      StringStruct('CompanyName', 'Takuma Miyamura'),
      StringStruct('FileDescription', '{NAME}'),
      StringStruct('FileVersion', '{VERSION}'),
      StringStruct('ProductName', '{NAME}'),
      StringStruct('ProductVersion', '{VERSION}'),
      StringStruct('LegalCopyright', 'Copyright (c) 2026 Takuma Miyamura, MIT licence')])]),
    VarFileInfo([VarStruct('Translation', [1033, 1200])])])
"""
    os.makedirs(BUILD, exist_ok=True)
    with open(os.path.join(BUILD, "version_info.txt"), "w") as fh:
        fh.write(text)


def bundle_xtb():
    from mol3d_bridge import core
    if MAC:
        print("xtb: no official macOS binary - not bundled (users can install it via conda).")
        return ""
    tmp = tempfile.mkdtemp(prefix="xtb_dl_")
    exe = core.download_xtb(tmp, lambda m: print(" ", m, end="\r"))
    print()
    root = os.path.dirname(os.path.dirname(exe))      # .../xtb-6.7.1
    dest = os.path.join(BUILD, "xtb_bundle")
    shutil.rmtree(dest, ignore_errors=True)
    os.makedirs(os.path.join(dest, "bin"))
    for f in os.listdir(os.path.join(root, "bin")):     # xtb + its runtime DLLs only
        if f in ("xtb", "xtb.exe") or f.lower().endswith(".dll"):
            shutil.copy2(os.path.join(root, "bin", f), os.path.join(dest, "bin", f))
    if os.path.isdir(os.path.join(root, "share", "xtb")):
        shutil.copytree(os.path.join(root, "share", "xtb"), os.path.join(dest, "share", "xtb"))
    for name in ("COPYING", "COPYING.LESSER", "README.md"):
        url = f"https://raw.githubusercontent.com/grimme-lab/xtb/v{core.XTB_VERSION}/{name}"
        try:
            urllib.request.urlretrieve(url, os.path.join(dest, name))
        except Exception as ex:  # noqa: BLE001
            print("  warning: could not fetch", name, ex)
    print("xtb bundled from", root)
    return dest


BUNDLED = ["rdkit", "numpy", "PySide6-Essentials", "PySide6", "shiboken6", "pyinstaller",
           "pyobjc-core", "pyobjc-framework-Cocoa"]


def collect_licences():
    """Copy the licence files shipped inside the installed (= bundled) packages + versions list."""
    from importlib import metadata
    dest = os.path.join(BUILD, "third_party_licenses")
    shutil.rmtree(dest, ignore_errors=True)
    os.makedirs(dest)
    versions = [f"Mol3D Bridge {VERSION}", f"Python {platform.python_version()}"]
    for name in BUNDLED:
        try:
            dist = metadata.distribution(name)
        except metadata.PackageNotFoundError:
            continue
        versions.append(f"{dist.metadata['Name']} {dist.version}  ({dist.metadata.get('License-Expression') or dist.metadata.get('License') or 'see files'})".replace("\n", " "))
        for f in dist.files or []:
            low = str(f).lower()
            base = os.path.basename(low)
            if ("licen" in low or "copying" in base or "notice" in base) and not low.endswith((".py", ".pyc", ".so", ".pyd", ".dll")):
                src = dist.locate_file(f)
                if os.path.isfile(src):
                    out = os.path.join(dest, dist.metadata["Name"], *str(f).replace("\\", "/").split("/")[-2:])
                    os.makedirs(os.path.dirname(out), exist_ok=True)
                    shutil.copy2(src, out)
    versions.append("Qt / PySide6 source: https://download.qt.io/official_releases/QtForPython/ and "
                    "https://download.qt.io/official_releases/qt/")
    with open(os.path.join(dest, "VERSIONS.txt"), "w") as fh:
        fh.write("\n".join(versions) + "\n")
    print("licences collected:", ", ".join(sorted(os.listdir(dest))))
    return dest


def pyinstaller(xtb_dir):
    env = dict(os.environ, M3B_XTB_DIR=xtb_dir or "")
    shutil.rmtree(os.path.join(DIST, NAME), ignore_errors=True)
    shutil.rmtree(os.path.join(DIST, NAME + ".app"), ignore_errors=True)
    run([sys.executable, "-m", "PyInstaller", "--noconfirm", "--clean",
         "--distpath", DIST, "--workpath", os.path.join(BUILD, "pyi"),
         os.path.join(ROOT, "packaging", "mol3d_bridge.spec")], env=env)


def frozen_exe():
    if MAC:
        return os.path.join(DIST, NAME + ".app", "Contents", "MacOS", NAME)
    return os.path.join(DIST, NAME, NAME + (".exe" if WIN else ""))


def selftest():
    out = os.path.join(BUILD, "selftest.txt")
    if os.path.exists(out):
        os.remove(out)
    env = dict(os.environ, QT_QPA_PLATFORM="offscreen", MOL3D_SELFTEST_OUT=out)
    subprocess.run([frozen_exe(), "--selftest"], env=env, check=True, timeout=600)
    with open(out) as fh:
        print(fh.read().strip())


def find_iscc():
    iscc = shutil.which("iscc") or shutil.which("ISCC")
    cands = [r"C:\Program Files (x86)\Inno Setup 6\ISCC.exe", r"C:\Program Files\Inno Setup 6\ISCC.exe",
             os.path.join(os.environ.get("LOCALAPPDATA", ""), "Programs", "Inno Setup 6", "ISCC.exe")]
    for cand in cands:
        if not iscc and os.path.exists(cand):
            iscc = cand
    return iscc


def windows_installer():
    iscc = find_iscc()
    if not iscc:
        print("Inno Setup 6 not found - skipping the installer (get it from jrsoftware.org).")
        return
    run([iscc, f"/DAppVersion={VERSION}", os.path.join(ROOT, "packaging", "windows", "installer.iss")])


def mac_dmg():
    arch = platform.machine()
    stage = os.path.join(BUILD, "dmg")
    shutil.rmtree(stage, ignore_errors=True)
    os.makedirs(stage)
    run(["ditto", os.path.join(DIST, NAME + ".app"), os.path.join(stage, NAME + ".app")])
    os.symlink("/Applications", os.path.join(stage, "Applications"))
    for f in ("LICENSE", "THIRD_PARTY_NOTICES.md"):
        shutil.copy2(os.path.join(ROOT, f), os.path.join(stage, f))
    shutil.copytree(os.path.join(ROOT, "licenses"), os.path.join(stage, "licenses"))
    dmg = os.path.join(DIST, f"Mol3DBridge-{VERSION}-macOS-{arch}.dmg")
    if os.path.exists(dmg):
        os.remove(dmg)
    run(["hdiutil", "create", "-volname", NAME, "-srcfolder", stage, "-ov", "-format", "UDZO", dmg])


def linux_tarball():
    arch = platform.machine()
    name = f"Mol3DBridge-{VERSION}-Linux-{arch}"
    out = os.path.join(DIST, name + ".tar.gz")
    with tarfile.open(out, "w:gz") as t:
        t.add(os.path.join(DIST, NAME), f"{name}/{NAME}")
        for f in ("install.sh", "uninstall.sh", "mol3d-bridge.desktop"):
            t.add(os.path.join(ROOT, "packaging", "linux", f), f"{name}/{f}")
        t.add(os.path.join(SRC, "mol3d_bridge", "resources", "icon_256.png"), f"{name}/mol3d-bridge.png")
        t.add(os.path.join(ROOT, "README.md"), f"{name}/README.md")
        t.add(os.path.join(ROOT, "LICENSE"), f"{name}/LICENSE")
        t.add(os.path.join(ROOT, "THIRD_PARTY_NOTICES.md"), f"{name}/THIRD_PARTY_NOTICES.md")
        t.add(os.path.join(ROOT, "licenses"), f"{name}/licenses")
    print("wrote", out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--with-xtb", action="store_true", help="bundle xtb 6.7.1 (Windows/Linux)")
    ap.add_argument("--no-installer", action="store_true")
    ap.add_argument("--no-selftest", action="store_true")
    ap.add_argument("--installer-only", action="store_true",
                    help="only package an existing dist/ build into the installer")
    a = ap.parse_args()
    if a.installer_only:
        if not os.path.exists(frozen_exe()):
            sys.exit("No finished build in dist/ - run the full build first.")
        windows_installer() if WIN else mac_dmg() if MAC else linux_tarball()
        return
    print(f"Building {NAME} {VERSION} for {sys.platform}/{platform.machine()}")
    version_info()
    xtb = bundle_xtb() if a.with_xtb else ""
    collect_licences()
    pyinstaller(xtb)
    if not a.no_selftest:
        selftest()
    if not a.no_installer:
        if WIN:
            windows_installer()
        elif MAC:
            mac_dmg()
        else:
            linux_tarball()
    print("done ->", DIST)


if __name__ == "__main__":
    main()
