# Mol3D Bridge

A Chem3D-style 3D companion for **ChemDraw** and **ReShiki** — desktop app (Windows, macOS,
Linux), Python library and command-line tool.

Paste a drawing → 3D model → adjust it by hand → minimize it (MMFF94, MMFF94s, UFF,
GFN2-xTB, GFN-FF) and watch it fall into the minimum → copy it back as a 2D drawing or as a
rotatable 3D structure, or write ORCA / Gaussian input files.

© 2026 Takuma Miyamura — MIT licence (see `LICENSE`). Third-party components and their
licences: `THIRD_PARTY_NOTICES.md` and `licenses/`.

> ChemDraw and Chem3D are trademarks of Revvity, Inc. Mol3D Bridge is an independent project and
> is not affiliated with, sponsored or endorsed by Revvity. ReShiki is an independent open-source
> project; Mol3D Bridge contains no ReShiki code.

## Quick installation (Windows)

1. **Download** the project as a zip file (on GitHub: green **Code** button ▸ **Download ZIP**).
2. **Extract** it: right-click the zip ▸ **Extract All…** (a GitHub download unpacks to
   `mol3d-bridge-main` – use that folder wherever `mol3d-bridge` is written below).
3. Make sure **Python 3.9 or newer** is installed (python.org – 3.12 recommended; tick
   *Add python.exe to PATH* in the installer). Check in a terminal with `py --version`.
4. Double-click **`mol3d-bridge\packaging\windows\build_windows.bat`** and wait until it says
   *Done* (5–10 minutes the first time; it downloads RDKit, Qt and xtb by itself).
5. Make a desktop shortcut: open `mol3d-bridge\dist\Mol3D Bridge\`, right-click
   **`Mol3D Bridge.exe`** ▸ *Show more options* ▸ *Send to* ▸ **Desktop (create shortcut)**.
6. Start Mol3D Bridge by double-clicking the shortcut (or `Mol3D Bridge.exe` itself).

Keep the whole `dist\Mol3D Bridge` folder together – the program needs the files next to the
`.exe`. Optional: `packaging\windows\make_installer.bat` turns it into a normal setup program
(Start-menu entry, uninstaller). For the TS search, pysisyphus is installed from inside the program
(*TS search… ▸ pysisyphus program… ▸ Install pysisyphus automatically*).

## Install

| System | File | How |
|---|---|---|
| Windows 10/11 (64-bit) | `Mol3DBridge-1.3.1-Windows-x64-Setup.exe` | double-click, *Next*… Start-menu entry, optional desktop icon and "Open with" for .cdx/.cdxml/.mol/.sdf/.xyz/.pdb. No admin rights needed. Uninstall in *Settings ▸ Apps*. |
| macOS 11+ | `Mol3DBridge-1.3.1-macOS-arm64.dmg` (Apple Silicon) / `-x86_64.dmg` (Intel) | open, drag **Mol3D Bridge** to *Applications*. First start: right-click ▸ *Open* (the app is not notarised). |
| Linux x86-64 | `Mol3DBridge-1.3.1-Linux-x86_64.tar.gz` | unpack, run `./install.sh` (menu entry + `mol3d-bridge` command), `./uninstall.sh` removes it. |

Everything is included: Python, RDKit, Qt – nothing else to install. The Windows and Linux
packages also include **xtb 6.7.1**, so GFN2-xTB works immediately. On macOS install xtb once
with `conda install -c conda-forge xtb` and pick it in *Calculate ▸ xtb settings…*.

Windows may show a SmartScreen notice for a new, unsigned program: *More info ▸ Run anyway*.

## Transition-state search (pysisyphus)

**Calculate ▸ Transition-state search…** (Ctrl+T, toolbar *TS search…*) opens a separate window:

* **Reactant | Product side by side** – each with *Paste* (ChemDraw / ReShiki), *From main window*,
  *Open…*, a SMILES box, *Pre-optimize* (method in the toolbar: MMFF94 … GFN2-xTB) and *Edit atoms*.
* **Atom numbers** (Ctrl+L) shows the numbering on every view. Below the structures the window
  says which bonds break and form with the current numbering.
* **Auto-number product** renumbers the product so each atom matches the reactant: common
  substructure + optimal assignment, choosing the numbering with the fewest broken/formed bonds,
  then the least atomic motion (bond-forming atoms close in the reactant). Tested on Diels–Alder,
  Claisen, Cope, SN2, keto–enol, ring opening and HCN→HNC with scrambled input numbering.
  **Other numbering** cycles through other, equally simple numberings (e.g. Claisen [3,3] vs [1,3]).
* **Manual relabelling:** *Relabel by click* (click an atom, type its new number – the two atoms
  swap numbers) or select two atoms and *Swap 2 selected*. *Align product onto reactant*
  superposes the two.
* **Calculation tab:** every setting of your example input (dlc, growing string with 18 nodes and
  climbing image, string optimizer, RS-I-RFO TS optimization with Hessian, Euler-PC IRC,
  preopt/endopt, xtb with pal/etemp/charge/multiplicity, ALPB solvent) as a form, with the
  generated `ts.yaml` shown live – or edited by hand.
* **Monitor tab:** the pysisyphus output live, the current stage (pre-optimization → string →
  TS optimization → IRC → end optimizations), and the newest reaction path animated next to its
  energy profile.
* **Results tab:** the **transition state** (forming/breaking bonds dashed, *Play IRC* animation,
  imaginary frequency) and the **two minima** from the IRC, each labelled reactant- or product-side,
  with barriers in kJ/mol. Every result can be copied to ChemDraw/ReShiki (2D or 3D), saved,
  turned into an ORCA/Gaussian input or opened in the main window. *Load a previous run folder…*
  reopens saved runs.
* **Results folder (automatic):** every run – finished, stopped or failed – is saved in its own
  folder, by default `Documents\Mol3D Bridge\TS results\2026-10-20_14-15_<run name or formula>\`:
  `summary.txt` (energies, relative energies and barriers in kJ/mol and kcal/mol, imaginary
  frequency, bonds formed/broken, the input), `TS.xyz`, `minimum_reactant-side.xyz`,
  `minimum_product-side.xyz`, `results.sdf` (all three with energies), `IRC_path.xyz` (the whole
  path as a movie), the input structures, `ts.yaml`, `pysis.log` and – optionally – every file
  pysisyphus wrote (`pysisyphus\`). Choose the folder, a run name, and whether to keep all files in
  the *Results* box of the Calculation tab; *TS results folder* in the toolbar opens it.

**Setup once:** *pysisyphus program… ▸ Install pysisyphus automatically* – sets up a private
Python 3.12 (checksum-verified download) with pysisyphus 1.0.0 inside Mol3D Bridge's data folder;
no administrator rights, nothing else on the PC changes (≈300 MB download, ≈1.5 GB on disk, a
few minutes). An existing `pip install pysisyphus` is detected automatically instead (PATH, the
Windows Python install manager, python.org and conda installations), or choose `pysis` with
*Browse…*. xtb is taken from *Calculate ▸ xtb settings*.
pysisyphus is free software under the GNU GPL v3 (not included in Mol3D Bridge; installed
separately and run as its own program – see `THIRD_PARTY_NOTICES.md`). Please cite it:
J. Steinmetzer, S. Kupfer, S. Gräfe, "pysisyphus: Exploring potential energy surfaces in ground and
excited states", *Int. J. Quantum Chem.* **2021**, 121, e26390. doi:10.1002/qua.26390

From Python: `from mol3d_bridge import ts` – `ts.map_atoms(reactant, product)`,
`ts.yaml_input(ts.TSSettings(...))`, `ts.write_inputs(...)`, `ts.run_pysis(...)`,
`ts.collect_results(...)`.

## Use it from Python (3.9 or newer)

```bash
pip install "git+https://github.com/Takuma-Miyamura/mol3d-bridge"          # API + m3b command
pip install "mol3d-bridge[gui] @ git+https://github.com/Takuma-Miyamura/mol3d-bridge"   # + desktop window
```

```python
from mol3d_bridge import Molecule3D

m = Molecule3D.from_smiles("CC(=O)N[C@@H](Cc1ccccc1)C(=O)OC")   # or from_file("x.cdx"), from_molblock(...)
m.minimize("MMFF94")              # "MMFF94s", "UFF", "GFN2-xTB", "GFN-FF" (xtb_path=... if xtb is not on PATH)
print(m.energy, m.stereo_ok())    # kcal/mol, stereochemistry unchanged?

conformers = m.conformer_search(n=30)        # unique conformers, lowest energy first
best = conformers[0]
best.set_internal([3, 4, 5, 6], 180.0)       # distance / angle / dihedral (0-based atoms)
best.minimize("MMFF94", fixed=[3, 4, 5, 6])  # hold atoms fixed
best.save("best.xyz")                        # .xyz .mol .sdf .pdb
best.save("best.inp", keywords="! r2SCAN-3c Opt Freq")              # ORCA
best.save("best.gjf", route="#p opt freq m062x/6-311+g(d,p)")       # Gaussian
frames = Molecule3D.from_smiles("CCO").trajectory()   # every minimization step (N x 3 arrays)
```

`m.mol` is the RDKit molecule (with hydrogens and the 3D conformer) for anything else.

### Command line: `m3b`

```bash
m3b "CC(C)[C@@H]1CC[C@@H](C)C[C@H]1O" -o menthol.xyz                # build + MMFF94, write XYZ
m3b drawing.cdx --method GFN2-xTB --orca job.inp --keywords "! r2SCAN-3c Opt Freq"
m3b input.mol --conformers 30 -o conformers.sdf                     # all conformers with energies
m3b "c1ccccc1O" --gaussian phenol.gjf --route "#p opt freq m062x/6-311+g(d,p)" --nprocs 8
m3b --help
```

### Python versions

| Python | Works | Notes |
|---|---|---|
| 3.9 | yes | RDKit 2025.09: reads ChemDraw files; copies to ChemDraw use MDL molfiles |
| 3.10 – 3.14 | yes | newest RDKit: also writes native ChemDraw (CDX) data — recommended |

## Building the installers

**Easiest – on GitHub (free, all platforms at once):** put this folder in a GitHub repository.
The included workflow (`.github/workflows/build.yml`) builds the Windows installer, both macOS
DMGs and the Linux package on every push (download them under *Actions ▸ latest run ▸
Artifacts*). Push a tag such as `v1.0.2` and they are published as a GitHub Release.

**On your own computer** (builds the installer for that system):

* Windows: install Python 3.9 or newer — 3.12 recommended (python.org, tick *Add to PATH*) — and Inno Setup 6
  (jrsoftware.org), then double-click `packaging\windows\build_windows.bat`.
* macOS: `packaging/macos/build_mac.sh`
* Linux: `packaging/linux/build_linux.sh`

Each build runs a self-test of the finished program (`"Mol3D Bridge" --selftest`) before
packaging. Output goes to `dist/`.

**For development:** `pip install -e ".[gui]"` then `mol3d-bridge` (or `python -m mol3d_bridge`).

## What's new in 1.3.1

* Quick installation guide; pysisyphus licence (GPL v3) and citations checked and stated in the
  README, `THIRD_PARTY_NOTICES.md`, the pysisyphus dialog and every `summary.txt`.

## What's new in 1.3.0

* **TS results are stored automatically** in a results folder you choose (default
  `Documents\Mol3D Bridge\TS results`), one clearly named folder per run with a summary, the TS,
  both minima, the IRC path, an SD file with energies, the inputs and the full pysisyphus output.
  Stopped and failed runs are kept too (marked in the folder name). *Open a saved run…* reopens
  any of them.

## What's new in 1.2.4

* TS window: the **active structure** (the one you last clicked) has a blue frame and a blue
  title; Copy, Paste, Swap and Relabel act on it.
* **Copy 3D / Paste** in both structures (also Ctrl+C / Ctrl+V): copies the exact 3D structure
  with hydrogens and atom numbers – between reactant and product, to the main window, or to
  ChemDraw / ReShiki. **→ Copy to product / ← Copy to reactant** does it in one click.
* Numbering check: the bonds that change are drawn dashed in orange (formed on the reactant side,
  broken on the product side) and the warning stays a warning – the calculation no longer asks
  to renumber; you swap atoms as you see fit.

## What's new in 1.2.3

* Windows: pysisyphus stopped with `OSError: [WinError 1314] A required privilege is not held`
  (it creates a symbolic link, which Windows only allows with administrator rights or Developer
  Mode). Mol3D Bridge now starts pysisyphus through a small launcher that skips that optional link.

## What's new in 1.2.2

* **One-click pysisyphus installation** in the TS window, like the xtb download, and automatic
  detection of pysisyphus installed with pip (also when its folder is not on PATH).

## What's new in 1.2.1

* Windows: pysisyphus crashed with `UnicodeEncodeError: 'charmap' codec` when it printed
  formulas such as C₅₄H₃₆. pysis is now always started in UTF-8 mode.

## What's new in 1.2.0

* **Transition-state search with pysisyphus** in its own window: reactant and product side by side,
  atom-number display, automatic and manual renumbering, the full pysisyphus input as a form and
  as YAML, live output monitor with the growing reaction path, and the TS plus both IRC minima as
  results.

## What's new in 1.1.0

* **Python API** (`from mol3d_bridge import Molecule3D`) and **`m3b` command-line tool**:
  build, minimize (MMFF94/MMFF94s/UFF/GFN2-xTB/GFN-FF), conformer search, set internal
  coordinates, export XYZ/MOL/SDF/PDB/ORCA/Gaussian — no window needed.
* Runs on **Python 3.9 and newer** (tested on 3.9, 3.12, 3.13). The desktop window is an optional
  extra (`pip install "mol3d-bridge[gui]"`).
* Licensing prepared for public release: author and copyright, full licence texts in `licenses/`,
  licence files of the exact bundled package versions collected into every installer
  (`third_party_licenses/`), Qt/xtb source information, trademark notice.

## What's new in 1.0.2

* **ChemDraw asked "Change Settings?" and, with *Don't Change Settings*, wanted a document 9212
  pages wide.** The ChemDraw data declared a default bond length of 0.0002 pt (RDKit writes the
  value in raw 1/65536-pt units), so ChemDraw scaled the structure up ~67,000x. Copies now declare
  the real bond length and are drawn at it.
* New *Settings > ChemDraw bond length for copies...*: 14.4 pt (ACS Document 1996, default), 30 pt
  (ChemDraw "New Document") or custom - set it to your document style.

## What's new in 1.0

* Real desktop application with installers, icon, About box, file associations, crash log
  (*Help ▸ Open log / data folder*).
* **Paste no longer minimizes.** The structure is built in 3D and shown oriented like your
  drawing; press **Ctrl+M** when you want to minimize (or switch on *Settings ▸ Minimize
  automatically after pasting*).
* **Strained / unusual structures fixed** (e.g. cycloparaphenylenes such as
  `C1(C2=CC=C(C3=CC=C(C4=CC=C5C=C4)C=C3)C=C2)=CC=C(C6=CC=C5C=C6)C=C1`): 3D building now tries
  ETKDG, then more tolerant distance geometry, each with a time limit – it can no longer hang –
  and as a last resort starts from the drawing itself.
* **ChemDraw / ReShiki pre-cleaned drawings:** coordinates in other units (points, 1/65536 pt
  from 3D CDX) are rescaled, drawings that only carry a few depth hints are no longer mistaken
  for 3D models, invalid structures give a clear message, and pasting while a calculation runs is
  refused instead of colliding with it.

## Workflow

| Step | ChemDraw / ReShiki | Mol3D Bridge |
|---|---|---|
| 1 | select the structure, **Ctrl+C** | **Ctrl+V** → hydrogens added, 3D model built and shown like your drawing - **not** minimized |
| 2 | | rotate / measure / **Ctrl+M** minimize again / **Ctrl+K** conformer search |
| 3 | | rotate to the view you want (or **Ctrl+0** = best flat view), **Ctrl+C** |
| 4 | **Ctrl+V** → 2D drawing in that orientation, wedges/hashes from real depth | |
| 5 | | **Ctrl+E** → ORCA input, **Ctrl+Shift+G** → Gaussian input (preview, copy, save) |

### Adjusting the starting structure (Edit atoms mode, key **E**)

| Action | How |
|---|---|
| Rotate view | left-drag (Rotate mode) or right-drag (Edit mode); wheel = zoom |
| Select atoms | click, Shift+click to add, drag a box on empty space, double-click = whole molecule, Ctrl+A |
| Select a substituent | select the atom that stays, then the atom of the group → **Ctrl+Shift+S** |
| Move | drag a selected atom (screen plane); **Alt+drag** = depth; arrow keys / PgUp/PgDn nudge |
| Rotate a group | **Shift+drag** on a selected atom (rotates around the group's centre) |
| Exact values | select 2/3/4 atoms in order → **Ctrl+G** sets distance / angle / dihedral (last atom's side moves) |
| Minimize from your geometry | **Ctrl+M** (Edit ▸ *Hold selected atoms fixed* freezes them; *Relax after each move* minimizes automatically) |
| Undo | **Ctrl+Z** (every drag is one undo step) |

Attached hydrogens move with their atoms (toggle in the Edit menu). After a move the
force-field energy of the unminimized geometry is shown, and you are warned if an edit inverted a
stereocentre or C=C relative to the drawing. *Rebuild* (Ctrl+B) throws the edits away and re-embeds
the pasted drawing.

### Copy keeping 3D (like Chem3D's copy to ChemDraw)

* **Ctrl+D – for ChemDraw:** the drawing is the current view (with wedges), and every atom also
  carries its 3D position (CDX `3DPosition`), as ChemDraw/Chem3D store 3D structures. Rotating it
  in 3D in ChemDraw keeps the real geometry. (*Settings ▸ ChemDraw 3D copy: use CDX 3D positions*
  off = a 3D molfile instead.)
* **Ctrl+Shift+D – for ReShiki:** a 3D molfile (x/y = current view, z = depth). ReShiki keeps the
  z value as atom depth; its 3D tool then starts from this geometry. (No CDX is put on the clipboard
  here on purpose: ReShiki reads CDX first but keeps depth only from molfiles.)
* Both are checked when copying: stereo decoded from the written 3D coordinates must match the model.
* If ChemDraw ever shows the mirror image after a 3D rotation, tick *Settings ▸ ChemDraw 3D copy:
  invert depth*.
* Pasting **into** Mol3D Bridge: 3D positions inside ChemDraw data are ignored (the drawing and its
  wedges are used and the model is rebuilt) unless a real 3D molfile is on the clipboard too.

Other copies: **Ctrl+Shift+C** clean 2D depiction (stereo from the 3D model),
**Ctrl+Alt+C** 3D molfile, Edit ▸ Copy XYZ, **Ctrl+Shift+I** picture of the 3D view.

## Clipboard formats

Read (best first) and written, matching what ChemDraw writes and ReShiki reads/writes:

| | Windows | macOS | Linux |
|---|---|---|---|
| CDX | `ChemDraw Interchange Format` | `com.revvity.chemdraw.cdx-clipboard` (+ PerkinElmer/CambridgeSoft aliases) | `chemical/x-cdx` |
| CDXML | `ChemDraw XML`, `chemical/x-cdxml` (read) | `com.revvity.cdxml` … (read) | `chemical/x-cdxml` (read) |
| Molfile | `MDLCT` (length-prefixed), `chemical/x-mdl-molfile` | `com.mdli.molfile` | `chemical/x-mdl-molfile` |
| SMILES | `SMILES` | `org.opensmiles.smiles` | `chemical/x-daylight-smiles` |
| Text | SMILES (or molfile, see Settings), also accepts pasted MOL / CDXML / InChI / XYZ text | | |

CDX reading/writing uses RDKit's ChemDraw support (`Chem.HasChemDrawCDXSupport()`).
If your RDKit lacks it, MOL + SMILES are still used, and in ChemDraw
*Edit ▸ Copy As ▸ MOL Text* always works.

## Methods

| Method | Needs | Typical time (30 atoms) | Notes |
|---|---|---|---|
| MMFF94 / MMFF94s | nothing (RDKit) | < 0.2 s | default; organic molecules |
| UFF | nothing | < 0.2 s | automatic fallback when MMFF lacks parameters (metals) |
| GFN2-xTB | `xtb` program | ~1–3 s | semi-empirical QM; much better starting geometries for DFT, charged species, metals |
| GFN-FF | `xtb` program | < 1 s | xtb's general force field |

**xtb setup:** Calculate ▸ *xtb settings…* → *Download xtb 6.7.1* (official GitHub release,
SHA-256 verified; Windows/Linux). On macOS: `conda install -c conda-forge xtb`, then *Browse…*.
Settings: threads, ALPB implicit solvent, convergence level, charge / unpaired electrons
(auto from the structure). Frozen atoms (*Hold selected atoms fixed*) are passed to xtb as `$fix`.
Cite: Bannwarth, Ehlert, Grimme, *J. Chem. Theory Comput.* 2019, 15, 1652 (doi:10.1021/acs.jctc.8b01176).

**Animated optimisation:** every minimization is shown step by step with a live energy curve
(ΔE vs step, log scale). Force fields use FIRE with a small step cap (0.02 Å per atom per step,
Settings ▸ *Force-field step size*) followed by a BFGS polish, so the result is the true minimum;
xtb steps are streamed from `xtbopt.log` and interpolated. The playback is timed by the wall
clock (Settings ▸ *Animation length*, default 5 s; geometries that barely move take ~1 s). You can
rotate/zoom while it runs; **Esc / ■ Stop** skips the animation (or stops xtb and keeps the last
step). *Settings ▸ Animate optimizations* off = instant.

**Conformer search** = N ETKDG starts minimized with MMFF (duplicates pruned at 0.25 Å RMSD,
sorted by energy). With GFN2-xTB selected you are offered *Refine all conformers with xtb*,
which re-optimizes and re-ranks them.

## Gaussian input

File ▸ *Gaussian input…* (Ctrl+Shift+G) or toolbar *Gaussian…*: route presets (B3LYP-D3BJ,
ωB97X-D, M06-2X, PBE0 + SMD, TS search, custom), charge/multiplicity check, `%nprocshared`,
`%mem`, `%chk` named after the file, text after the coordinates (ModRedundant, gen basis),
editable preview, copy or save as `.gjf` / `.com`.

## Notes / limits

* Atoms can be moved but not added/deleted: change connectivity in ChemDraw/ReShiki and paste again.
* Stereo is checked both ways (drawing → 3D, 3D → drawing); a warning appears in the status bar
  if a projection makes a C=C or stereocentre ambiguous — rotate or use *clean 2D*.
* Atom indices (View ▸ Labels ▸ Atom index) are 0-based, the same numbering ORCA uses
  for `%geom Constraints`.
* MM geometries are a starting point; optimise with ORCA (e.g. r2SCAN-3c or GFN2-xTB) afterwards.

Files: `mol3d_bridge.py` (GUI), `m3b_core.py` (chemistry, scriptable without GUI),
`m3b_clipboard.py` (ChemDraw/ReShiki clipboard bridge).
