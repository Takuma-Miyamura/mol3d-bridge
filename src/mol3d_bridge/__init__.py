"""
Mol3D Bridge – a Chem3D-style 3D companion for ChemDraw and ReShiki.

Python API (no window needed):

    from mol3d_bridge import Molecule3D
    m = Molecule3D.from_smiles("CCO"); m.minimize("MMFF94"); m.save("ethanol.xyz")
"""
__version__ = "1.3.1"
__author__ = "Takuma Miyamura"
APP_NAME = "Mol3D Bridge"
APP_ID = "Mol3DBridge.Mol3DBridge"

__all__ = ["Molecule3D", "METHODS", "__version__"]


def __getattr__(name):  # lazy: importing the package (e.g. for the GUI) stays light
    if name in ("Molecule3D", "METHODS"):
        from . import api
        return getattr(api, name)
    raise AttributeError(name)
