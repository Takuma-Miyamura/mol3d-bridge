"""
Clipboard bridge for ChemDraw and ReShiki.

Format names follow what ChemDraw writes and what ReShiki reads/writes
(see ReShiki native/windows/src/clipboard.rs and src/clipboard/priority.rs):

  kind      Windows (registered format)          macOS (UTI)                               Linux (MIME)
  cdx       "ChemDraw Interchange Format"        com.revvity.chemdraw.cdx-clipboard, ...   chemical/x-cdx
  cdxml     "ChemDraw XML", chemical/x-cdxml     com.revvity.cdxml, ...                    chemical/x-cdxml
  molfile   "MDLCT" (length-prefixed lines),     com.mdli.molfile                          chemical/x-mdl-molfile
            chemical/x-mdl-molfile (plain text)
  smiles    "SMILES"                             org.opensmiles.smiles                     chemical/x-daylight-smiles
  text      CF_UNICODETEXT                       public.utf8-plain-text                    text/plain
"""
from __future__ import annotations

import sys

from PySide6.QtCore import QMimeData
from PySide6.QtWidgets import QApplication

from .core import mdlct_encode

WIN = sys.platform.startswith("win")
MAC = sys.platform == "darwin"

READ_ORDER = ("cdx", "cdxml", "molfile", "smiles")

WIN_FORMATS = {
    "cdx": ["ChemDraw Interchange Format"],
    "cdxml": ["ChemDraw XML", "chemical/x-cdxml"],
    "molfile": ["MDLCT", "chemical/x-mdl-molfile"],
    "smiles": ["SMILES"],
}
MAC_TYPES = {
    "cdx": ["com.revvity.chemdraw.cdx-clipboard", "com.perkinelmer.chemdraw.cdx-clipboard",
            "com.cambridgesoft.cdx", "com.revvity.cdx", "com.perkinelmer.cdx"],
    "cdxml": ["com.revvity.cdxml", "com.perkinelmer.cdxml", "com.cambridgesoft.cdxml",
              "com.perkinelmer.chemdraw.cdxml", "public.cdxml"],
    "molfile": ["com.mdli.molfile"],
    "smiles": ["org.opensmiles.smiles"],
}
LINUX_MIMES = {
    "cdx": ["chemical/x-cdx", "application/x-chemdraw"],
    "cdxml": ["chemical/x-cdxml"],
    "molfile": ["chemical/x-mdl-molfile"],
    "smiles": ["chemical/x-daylight-smiles"],
}

# ChemDraw on macOS writes the first CDX type; offer the three ReShiki offers.
MAC_CDX_WRITE = MAC_TYPES["cdx"][:3]


def _wmime(name: str) -> str:
    return f'application/x-qt-windows-mime;value="{name}"'


def _appkit():
    try:
        from AppKit import NSPasteboard  # pyobjc-framework-Cocoa
        from Foundation import NSData
        return NSPasteboard, NSData
    except Exception:
        return None, None


# ------------------------------------------------------------------------------ read
def read_items():
    """Return [(kind, bytes, source_name)] in preferred order (best first)."""
    items = []
    if MAC:
        NSPasteboard, _ = _appkit()
        if NSPasteboard is not None:
            pb = NSPasteboard.generalPasteboard()
            for kind in READ_ORDER:
                for t in MAC_TYPES[kind]:
                    d = pb.dataForType_(t)
                    if d is not None and d.length() > 0:
                        items.append((kind, bytes(d), t))
                        break
            s = pb.stringForType_("public.utf8-plain-text")
            if s:
                items.append(("text", str(s).encode("utf-8"), "public.utf8-plain-text"))
            return items

    md = QApplication.clipboard().mimeData()
    if md is None:
        return items
    for kind in READ_ORDER:
        names = [_wmime(n) for n in WIN_FORMATS[kind]] if WIN else LINUX_MIMES[kind]
        for name in names:
            try:
                data = bytes(md.data(name))
            except Exception:
                data = b""
            if data and data.strip(b"\x00"):
                items.append((kind, data, name))
                break
    if md.hasText() and md.text().strip():
        items.append(("text", md.text().encode("utf-8"), "text"))
    return items


def available_summary() -> str:
    return ", ".join(f"{k} ({s})" for k, _, s in read_items()) or "nothing usable"


# ----------------------------------------------------------------------------- write
def write_items(cdx: bytes | None = None, molfile: str | None = None,
                smiles: str | None = None, text: str | None = None):
    """Put one structure on the clipboard in every format ChemDraw/ReShiki understand."""
    if MAC:
        NSPasteboard, NSData = _appkit()
        if NSPasteboard is not None:
            payload = []
            if cdx:
                payload += [(t, cdx) for t in MAC_CDX_WRITE]
            if molfile:
                payload.append(("com.mdli.molfile", molfile.encode("ascii", "replace")))
            if smiles:
                payload.append(("org.opensmiles.smiles", smiles.encode("utf-8")))
            if text:
                payload.append(("public.utf8-plain-text", text.encode("utf-8")))
            pb = NSPasteboard.generalPasteboard()
            pb.clearContents()
            pb.declareTypes_owner_([t for t, _ in payload], None)
            for t, b in payload:
                pb.setData_forType_(NSData.dataWithBytes_length_(b, len(b)), t)
            return "macOS pasteboard"

    md = QMimeData()
    if WIN:
        if cdx:
            md.setData(_wmime("ChemDraw Interchange Format"), cdx)
        if molfile:
            md.setData(_wmime("MDLCT"), mdlct_encode(molfile))
            md.setData(_wmime("chemical/x-mdl-molfile"), molfile.encode("ascii", "replace"))
        if smiles:
            md.setData(_wmime("SMILES"), smiles.encode("utf-8"))
    else:
        if cdx:
            md.setData("chemical/x-cdx", cdx)
        if molfile:
            md.setData("chemical/x-mdl-molfile", molfile.encode("ascii", "replace"))
        if smiles:
            md.setData("chemical/x-daylight-smiles", smiles.encode("utf-8"))
    if text:
        md.setText(text)
    QApplication.clipboard().setMimeData(md)
    return "clipboard"


def write_text(text: str):
    QApplication.clipboard().setText(text)
