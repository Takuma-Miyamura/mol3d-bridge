#!/usr/bin/env bash
# Installs Mol3D Bridge for the current user (no root needed). Run from the unpacked folder.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
DEST="${HOME}/.local/opt/mol3d-bridge"
rm -rf "$DEST"; mkdir -p "$DEST"
cp -a "$HERE/Mol3D Bridge/." "$DEST/"
mkdir -p "$HOME/.local/share/applications" "$HOME/.local/share/icons/hicolor/256x256/apps" "$HOME/.local/bin"
cp "$HERE/mol3d-bridge.png" "$HOME/.local/share/icons/hicolor/256x256/apps/mol3d-bridge.png"
sed "s#@INSTALL_DIR@#$DEST#" "$HERE/mol3d-bridge.desktop" > "$HOME/.local/share/applications/mol3d-bridge.desktop"
ln -sf "$DEST/Mol3D Bridge" "$HOME/.local/bin/mol3d-bridge"
command -v update-desktop-database >/dev/null && update-desktop-database "$HOME/.local/share/applications" || true
echo "Installed. Start 'Mol3D Bridge' from your application menu, or run: mol3d-bridge"
