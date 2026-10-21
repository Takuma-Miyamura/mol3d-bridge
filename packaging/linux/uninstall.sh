#!/usr/bin/env bash
rm -rf "$HOME/.local/opt/mol3d-bridge" "$HOME/.local/bin/mol3d-bridge" \
       "$HOME/.local/share/applications/mol3d-bridge.desktop" \
       "$HOME/.local/share/icons/hicolor/256x256/apps/mol3d-bridge.png"
echo "Mol3D Bridge removed (settings stay in ~/.config/Mol3DBridge and ~/.local/share/Mol3DBridge)."
