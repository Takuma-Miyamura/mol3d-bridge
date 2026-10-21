"""Starts the desktop window; explains what to install if the GUI libraries are missing."""
import sys


def main():
    try:
        from .app import main as app_main
    except ImportError as ex:
        if "PySide6" in str(ex) or "shiboken" in str(ex):
            sys.stderr.write("The Mol3D Bridge window needs Qt for Python:\n"
                             '    pip install "mol3d-bridge[gui]"\n'
                             "(The Python API and the m3b command work without it.)\n")
            sys.exit(1)
        raise
    app_main()


if __name__ == "__main__":
    main()
