@echo off
REM Builds "Mol3D Bridge" and its installer on Windows.
REM Needs: Python 3.9 or newer (python.org) and Inno Setup 6 (jrsoftware.org).
REM Python 3.12 is preferred when several are installed (newest RDKit = best ChemDraw support).
setlocal
cd /d "%~dp0\..\.."

REM ---- find Python 3.9 - 3.14 (3.12 preferred) -------------------------------------
set "PY="
for %%V in (3.12 3.13 3.11 3.10 3.14 3.9) do (
  if not defined PY (
    py -%%V -c "import sys" >nul 2>&1 && set "PY=py -%%V"
  )
)
if not defined PY (
  python -c "import sys; sys.exit(0 if (3,9) <= sys.version_info[:2] <= (3,14) else 1)" >nul 2>&1 && set "PY=python"
)
if not defined PY (
  echo.
  echo  Python 3.9 or newer is needed, but no suitable version was found.
  echo  Install Python 3.12 from https://www.python.org/downloads/
  echo  ^(or in a terminal:  winget install Python.Python.3.12 ^)
  echo  and run this file again.
  echo.
  pause
  exit /b 1
)
echo Using %PY%:
%PY% --version

REM ---- (re)create the build environment if it was made with an old Python -----------
if exist .buildenv\Scripts\python.exe (
  .buildenv\Scripts\python -c "import sys; sys.exit(0 if sys.version_info >= (3,9) else 1)" >nul 2>&1 || (
    echo Removing the old build environment ^(made with an older Python^)...
    rmdir /s /q .buildenv
  )
)
if not exist .buildenv\Scripts\python.exe (
  %PY% -m venv .buildenv || goto :err
)

.buildenv\Scripts\python -m pip install --upgrade pip
.buildenv\Scripts\python -m pip install ".[build]" || goto :err
.buildenv\Scripts\python packaging\build.py --with-xtb || goto :err
echo.
if exist "dist\*Setup.exe" (
  echo  Done. Installer:
  dir /b dist\*Setup.exe
  echo  ^(in the "dist" folder^)
) else (
  echo  Done. The program is built and tested:  dist\Mol3D Bridge\Mol3D Bridge.exe
  echo  You can already use it from there ^(portable^).
  echo  To make the setup .exe, run packaging\windows\make_installer.bat
  echo  ^(it installs Inno Setup for you if needed - no rebuild^).
)
pause
exit /b 0

:err
echo.
echo  Build failed - see the messages above.
pause
exit /b 1
