@echo off
REM Turns an existing build (dist\Mol3D Bridge\) into the Windows setup .exe - no rebuild.
REM Installs Inno Setup 6 for you (via winget) if it is missing.
setlocal
cd /d "%~dp0\..\.."
if not exist ".buildenv\Scripts\python.exe" (
  echo No build environment found - run build_windows.bat first.
  pause & exit /b 1
)
if not exist "dist\Mol3D Bridge\Mol3D Bridge.exe" (
  echo No finished build in dist\Mol3D Bridge - run build_windows.bat first.
  pause & exit /b 1
)
call :find_iscc
if not defined ISCC (
  echo Inno Setup 6 is not installed. Installing it now with winget...
  winget install -e --id JRSoftware.InnoSetup --scope user --accept-package-agreements --accept-source-agreements
  call :find_iscc
)
if not defined ISCC (
  echo.
  echo  Could not install Inno Setup automatically.
  echo  Download it from https://jrsoftware.org/isdl.php , install it, then run this file again.
  pause & exit /b 1
)
.buildenv\Scripts\python packaging\build.py --installer-only || goto :err
echo.
echo  Installer created:
dir /b dist\*Setup.exe
pause
exit /b 0

:find_iscc
set "ISCC="
for %%P in ("%ProgramFiles(x86)%\Inno Setup 6\ISCC.exe" "%ProgramFiles%\Inno Setup 6\ISCC.exe" "%LOCALAPPDATA%\Programs\Inno Setup 6\ISCC.exe") do (
  if not defined ISCC if exist %%P set "ISCC=%%~P"
)
exit /b 0

:err
echo.
echo  Packaging failed - see the messages above.
pause
exit /b 1
