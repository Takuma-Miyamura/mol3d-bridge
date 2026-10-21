; Inno Setup script for Mol3D Bridge  (compiled by packaging/build.py)
#ifndef AppVersion
  #define AppVersion "1.0.0"
#endif
#define AppName "Mol3D Bridge"
#define AppExe "Mol3D Bridge.exe"

[Setup]
AppId={{DDD1E946-5027-454E-8C55-C1299630FD7F}
AppName={#AppName}
AppVersion={#AppVersion}
AppVerName={#AppName} {#AppVersion}
AppPublisher=Takuma Miyamura
AppCopyright=Copyright (c) 2026 Takuma Miyamura
DefaultDirName={autopf}\{#AppName}
DefaultGroupName={#AppName}
DisableProgramGroupPage=yes
; installs for the current user without admin rights unless the user chooses "all users"
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
OutputDir=..\..\dist
OutputBaseFilename=Mol3DBridge-{#AppVersion}-Windows-x64-Setup
SetupIconFile=..\..\src\mol3d_bridge\resources\icon.ico
UninstallDisplayIcon={app}\{#AppExe}
UninstallDisplayName={#AppName}
LicenseFile=..\..\LICENSE
Compression=lzma2/ultra64
SolidCompression=yes
WizardStyle=modern
ChangesAssociations=yes
CloseApplications=yes

[Tasks]
Name: "desktopicon"; Description: "Create a &desktop shortcut"; GroupDescription: "Shortcuts:"
Name: "assoc"; Description: "Add Mol3D Bridge to 'Open with' for .cdx .cdxml .mol .sdf .xyz .pdb"; GroupDescription: "File types:"; Flags: unchecked

[Files]
Source: "..\..\dist\Mol3D Bridge\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs
Source: "..\..\LICENSE"; DestDir: "{app}"; DestName: "LICENSE.txt"; Flags: ignoreversion
Source: "..\..\THIRD_PARTY_NOTICES.md"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\..\licenses\*"; DestDir: "{app}\licenses"; Flags: ignoreversion

[Icons]
Name: "{autoprograms}\{#AppName}"; Filename: "{app}\{#AppExe}"
Name: "{autodesktop}\{#AppName}"; Filename: "{app}\{#AppExe}"; Tasks: desktopicon

[Registry]
Root: HKA; Subkey: "Software\Classes\Mol3DBridge.Structure"; ValueType: string; ValueName: ""; ValueData: "Chemical structure"; Flags: uninsdeletekey; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\Mol3DBridge.Structure\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: """{app}\{#AppExe}"",0"; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\Mol3DBridge.Structure\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\{#AppExe}"" ""%1"""; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.cdx\OpenWithProgids"; ValueType: string; ValueName: "Mol3DBridge.Structure"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.cdxml\OpenWithProgids"; ValueType: string; ValueName: "Mol3DBridge.Structure"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.mol\OpenWithProgids"; ValueType: string; ValueName: "Mol3DBridge.Structure"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.sdf\OpenWithProgids"; ValueType: string; ValueName: "Mol3DBridge.Structure"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.xyz\OpenWithProgids"; ValueType: string; ValueName: "Mol3DBridge.Structure"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc
Root: HKA; Subkey: "Software\Classes\.pdb\OpenWithProgids"; ValueType: string; ValueName: "Mol3DBridge.Structure"; ValueData: ""; Flags: uninsdeletevalue; Tasks: assoc

[Run]
Filename: "{app}\{#AppExe}"; Description: "Start {#AppName}"; Flags: nowait postinstall skipifsilent
