//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************

#ifndef G4GVTKMESSENGER_HH
#define G4GVTKMESSENGER_HH 1

#include "G4UImessenger.hh"

#include <vector>

class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcommand;
class G4UIcmdWithoutParameter;

class G4VtkMessenger : public G4UImessenger
{
  public:
    static G4VtkMessenger* GetInstance();
    ~G4VtkMessenger() override;

    G4String GetCurrentValue(G4UIcommand* command) override;
    void SetNewValue(G4UIcommand* command, G4String newValue) override;

  private:
    static G4VtkMessenger* fpInstance;
    G4VtkMessenger();

    G4UIdirectory* fpDirectory;
    G4UIcommand* fpCommandClearNonG4;
    G4UIcommand* fpCommandExport;
    G4UIcommand* fpCommandExportCutter;
    G4UIcommand* fpCommandDebugPrint;
    G4UIcommand* fpCommandInteractorStart;

    G4UIdirectory* fpDirectorySet;
    G4UIcommand* fpCommandPolyhedronPipeline;
    G4UIcmdWithABool* fpCommandWarnings;
    G4UIcmdWithABool* fpCommandHUD;
    G4UIcmdWithABool* fpCameraOrientation;
    G4UIcommand* fpCommandClipper;
    G4UIcommand* fpCommandCutter;
    G4UIcommand* fpCommandShadow;

    G4UIdirectory* fpDirectoryAdd;
    G4UIcommand* fpCommandImageOverlay;
    G4UIcommand* fpCommandGeometryOverlay;
};

#endif