//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: XrayFluoDetectorMessenger.cc
// GEANT4 tag $Name: xray_fluo-V03-02-00
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
//
// -------------------------------------------------------------------

#include "XrayFluoDetectorMessenger.hh"
#include "XrayFluoDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoDetectorMessenger::XrayFluoDetectorMessenger(XrayFluoDetectorConstruction * Det)
:Detector(Det)
{ 
   detDir = new G4UIdirectory("/apparate/");
  detDir->SetGuidance("detector control.");

  UpdateCmd = new G4UIcmdWithoutParameter("/apparate/update",this);
  UpdateCmd->SetGuidance("Update apparate geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(Idle);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

 XrayFluoDetectorMessenger::~XrayFluoDetectorMessenger()
{
  delete UpdateCmd;
  delete detDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
 if( command == UpdateCmd )
   { Detector->UpdateGeometry(); }
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


