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
//

#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorMessenger::DetectorMessenger(DetectorConstruction* geom)
  : fGeom(geom)
{
  fGeomDir = new G4UIdirectory("/my_geom/");
  fGeomDir->SetGuidance("Commands to set the geometry");

  fWorldXYCmd = new G4UIcmdWithADoubleAndUnit("/my_geom/worldXY",this);
  fWorldXYCmd->SetGuidance("Set XY half-length of the world volume.");
  fWorldXYCmd->SetParameterName("worldXY",true);
  fWorldXYCmd->SetUnitCategory("Length");
  fWorldXYCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fWorldZCmd = new G4UIcmdWithADoubleAndUnit("/my_geom/worldZ",this);
  fWorldZCmd->SetGuidance("Set Z half-length of the world volume.");
  fWorldZCmd->SetParameterName("worldZ",true);
  fWorldZCmd->SetUnitCategory("Length");
  fWorldZCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorMessenger::~DetectorMessenger()
{
  delete fGeomDir;
  delete fWorldXYCmd;
  delete fWorldZCmd;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fWorldXYCmd)
    fGeom->SetWorldXY(fWorldXYCmd->GetNewDoubleValue(newValue));

  else if (command == fWorldZCmd)
    fGeom->SetWorldZ(fWorldZCmd->GetNewDoubleValue(newValue));
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

