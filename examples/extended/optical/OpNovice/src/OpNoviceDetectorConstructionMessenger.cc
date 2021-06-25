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
// Geant4 headers
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
// project headers
#include "OpNoviceDetectorConstructionMessenger.hh"
#include "OpNoviceDetectorConstruction.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceDetectorConstructionMessenger::OpNoviceDetectorConstructionMessenger(
  OpNoviceDetectorConstruction* detcon)
  : G4UImessenger()
  , fOpNoviceDetCon(detcon)
{
  fDetConDir = new G4UIdirectory("/OpNovice/DetectorConstruction/");
  fDetConDir->SetGuidance("Configuring Detector Construction");
  //
  fverboseCmd =
    new G4UIcmdWithABool("/OpNovice/DetectorConstruction/enable_verbose", this);
  fverboseCmd->SetGuidance("Set flag for enabling verbose diagnostic printout");
  fverboseCmd->SetParameterName("enable_verbose", false);
  fverboseCmd->SetDefaultValue(false);
  fverboseCmd->AvailableForStates(G4State_PreInit);

  fdumpGdmlCmd =
    new G4UIcmdWithABool("/OpNovice/DetectorConstruction/dumpgdml", this);
  fdumpGdmlCmd->SetGuidance(
    "Set flag for enabling dumping the detector to a gdml file");
  fdumpGdmlCmd->SetParameterName("dumpgdml", false);
  fdumpGdmlCmd->SetDefaultValue(false);
  fdumpGdmlCmd->AvailableForStates(G4State_PreInit);

  fdumpGdmlFileNameCmd = new G4UIcmdWithAString(
    "/OpNovice/DetectorConstruction/DumpGDMLFileName", this);
  fdumpGdmlFileNameCmd->SetGuidance("Enter file name to dump gdml file ");
  fdumpGdmlFileNameCmd->SetParameterName("GDMLFileName", true);
  fdumpGdmlFileNameCmd->SetDefaultValue("OpNovice_dump.gdml");
  fdumpGdmlFileNameCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceDetectorConstructionMessenger::~OpNoviceDetectorConstructionMessenger()
{
  delete fDetConDir;
  delete fverboseCmd;
  delete fdumpGdmlCmd;
  delete fdumpGdmlFileNameCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceDetectorConstructionMessenger::SetNewValue(G4UIcommand* command,
                                                        G4String newValue)
{
  if(command == fverboseCmd)
    fOpNoviceDetCon->SetVerbose(fverboseCmd->GetNewBoolValue(newValue));
  if(command == fdumpGdmlCmd)
    fOpNoviceDetCon->SetDumpgdml(fdumpGdmlCmd->GetNewBoolValue(newValue));
  if(command == fdumpGdmlFileNameCmd)
    fOpNoviceDetCon->SetDumpgdmlFile(newValue);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
