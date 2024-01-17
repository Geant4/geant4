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

#include "OpNoviceDetectorMessenger.hh"
#include "OpNoviceDetectorConstruction.hh"
#include "OpNoviceGDMLDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceDetectorMessenger::OpNoviceDetectorMessenger(
  G4VUserDetectorConstruction* detcon)
  : G4UImessenger()
  , fOpNoviceDetCon(detcon)
{
  fDetConDir = new G4UIdirectory("/OpNovice/DetectorConstruction/");
  fDetConDir->SetGuidance("Configuring Detector Construction");

  fVerboseCmd =
    new G4UIcmdWithABool("/OpNovice/DetectorConstruction/enableVerbose", this);
  fVerboseCmd->SetGuidance("Set flag for enabling verbose diagnostic printout");
  fVerboseCmd->SetDefaultValue(false);
  fVerboseCmd->AvailableForStates(G4State_PreInit);

  fDumpGdmlCmd =
    new G4UIcmdWithABool("/OpNovice/DetectorConstruction/dumpGdml", this);
  fDumpGdmlCmd->SetGuidance(
    "Set flag for enabling dumping the detector to a gdml file");
  fDumpGdmlCmd->SetDefaultValue(false);
  fDumpGdmlCmd->AvailableForStates(G4State_PreInit);

  fDumpGdmlFileNameCmd = new G4UIcmdWithAString(
    "/OpNovice/DetectorConstruction/dumpGdmlFileName", this);
  fDumpGdmlFileNameCmd->SetGuidance("Enter file name to dump gdml file ");
  fDumpGdmlFileNameCmd->SetDefaultValue("OpNovice_dump.gdml");
  fDumpGdmlFileNameCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceDetectorMessenger::~OpNoviceDetectorMessenger()
{
  delete fDetConDir;
  delete fVerboseCmd;
  delete fDumpGdmlCmd;
  delete fDumpGdmlFileNameCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceDetectorMessenger::SetNewValue(G4UIcommand* command,
                                            G4String newValue)
{
  auto dc1 = dynamic_cast<OpNoviceDetectorConstruction*>(fOpNoviceDetCon);
  if(dc1 != nullptr)
  {
    if(command == fVerboseCmd)
      dc1->SetVerbose(fVerboseCmd->GetNewBoolValue(newValue));
    if(command == fDumpGdmlCmd)
      dc1->SetDumpGdml(fDumpGdmlCmd->GetNewBoolValue(newValue));
    if(command == fDumpGdmlFileNameCmd)
      dc1->SetDumpGdmlFile(newValue);
  }
  else
  {
    auto dc2 = dynamic_cast<OpNoviceGDMLDetectorConstruction*>(fOpNoviceDetCon);
    if(command == fVerboseCmd)
      dc2->SetVerbose(fVerboseCmd->GetNewBoolValue(newValue));
    if(command == fDumpGdmlCmd)
      dc2->SetDumpGdml(fDumpGdmlCmd->GetNewBoolValue(newValue));
    if(command == fDumpGdmlFileNameCmd)
      dc2->SetDumpGdmlFile(newValue);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
