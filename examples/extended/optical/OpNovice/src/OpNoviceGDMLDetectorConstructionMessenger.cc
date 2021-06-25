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
#include "OpNoviceGDMLDetectorConstructionMessenger.hh"
#include "OpNoviceGDMLDetectorConstruction.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
OpNoviceGDMLDetectorConstructionMessenger::
  OpNoviceGDMLDetectorConstructionMessenger(
    OpNoviceGDMLDetectorConstruction* detcon)
  : G4UImessenger()
  , fOpNoviceDetCon(detcon)
{
  fDetConDir = new G4UIdirectory("/OpNovice/DetectorConstruction/");
  fDetConDir->SetGuidance("Configuring Detector Construction");
  //
  verboseCmd =
    new G4UIcmdWithABool("/OpNovice/DetectorConstruction/enable_verbose", this);
  verboseCmd->SetGuidance("Set flag for enabling verbose diagnostic printout");
  verboseCmd->SetParameterName("enable_verbose", false);
  verboseCmd->SetDefaultValue(false);
  verboseCmd->AvailableForStates(G4State_PreInit);

  dumpGdmlCmd =
    new G4UIcmdWithABool("/OpNovice/DetectorConstruction/dumpgdml", this);
  dumpGdmlCmd->SetGuidance(
    "Set flag for enabling dumping the detector to a gdml file");
  dumpGdmlCmd->SetParameterName("dumpgdml", false);
  dumpGdmlCmd->SetDefaultValue(false);
  dumpGdmlCmd->AvailableForStates(G4State_PreInit);

  dumpGdmlFileNameCmd = new G4UIcmdWithAString(
    "/OpNovice/DetectorConstruction/DumpGDMLFileName", this);
  dumpGdmlFileNameCmd->SetGuidance("Enter file name to dump gdml file ");
  dumpGdmlFileNameCmd->SetParameterName("GDMLFileName", true);
  dumpGdmlFileNameCmd->SetDefaultValue("OpNovice_dump.gdml");
  dumpGdmlFileNameCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
OpNoviceGDMLDetectorConstructionMessenger::
  ~OpNoviceGDMLDetectorConstructionMessenger()
{
  delete fDetConDir;
  delete verboseCmd;
  delete dumpGdmlCmd;
  delete dumpGdmlFileNameCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void OpNoviceGDMLDetectorConstructionMessenger::SetNewValue(
  G4UIcommand* command, G4String newValue)
{
  if(command == verboseCmd)
    fOpNoviceDetCon->SetVerbose(verboseCmd->GetNewBoolValue(newValue));
  if(command == dumpGdmlCmd)
    fOpNoviceDetCon->SetDumpgdml(dumpGdmlCmd->GetNewBoolValue(newValue));
  if(command == dumpGdmlFileNameCmd)
    fOpNoviceDetCon->SetDumpgdmlFile(newValue);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
