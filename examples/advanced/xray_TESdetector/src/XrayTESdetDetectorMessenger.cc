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
//
/// \file XrayTESdetDetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
// Authors: P.Dondero (paolo.dondero@cern.ch), R.Stanzani (ronny.stanzani@cern.ch)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "XrayTESdetDetectorMessenger.hh"
#include "XrayTESdetDetectorConstruction.hh"

// #include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
// #include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
// #include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

XrayTESdetDetectorMessenger::XrayTESdetDetectorMessenger(XrayTESdetDetectorConstruction *Det) :
    fpDetector(Det), fTheReadCommand(nullptr)
{
  fTheReadCommand = new G4UIcmdWithAString("/XrayTESdet/det/readGDML", this);
  fTheReadCommand ->SetGuidance("READ GDML file with given name");
  fTheReadCommand ->SetParameterName("FileRead", false);
  fTheReadCommand ->SetDefaultValue("xray_TESdetector.gdml");
  fTheReadCommand ->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

XrayTESdetDetectorMessenger::~XrayTESdetDetectorMessenger()
{
  delete fTheReadCommand;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void XrayTESdetDetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fTheReadCommand) fpDetector->SetReadFile(newValue);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
