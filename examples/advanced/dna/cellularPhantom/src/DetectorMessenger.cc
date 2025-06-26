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
// -----------------------------------------------------------------------------
//       MONTE CARLO SIMULATION OF REALISTIC GEOMETRY FROM MICROSCOPES IMAGES
//
// Authors and contributors:
// P. Barberet (a), S. Incerti (a), N. H. Tran (a), L. Morelli (a,b)
//
// a) University of Bordeaux, CNRS, LP2i, UMR5797, Gradignan, France
// b) Politecnico di Milano, Italy
//
// If you use this code, please cite the following publication:
// P. Barberet et al.,
// "Monte-Carlo dosimetry on a realistic cell monolayer
// geometry exposed to alpha particles."
// Ph. Barberet et al 2012 Phys. Med. Biol. 57 2189
// doi: 110.1088/0031-9155/57/8/2189
// -----------------------------------------------------------------------------

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * det)
:G4UImessenger(), fDetector(det)
{
  fPhantomDir = new G4UIdirectory("/phantom/");
  fPhantomDir->SetGuidance(" Cell phantom settings");

  fNameCmd = new G4UIcmdWithAString("/phantom/fileName",this);
  fNameCmd->SetGuidance("Select phantom file name");
  fNameCmd->SetParameterName("fileName",true);
  fNameCmd->SetDefaultValue("phantom.dat");
  fNameCmd->AvailableForStates(G4State_PreInit);

  fMatCmd = new G4UIcmdWithAString("/phantom/mediumMat",this);
  fMatCmd->SetGuidance("Select material for the phantom medium");
  fMatCmd->SetParameterName("mediumMat",true);
  fMatCmd->AvailableForStates(G4State_PreInit);

  fDenRedCmd = new G4UIcmdWithADoubleAndUnit("/phantom/redDen",this);
  fDenRedCmd->SetGuidance("Select density for the red volume");
  fDenRedCmd->SetParameterName("redDen",true);
  fDenRedCmd->SetDefaultValue(1.);
  fDenRedCmd->SetDefaultUnit("g/cm3");
  fDenRedCmd->AvailableForStates(G4State_PreInit);

  fDenGreenCmd = new G4UIcmdWithADoubleAndUnit("/phantom/greenDen",this);
  fDenGreenCmd->SetGuidance("Select density for the green volume");
  fDenGreenCmd->SetParameterName("greenDen",true);
  fDenGreenCmd->SetDefaultValue(1.);
  fDenGreenCmd->SetDefaultUnit("g/cm3");
  fDenGreenCmd->AvailableForStates(G4State_PreInit);

  fDenBlueCmd = new G4UIcmdWithADoubleAndUnit("/phantom/blueDen",this);
  fDenBlueCmd->SetGuidance("Select density for the blue volume");
  fDenBlueCmd->SetParameterName("blueDen",true);
  fDenBlueCmd->SetDefaultValue(1.);
  fDenBlueCmd->SetDefaultUnit("g/cm3");
  fDenBlueCmd->AvailableForStates(G4State_PreInit);

  fShiftXCmd = new G4UIcmdWithADoubleAndUnit("/phantom/shiftX",this);
  fShiftXCmd->SetGuidance("Set phantom X shift");
  fShiftXCmd->SetParameterName("shiftX",true);
  fShiftXCmd->SetDefaultValue(0.);
  fShiftXCmd->SetDefaultUnit("um");
  fShiftXCmd->AvailableForStates(G4State_PreInit);

  fShiftYCmd = new G4UIcmdWithADoubleAndUnit("/phantom/shiftY",this);
  fShiftYCmd->SetGuidance("Set phantom Y shift");
  fShiftYCmd->SetParameterName("shiftY",true);
  fShiftYCmd->SetDefaultValue(0.);
  fShiftYCmd->SetDefaultUnit("um");
  fShiftYCmd->AvailableForStates(G4State_PreInit);

  fShiftZCmd = new G4UIcmdWithADoubleAndUnit("/phantom/shiftZ",this);
  fShiftZCmd->SetGuidance("Set phantom Z shift");
  fShiftZCmd->SetParameterName("shiftZ",true);
  fShiftZCmd->SetDefaultValue(0.);
  fShiftZCmd->SetDefaultUnit("um");
  fShiftZCmd->AvailableForStates(G4State_PreInit);

  fMediumSizeXYCmd = new G4UIcmdWithADoubleAndUnit("/phantom/mediumSizeXY",this);
  fMediumSizeXYCmd->SetGuidance("Set cellular medium size XY");
  fMediumSizeXYCmd->SetParameterName("mediumSizeXY",false);
  fMediumSizeXYCmd->SetDefaultUnit("um");
  fMediumSizeXYCmd->AvailableForStates(G4State_PreInit);

  fMediumSizeZCmd = new G4UIcmdWithADoubleAndUnit("/phantom/mediumSizeZ",this);
  fMediumSizeZCmd->SetGuidance("Set cellular medium size Z");
  fMediumSizeZCmd->SetParameterName("mediumSizeZ",false);
  fMediumSizeZCmd->SetDefaultUnit("um");
  fMediumSizeZCmd->AvailableForStates(G4State_PreInit);

  fWorldDir = new G4UIdirectory("/world/");
  fWorldDir->SetGuidance(" World volume settings");

  fWorldSizeXYCmd = new G4UIcmdWithADoubleAndUnit("/world/sizeXY",this);
  fWorldSizeXYCmd->SetGuidance("Set world size XY");
  fWorldSizeXYCmd->SetParameterName("sizeXY",false);
  fWorldSizeXYCmd->SetDefaultUnit("um");
  fWorldSizeXYCmd->AvailableForStates(G4State_PreInit);

  fWorldSizeZCmd = new G4UIcmdWithADoubleAndUnit("/world/sizeZ",this);
  fWorldSizeZCmd->SetGuidance("Set world size Z");
  fWorldSizeZCmd->SetParameterName("sizeZ",false);
  fWorldSizeZCmd->SetDefaultUnit("um");
  fWorldSizeZCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fWorldDir;
  delete fPhantomDir;
  delete fNameCmd;
  delete fMatCmd;
  delete fDenRedCmd;
  delete fDenGreenCmd;
  delete fDenBlueCmd;
  delete fShiftXCmd;
  delete fShiftYCmd;
  delete fShiftZCmd;
  delete fMediumSizeXYCmd;
  delete fMediumSizeZCmd;
  delete fWorldSizeXYCmd;
  delete fWorldSizeZCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if( command == fMatCmd ) {
    fDetector->SetTargetMaterial(newValue);
  }
  else if(command == fDenRedCmd) {
    fDetector->SetRedDensity(fDenRedCmd->GetNewDoubleValue(newValue));
  }
  else if(command == fDenGreenCmd) {
    fDetector->SetGreenDensity(fDenGreenCmd->GetNewDoubleValue(newValue));
  }
  else if(command == fDenBlueCmd) {
    fDetector->SetBlueDensity(fDenBlueCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fShiftXCmd) {
    fDetector->SetShiftX(fShiftXCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fShiftYCmd) {
    fDetector->SetShiftY(fShiftYCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fShiftZCmd) {
    fDetector->SetShiftZ(fShiftZCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fMediumSizeXYCmd) {
    fDetector->SetMediumSizeXY(fMediumSizeXYCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fMediumSizeZCmd) {
    fDetector->SetMediumSizeZ(fMediumSizeZCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fWorldSizeXYCmd) {
    fDetector->SetWorldSizeXY(fWorldSizeXYCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fWorldSizeZCmd) {
    fDetector->SetWorldSizeZ(fWorldSizeZCmd->GetNewDoubleValue(newValue));
  }
  else if(command == fNameCmd) {
    fDetector->SetPhantomFileName(newValue);
  }
}
