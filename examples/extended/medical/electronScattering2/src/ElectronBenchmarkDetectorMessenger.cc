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
/// \file medical/electronScattering2/src/ElectronBenchmarkDetectorMessenger.cc
/// \brief Implementation of the ElectronBenchmarkDetectorMessenger class

#include "ElectronBenchmarkDetectorMessenger.hh"
#include "ElectronBenchmarkDetector.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ElectronBenchmarkDetectorMessenger::ElectronBenchmarkDetectorMessenger(
                                ElectronBenchmarkDetector* det)
: G4UImessenger(), fDetector(det),
  fListDir(0), fPrimFoilMatCmd(0), fPrimFoilThickCmd(0)
{
    fListDir = new G4UIdirectory("/primFoil/");
    fListDir->SetGuidance("Primary foil commands");
    
    // Foil Material
    fPrimFoilMatCmd = new G4UIcmdWithAString("/primFoil/material",this);
    fPrimFoilMatCmd->SetGuidance("Material of primary foil");
    fPrimFoilMatCmd->SetParameterName("primFoilMat",false);
    fPrimFoilMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    fPrimFoilMatCmd->SetToBeBroadcasted(false);
    
    // Foil Thickness
    fPrimFoilThickCmd=new G4UIcmdWithADoubleAndUnit("/primFoil/thickness",this);
    fPrimFoilThickCmd->SetGuidance("Thickness of primary foil");
    fPrimFoilThickCmd->SetParameterName("thickness",false);
    fPrimFoilThickCmd->SetDefaultUnit("cm");
    fPrimFoilThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    fPrimFoilThickCmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ElectronBenchmarkDetectorMessenger::~ElectronBenchmarkDetectorMessenger()
{
    delete fListDir;
    delete fPrimFoilMatCmd;
    delete fPrimFoilThickCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ElectronBenchmarkDetectorMessenger::SetNewValue(G4UIcommand* command,
                                                     G4String newValue)
{
    if ( command == fPrimFoilMatCmd  ){
        fDetector->SetPrimFoilMaterial(newValue);
    } else if ( command == fPrimFoilThickCmd  ){
        fDetector->SetPrimFoilThickness(
                                fPrimFoilThickCmd->GetNewDoubleValue(newValue));
    } else {
        G4cerr << "***** Command is not found !!! " << newValue << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
