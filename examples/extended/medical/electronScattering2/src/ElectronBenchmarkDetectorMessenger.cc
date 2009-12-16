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

ElectronBenchmarkDetectorMessenger::ElectronBenchmarkDetectorMessenger(ElectronBenchmarkDetector* det)
:detector(det)
{  
 listDir = new G4UIdirectory("/primFoil/");
 listDir->SetGuidance("Primary foil commands");

 // Foil Material
 primFoilMatCmd = new G4UIcmdWithAString("/primFoil/material",this);
 primFoilMatCmd->SetGuidance("Material of primary foil");
 primFoilMatCmd->SetParameterName("primFoilMat",false);
 primFoilMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

 // Foil Thickness
 primFoilThickCmd = new G4UIcmdWithADoubleAndUnit("/primFoil/thickness",this);
 primFoilThickCmd->SetGuidance("Thickness of primary foil");
 primFoilThickCmd->SetParameterName("thickness",false);
 primFoilThickCmd->SetDefaultUnit("cm");
 primFoilThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

ElectronBenchmarkDetectorMessenger::~ElectronBenchmarkDetectorMessenger()
{
  delete listDir;
  delete primFoilMatCmd;
  delete primFoilThickCmd;
}

void ElectronBenchmarkDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
 if ( command == primFoilMatCmd  ){
   detector->SetPrimFoilMaterial(newValue);
 } else if ( command == primFoilThickCmd  ){
   detector->SetPrimFoilThickness(primFoilThickCmd->GetNewDoubleValue(newValue));
 } else {
   G4cerr << "***** Command is not found !!! " << newValue << G4endl;
 }
}

