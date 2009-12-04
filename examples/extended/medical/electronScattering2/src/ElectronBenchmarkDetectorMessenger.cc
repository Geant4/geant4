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

