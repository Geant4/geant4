// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02DetectorMessenger.cc,v 1.3 2000-12-04 16:24:07 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#include "ExN02DetectorMessenger.hh"

#include "ExN02DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ExN02DetectorMessenger::ExN02DetectorMessenger(ExN02DetectorConstruction * myDet)
:myDetector(myDet)
{ 

  mydetDir = new G4UIdirectory("/mydet/");
  mydetDir->SetGuidance("ExN02 detector control.");
  
  TargMatCmd = new G4UIcmdWithAString("/mydet/setTargetMate",this);
  TargMatCmd->SetGuidance("Select Material of the Target.");
  TargMatCmd->SetParameterName("choice",false);
  TargMatCmd->AvailableForStates(PreInit,Idle);
  
  ChamMatCmd = new G4UIcmdWithAString("/mydet/setChamberMate",this);
  ChamMatCmd->SetGuidance("Select Material of the Target.");
  ChamMatCmd->SetParameterName("choice",false);
  ChamMatCmd->AvailableForStates(PreInit,Idle);  
  
  FieldCmd = new G4UIcmdWithADoubleAndUnit("/mydet/setField",this);  
  FieldCmd->SetGuidance("Define magnetic field.");
  FieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  FieldCmd->SetParameterName("Bz",false);
  FieldCmd->SetDefaultUnit("tesla");
  FieldCmd->SetUnitCategory("Magnetic flux density");
  FieldCmd->AvailableForStates(PreInit,Idle);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ExN02DetectorMessenger::~ExN02DetectorMessenger()
{
  delete TargMatCmd;
  delete ChamMatCmd;
  delete FieldCmd;
  delete mydetDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ExN02DetectorMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{ 
  if( command == TargMatCmd )
   { myDetector->setTargetMaterial(newValues);}
   
  if( command == ChamMatCmd )
   { myDetector->setChamberMaterial(newValues);}  
  
  if( command == FieldCmd )
   { myDetector->SetMagField(FieldCmd->GetNewDoubleValue(newValues));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
