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
/// \file errProp/src/ExErrorDetectorMessenger.cc
/// \brief Implementation of the ExErrorDetectorMessenger class
//

#include "ExErrorDetectorMessenger.hh"

#include "ExErrorDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExErrorDetectorMessenger::ExErrorDetectorMessenger(ExErrorDetectorConstruction* myDet)
:G4UImessenger(),fMyDetector(myDet),fMydetDir(0),fFieldCmd(0)
{ 

  fMydetDir = new G4UIdirectory("/exerror/");
  fMydetDir->SetGuidance("ExError control.");
  
  fFieldCmd = new G4UIcmdWithADoubleAndUnit("/exerror/setField",this);  
  fFieldCmd->SetGuidance("Define magnetic field.");
  fFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  fFieldCmd->SetParameterName("Bz",false);
  fFieldCmd->SetDefaultUnit("tesla");
  fFieldCmd->SetUnitCategory("Magnetic flux density");
  fFieldCmd->AvailableForStates(G4State_PreInit,G4State_Init,G4State_Idle);  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExErrorDetectorMessenger::~ExErrorDetectorMessenger()
{
  delete fFieldCmd;
  delete fMydetDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExErrorDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{   
  if( command == fFieldCmd ) { 
    fMyDetector->SetMagField(fFieldCmd->GetNewDoubleValue(newValue));
  }

}

