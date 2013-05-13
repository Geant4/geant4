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
// $Id: B2FieldMessenger.cc 66536 2012-12-19 14:32:36Z ihrivnac $
// 
/// \file B2FieldMessenger.cc
/// \brief Implementation of the B2FieldMessenger class

#include "B2FieldMessenger.hh"
#include "B2MagneticField.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2FieldMessenger::B2FieldMessenger(B2MagneticField* field)
 : G4UImessenger(),
   fField(field)
{
  fSetFieldCmd = new G4UIcmdWithADoubleAndUnit("/B2/det/setField",this);
  fSetFieldCmd->SetGuidance("Define magnetic field.");
  fSetFieldCmd->SetGuidance("Magnetic field will be in X direction.");
  fSetFieldCmd->SetParameterName("Bx",false);
  fSetFieldCmd->SetUnitCategory("Magnetic flux density");
  fSetFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2FieldMessenger::~B2FieldMessenger()
{
  delete fSetFieldCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2FieldMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == fSetFieldCmd ) {
    fField
      ->SetMagFieldValue(fSetFieldCmd->GetNewDoubleValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
