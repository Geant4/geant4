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
/// \file electromagnetic/TestEm1/src/FieldMessenger.cc
/// \brief Implementation of the FieldMessenger class
//
// $Id: FieldMessenger.cc $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "FieldMessenger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
FieldMessenger::FieldMessenger(DetectorConstruction* Det)
  :G4UImessenger(),
   fMagFieldCmd(0),
   fDetector(Det)
{
  fMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setField",this);
  fMagFieldCmd->SetGuidance("Define magnetic field.");
  fMagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  fMagFieldCmd->SetParameterName("Bz",false);
  fMagFieldCmd->SetUnitCategory("Magnetic flux density");
    //fMagFieldCmd->SetDefaultUnit("tesla");
  fMagFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
} 

FieldMessenger::~FieldMessenger()
{
  delete fMagFieldCmd;
}

void FieldMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if ( command == fMagFieldCmd )
    { fDetector->SetMagField(fMagFieldCmd->GetNewDoubleValue(newValue)); }
}
