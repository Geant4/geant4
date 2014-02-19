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
// $Id: G4GlobalMagFieldMessenger.cc 66536 2012-12-19 14:32:36Z ihrivnac $
// 
// class G4GlobalMagFieldMessenger
//
// Implementation
//
// Implementation of the G4GlobalMagFieldMessenger class
//
// Author: Ivana Hrivnacova, 28/08/2013  (ivana@ipno.in2p3.fr)
// --------------------------------------------------------------------

#include "G4GlobalMagFieldMessenger.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//______________________________________________________________________________

G4GlobalMagFieldMessenger::G4GlobalMagFieldMessenger(const G4ThreeVector& value)
 : G4UImessenger(),
   fMagField(0),
   fVerboseLevel(0),
   fDirectory(0),
   fSetValueCmd(0),
   fSetVerboseCmd(0)
{
  fDirectory = new G4UIdirectory("/globalField/");
  fDirectory->SetGuidance("Global uniform magnetic field UI commands");
  
  fSetValueCmd = new G4UIcmdWith3VectorAndUnit("/globalField/setValue",this);
  fSetValueCmd->SetGuidance("Set uniform magnetic field value.");
  fSetValueCmd->SetParameterName("Bx", "By", "By", false);
  fSetValueCmd->SetUnitCategory("Magnetic flux density");
  fSetValueCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fSetVerboseCmd = new G4UIcmdWithAnInteger("/globalField/verbose",this);
  fSetVerboseCmd->SetGuidance("Set verbose level: ");
  fSetVerboseCmd->SetGuidance("  0: no output");
  fSetVerboseCmd->SetGuidance("  1: printing new field value");
  fSetVerboseCmd->SetParameterName("globalFieldVerbose", false);
  fSetVerboseCmd->SetRange("globalFieldVerbose>=0");
  fSetVerboseCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  // Create field
  fMagField = new G4UniformMagField(value);
  
  // Set field value (the field is not activated if value is zero)
  SetField(value, "G4GlobalMagFieldMessenger::G4GlobalMagFieldMessenger");
}

//______________________________________________________________________________

G4GlobalMagFieldMessenger::~G4GlobalMagFieldMessenger()
{
  delete fMagField;
  delete fSetValueCmd;
  delete fSetVerboseCmd;
  delete fDirectory;
}

//______________________________________________________________________________

void G4GlobalMagFieldMessenger::SetField(const G4ThreeVector& value,
                                         const G4String& /*inFunction*/)
{
  // Get field manager (or create it if it does not yet exist)
  G4FieldManager* fieldManager
    = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  // Inactivate field if its value is zero
  if ( value == G4ThreeVector() )
  {
    fieldManager->SetDetectorField(0);
    fieldManager->CreateChordFinder(0);
    
    if ( fVerboseLevel > 0 )
    {
      G4cout << "Magnetic field is inactive, fieldValue = (0,0,0)." << G4endl;
    }
  } 
  else
  { 
    fMagField->SetFieldValue(value);
    fieldManager->SetDetectorField(fMagField);
    fieldManager->CreateChordFinder(fMagField);
    
    if ( fVerboseLevel > 0 )
    {
      G4cout << "Magnetic field is active, fieldValue = (" 
             << G4BestUnit(value, "Magnetic flux density") 
             << ")." << G4endl;
    }
  }  
}

//______________________________________________________________________________

void G4GlobalMagFieldMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if ( command == fSetValueCmd )
  {
    SetField(fSetValueCmd->GetNew3VectorValue(newValue),
             "G4GlobalMagFieldMessenger::SetNewValue");
  }
  else if ( command == fSetVerboseCmd )
  {
    SetVerboseLevel(fSetVerboseCmd->GetNewIntValue(newValue));
  }  
}

//______________________________________________________________________________

void G4GlobalMagFieldMessenger::SetFieldValue(const G4ThreeVector& value)
{
  SetField(value, "G4GlobalMagFieldMessenger::SetFieldValue");
}  

//______________________________________________________________________________

G4ThreeVector G4GlobalMagFieldMessenger::GetFieldValue() const
{
  if ( fMagField ) return fMagField->GetConstantFieldValue();
  
  return G4ThreeVector();
}  
