// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02MagneticField.cc,v 1.4 2000-12-04 16:24:07 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//  
//   User Field class implementation.
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#include "ExN02MagneticField.hh"
#include "G4FieldManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ExN02MagneticField::ExN02MagneticField()
  : G4UniformMagField(G4ThreeVector())
{
  GetGlobalFieldManager()->SetDetectorField(this);
  GetGlobalFieldManager()->CreateChordFinder(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ExN02MagneticField::ExN02MagneticField(G4ThreeVector fieldVector)
  : G4UniformMagField(fieldVector)
{
  GetGlobalFieldManager()->SetDetectorField(this);    
  GetGlobalFieldManager()->CreateChordFinder(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

// Set the value of the Global Field to fieldValue along X
//
void ExN02MagneticField::SetFieldValue(G4double fieldValue)
{
   G4UniformMagField::SetFieldValue(G4ThreeVector(fieldValue,0,0));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

// Set the value of the Global Field
//
void ExN02MagneticField::SetFieldValue(G4ThreeVector fieldVector)
{
  // Find the Field Manager for the global field
  G4FieldManager* fieldMgr= GetGlobalFieldManager();
    
  if(fieldVector!=G4ThreeVector(0.,0.,0.))
  { 
    G4UniformMagField::SetFieldValue(fieldVector);
    fieldMgr->SetDetectorField(this);
  } else {
    // If the new field's value is Zero, then it is best to
    //  insure that it is not used for propagation.
    G4MagneticField* magField = NULL;
    fieldMgr->SetDetectorField(magField);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ExN02MagneticField::~ExN02MagneticField()
{
  // GetGlobalFieldManager()->SetDetectorField(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#include "G4TransportationManager.hh"

G4FieldManager*  ExN02MagneticField::GetGlobalFieldManager()
{
  return G4TransportationManager::GetTransportationManager()->GetFieldManager();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....    
