// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T08MagneticField.cc,v 1.1 1999-01-08 16:35:20 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
//   User Field class implementation.
//
#include "T08MagneticField.hh"
#include "G4FieldManager.hh"

//  Constructors:

T08MagneticField::T08MagneticField()
  : G4UniformMagField(G4ThreeVector())
{
  GetGlobalFieldManager()->CreateChordFinder(this);
}

T08MagneticField::T08MagneticField(G4ThreeVector fieldVector)
  : G4UniformMagField(fieldVector)
{    
  GetGlobalFieldManager()->CreateChordFinder(this);
}

// Set the value of the Global Field to fieldValue along Z
//
void T08MagneticField::SetFieldValue(G4double fieldValue)
{
   this->SetFieldValue(G4ThreeVector(0,0,fieldValue));
}

// Set the value of the Global Field
//
void T08MagneticField::SetFieldValue(G4ThreeVector fieldVector)
{
  // Find the Field Manager for the global field
  G4FieldManager* fieldMgr= GetGlobalFieldManager();
    
  if(fieldVector!=G4ThreeVector(0.,0.,0.))
  { 
    this->SetFieldValue(fieldVector);
    fieldMgr->SetDetectorField(this);
  } else {
    // If the new field's value is Zero, then it is best to
    //  insure that it is not used for propagation.
    G4MagneticField* magField = NULL;
    fieldMgr->SetDetectorField(magField);
  }
}

T08MagneticField::~T08MagneticField()
{
  // GetGlobalFieldManager()->SetDetectorField(0);
}

//  Utility method
#include "G4TransportationManager.hh"

G4FieldManager*  T08MagneticField::GetGlobalFieldManager()
{
  return G4TransportationManager::GetTransportationManager()
	     ->GetFieldManager();
}
    
