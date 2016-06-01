// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FieldManager.cc,v 2.1 1998/07/12 02:55:15 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//
#include "G4FieldManager.hh"

G4FieldManager::G4FieldManager()
{ 
   fDetectorField= 0;
   fAllocatedChordFinder= false;
}

G4FieldManager::G4FieldManager(G4MagneticField *detectorField)
{
   fDetectorField= detectorField;

   this->CreateChordFinder(detectorField);
}


G4FieldManager::~G4FieldManager()
{
   if( fAllocatedChordFinder ){
      delete fChordFinder;
   }
}

void
G4FieldManager::CreateChordFinder(G4MagneticField *detectorMagField)
{
   fChordFinder= new G4ChordFinder( detectorMagField );
   fAllocatedChordFinder= true;
}

