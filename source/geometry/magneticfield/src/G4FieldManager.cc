// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FieldManager.cc,v 1.3 2000-11-09 18:06:37 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
   if ( fAllocatedChordFinder )
      delete fChordFinder;
   fChordFinder= new G4ChordFinder( detectorMagField );
   fAllocatedChordFinder= true;
}

