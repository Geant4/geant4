// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ElectronOccupancy.cc,v 1.1 1999-08-18 09:15:23 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      Hisaya Kurashige, 17 Aug 1999
// ----------------------------------------------------------------
//     This class has information of occupation of electrons 
//     in atomic orbits
// ---------------------------------------------------------------


#include "G4ElectronOccupancy.hh"

G4ElectronOccupancy::G4ElectronOccupancy(G4int sizeOrbit ):
                  theSizeOfOrbit(sizeOrbit)
{
  // check size
  if  ( (theSizeOfOrbit <1 ) || (theSizeOfOrbit > MaxSizeOfOrbit) ) {
    theSizeOfOrbit = MaxSizeOfOrbit;
  }

  // allocate and clear the array of theOccupancies 
  G4int*  theOccupancies = new G4int[theSizeOfOrbit];
  G4int   index =0;
  for (index = 0; index <  theSizeOfOrbit; index++) {
    theOccupancies[index] =0;
  }

   theTotalOccupancy =0;
}

G4ElectronOccupancy::~G4ElectronOccupancy()
{
   theSizeOfOrbit = -1;

   delete [] theOccupancies;
   theOccupancies =0;
   theTotalOccupancy =0;
}


G4ElectronOccupancy::G4ElectronOccupancy(const G4ElectronOccupancy& right)
{
  theSizeOfOrbit = right.theSizeOfOrbit;

  // allocate and clear the array of theOccupancies 
  G4int*  theOccupancies = new G4int[theSizeOfOrbit];
  G4int   index =0;
  for (index = 0; index <  theSizeOfOrbit; index++) {
    theOccupancies[index] = right.theOccupancies[index];
  }

  theTotalOccupancy = right.theTotalOccupancy;
}

G4ElectronOccupancy& G4ElectronOccupancy::operator=(const G4ElectronOccupancy& right)
{
  if ( this != &right) {
    theSizeOfOrbit = right.theSizeOfOrbit;
    
    // allocate and clear the array of theOccupancies 
    G4int*  theOccupancies = new G4int[theSizeOfOrbit];
    G4int   index =0;
    for (index = 0; index <  theSizeOfOrbit; index++) {
      theOccupancies[index] = right.theOccupancies[index];
    }
    
    theTotalOccupancy = right.theTotalOccupancy;
  }
  return *this;
}

G4int G4ElectronOccupancy::operator==(const G4ElectronOccupancy& right) const
{
  G4int index;
  G4bool value = true;
  for (index = 0; index < MaxSizeOfOrbit; index++) {
    if ( (index < theSizeOfOrbit ) && ( index < right.theSizeOfOrbit) ) {
      value = value && 
         (theOccupancies[index] == right.theOccupancies[index]) ;
    } else if ((index < theSizeOfOrbit ) && ( index >= right.theSizeOfOrbit)) {
      value = value && (theOccupancies[index] == 0);
    } else if ((index >= theSizeOfOrbit ) && ( index <right.theSizeOfOrbit)) {
      value = value && (right.theOccupancies[index] == 0);
    }
  }
  return value;
}

G4int G4ElectronOccupancy::operator!=(const G4ElectronOccupancy& right) const
{
  return !(*this == right);
}


void G4ElectronOccupancy::DumpInfo() const
{
  G4cout << "  -- Electron Occupancy -- " << endl;
  G4int index;
  for (index = 0; index < theSizeOfOrbit; index++) {
    G4cout << "   " << index << "-th orbit       " 
           <<  theOccupancies[index] << endl;
  }
}
