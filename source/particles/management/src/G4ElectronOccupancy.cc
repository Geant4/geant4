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
// $Id: G4ElectronOccupancy.cc 79357 2014-02-25 10:06:54Z gcosmo $
//
// 
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, based on object model of
//      Hisaya Kurashige, 17 Aug 1999
// ----------------------------------------------------------------
//     This class has information of occupation of electrons 
//     in atomic orbits
// ---------------------------------------------------------------

#include "G4ElectronOccupancy.hh"
#include <sstream>

G4ThreadLocal G4Allocator<G4ElectronOccupancy> *aElectronOccupancyAllocator = 0;

G4ElectronOccupancy::G4ElectronOccupancy(G4int sizeOrbit )
  : theSizeOfOrbit(sizeOrbit)
{
  // check size
  if  ( (theSizeOfOrbit <1 ) || (theSizeOfOrbit > MaxSizeOfOrbit) ) {
    theSizeOfOrbit = MaxSizeOfOrbit;
  }

  // allocate and clear the array of theOccupancies 
  theOccupancies = new G4int[theSizeOfOrbit];
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
  theOccupancies = new G4int[theSizeOfOrbit];
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
    if ( theOccupancies != 0 ) delete [] theOccupancies;
    theOccupancies = new G4int[theSizeOfOrbit];
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
  G4cout << "  -- Electron Occupancy -- " << G4endl;
  G4int index;
  for (index = 0; index < theSizeOfOrbit; index++) {
    G4cout << "   " << index << "-th orbit       " 
           <<  theOccupancies[index] << G4endl;
  }
}

G4int  G4ElectronOccupancy::AddElectron(G4int orbit, G4int number)
{
  G4int value =0;
  if (orbit>=theSizeOfOrbit){
    std::ostringstream smsg;
    smsg<<  "Orbit (" << orbit 
	<<") exceeds the maximum("
	<<theSizeOfOrbit-1<<")  ";
    G4String msg = smsg.str();
    G4Exception("G4ElectronOccupancy::AddElectron()","PART131",
		JustWarning, msg);
  } else if (orbit >=0) {
    theOccupancies[orbit] += number;
    theTotalOccupancy += number; 
    value = number;   
  }
  return value;
}

G4int  G4ElectronOccupancy::RemoveElectron(G4int orbit, G4int number)
{
  G4int value =0;
  if (orbit>=theSizeOfOrbit){
    std::ostringstream smsg;
    smsg<<  "Orbit (" << orbit 
	<<") exceeds the maximum("
	<<theSizeOfOrbit-1 <<") ";
    G4String msg = smsg.str();
    G4Exception("G4ElectronOccupancy::RemoveElectron()","PART131",
		JustWarning, msg);
  } else if (orbit >=0) {
    if ( theOccupancies[orbit] < number ) number = theOccupancies[orbit];
    theOccupancies[orbit] -= number;
    theTotalOccupancy -= number;    
    value = number;
  }
  return value;
}
