// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PVHitsCollection.cc,v 1.8 1999/11/28 21:54:16 morita Exp $
// GEANT4 tag $Name: geant4-03-01 $
//

// G4PVHitsCollection

#include <assert.h>
#include "G4PVHitsCollection.hh"
#include "G4PersistentHitMan.hh"
#include "G4PHCofThisEvent.hh"

G4PVHitsCollection::G4PVHitsCollection(G4String detName,G4String colNam)
 : pcollectionName(colNam), pSDname(detName)
{;}

G4PVHitsCollection::~G4PVHitsCollection()
{;}

int G4PVHitsCollection::operator==(const G4PVHitsCollection &right) const
{ 
  return ((pcollectionName==right.pcollectionName)
        &&(pSDname==right.pSDname));
}

