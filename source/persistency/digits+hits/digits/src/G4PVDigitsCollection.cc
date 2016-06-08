// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PVDigitsCollection.cc,v 1.3 1999/11/29 20:11:17 morita Exp $
// GEANT4 tag $Name: geant4-02-00 $
//

// G4PVDigitsCollection

#include <assert.h>
#include "G4PVDigitsCollection.hh"
#include "G4PersistentDigitMan.hh"
#include "G4PDCofThisEvent.hh"

G4PVDigitsCollection::G4PVDigitsCollection(G4String detName,G4String colNam)
 : pcollectionName(colNam), pDMname(detName)
{;}

G4PVDigitsCollection::~G4PVDigitsCollection()
{;}

int G4PVDigitsCollection::operator==(const G4PVDigitsCollection &right) const
{ 
  return ((pcollectionName==right.pcollectionName)
        &&(pDMname==right.pDMname));
}

