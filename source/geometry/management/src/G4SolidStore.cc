// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SolidStore.cc,v 1.1 1999-01-07 16:07:21 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4SolidStore
//
// Implementation for singleton container
//
// History:
// 10.07.95 P.Kent Initial version

#include "G4SolidStore.hh"
#include "globals.hh"

// Protected constructor: Construct underlying container with
// initial size of 100 entries
G4SolidStore::G4SolidStore() : RWTPtrOrderedVector<G4VSolid>(100)
{
}

// Destructor
G4SolidStore::~G4SolidStore() 
{
  while (!isEmpty()) delete first();
}

// Static class variable
G4SolidStore* G4SolidStore::fgInstance = 0;

// Add Solid to container
void G4SolidStore::Register(G4VSolid* pSolid)
{
    GetInstance()->insert(pSolid);
}

// Remove Solid from container
void G4SolidStore::DeRegister(G4VSolid* pSolid)
{
    GetInstance()->remove(pSolid);
}

// Return ptr to Store, setting if necessary
G4SolidStore* G4SolidStore::GetInstance()
{
    static G4SolidStore worldStore;
    if (!fgInstance)
	{
	    fgInstance = &worldStore;
	}
    return fgInstance;
}
