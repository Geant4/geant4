// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicalVolumeStore.cc,v 1.3 1999-11-11 15:35:53 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4PhysicalVolumeStore
//
// Implementation for singleton container
//
// History:
// 25.07.95 P.Kent Initial version

#include "G4PhysicalVolumeStore.hh"
#include "globals.hh"

// Protected constructor: Construct underlying container with
// initial size of 100 entries
G4PhysicalVolumeStore::G4PhysicalVolumeStore() : G4RWTPtrOrderedVector<G4VPhysicalVolume>(100)
{
}

// Destructor
G4PhysicalVolumeStore::~G4PhysicalVolumeStore()
{
  while (!isEmpty()) removeFirst();
}

// Static class variable
G4PhysicalVolumeStore* G4PhysicalVolumeStore::fgInstance = 0;

// Add Solid to container
void G4PhysicalVolumeStore::Register(G4VPhysicalVolume* pVolume)
{
    GetInstance()->insert(pVolume);
}

// Remove Solid from container
void G4PhysicalVolumeStore::DeRegister(G4VPhysicalVolume* pVolume)
{
    GetInstance()->remove(pVolume);
}

// Return ptr to Store, setting if necessary
G4PhysicalVolumeStore* G4PhysicalVolumeStore::GetInstance()
{
    static G4PhysicalVolumeStore worldStore;
    if (!fgInstance)
	{
	    fgInstance = &worldStore;
	}
    return fgInstance;
}


