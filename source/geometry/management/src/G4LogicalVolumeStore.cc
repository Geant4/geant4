// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LogicalVolumeStore.cc,v 1.3 1999-11-11 15:35:53 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4LogicalVolumeStore
//
// Implementation for singleton container
//
// History:
// 10.07.95 P.Kent Initial version

#include "G4LogicalVolumeStore.hh"
#include "globals.hh"

// Protected constructor: Construct underlying container with
// initial size of 100 entries
G4LogicalVolumeStore::G4LogicalVolumeStore() : G4RWTPtrOrderedVector<G4LogicalVolume>(100)
{
}

// Destructor
G4LogicalVolumeStore::~G4LogicalVolumeStore()
{
 while (!isEmpty()) removeFirst();
}

// Static class variable
G4LogicalVolumeStore* G4LogicalVolumeStore::fgInstance = 0;

// Add volume to container
void G4LogicalVolumeStore::Register(G4LogicalVolume* pVolume)
{
    GetInstance()->insert(pVolume);
}

// Remove volume from container
void G4LogicalVolumeStore::DeRegister(G4LogicalVolume* pVolume)
{
    GetInstance()->remove(pVolume);
}

// Return ptr to Store, setting if necessary
G4LogicalVolumeStore* G4LogicalVolumeStore::GetInstance()
{
    static G4LogicalVolumeStore worldStore;
    if (!fgInstance)
	{
	    fgInstance = &worldStore;
	}
    return fgInstance;
}
