//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4PhysicalVolumeStore.cc,v 1.8 2002-04-19 08:20:22 gcosmo Exp $
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

// Static class variables
G4PhysicalVolumeStore* G4PhysicalVolumeStore::fgInstance = 0;
G4bool G4PhysicalVolumeStore::locked = false;

// Protected constructor: Construct underlying container with
// initial size of 100 entries
G4PhysicalVolumeStore::G4PhysicalVolumeStore()
  : G4std::vector<G4VPhysicalVolume*>()
{
  reserve(100);
}

// Destructor
G4PhysicalVolumeStore::~G4PhysicalVolumeStore()
{
  Clean();
}

// Delete all elements from the store
void G4PhysicalVolumeStore::Clean()
{
  size_t i=0;
  locked = true;
  G4PhysicalVolumeStore* store = GetInstance();
  G4std::vector<G4VPhysicalVolume*>::iterator pos;

#ifdef G4GEOMETRY_VOXELDEBUG
  G4cout << "Deleting Physical Volumes ... ";
#endif

  for(pos=store->begin(); pos!=store->end(); pos++)
  {
    if (*pos) delete *pos; i++;
  }

#ifdef G4GEOMETRY_VOXELDEBUG
  G4cout << i-1 << " volumes deleted !" << G4endl;
#endif

  locked = false;
  store->clear();
}

// Add Solid to container
void G4PhysicalVolumeStore::Register(G4VPhysicalVolume* pVolume)
{
    GetInstance()->push_back(pVolume);
}

// Remove Solid from container
void G4PhysicalVolumeStore::DeRegister(G4VPhysicalVolume* pVolume)
{
  if (!locked)
  {
    for (iterator i=GetInstance()->begin(); i!=GetInstance()->end(); i++)
    {
      if (**i==*pVolume)
      {
        GetInstance()->erase(i);
        break;
      }
    }
  }
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
