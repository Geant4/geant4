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
// $Id: G4LogicalVolumeStore.cc,v 1.8 2002-04-26 16:24:36 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4LogicalVolumeStore
//
// Implementation for singleton container
//
// History:
// 10.07.95 P.Kent Initial version
// ********************************************************************

#include "globals.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4GeometryManager.hh"

// ***************************************************************************
// Static class variables
// ***************************************************************************
//
G4LogicalVolumeStore* G4LogicalVolumeStore::fgInstance = 0;
G4bool G4LogicalVolumeStore::locked = false;

// ***************************************************************************
// Protected constructor: Construct underlying container with
// initial size of 100 entries
// ***************************************************************************
//
G4LogicalVolumeStore::G4LogicalVolumeStore()
 : G4std::vector<G4LogicalVolume*>()
{
  reserve(100);
}

// ***************************************************************************
// Destructor
// ***************************************************************************
//
G4LogicalVolumeStore::~G4LogicalVolumeStore()
{
  Clean();
}

// ***************************************************************************
// Delete all elements from the store
// ***************************************************************************
//
void G4LogicalVolumeStore::Clean()
{
  // Do nothing if geometry is closed
  //
  if (G4GeometryManager::GetInstance()->IsGeometryClosed())
  {
    G4cout << "WARNING - Attempt to delete the logical volume store"
           << " while geometry closed !" << G4endl;
    return;
  }

  // Locks store for deletion of volumes. De-registration will be
  // performed at this stage. G4LogicalVolumes will not de-register themselves.
  //
  locked = true;  

  size_t i=0;
  G4LogicalVolumeStore* store = GetInstance();
  G4std::vector<G4LogicalVolume*>::iterator pos;

#ifdef G4GEOMETRY_VOXELDEBUG
  G4cout << "Deleting Logical Volumes ... ";
#endif

  for(pos=store->begin(); pos!=store->end(); pos++)
  {
    if (*pos) delete *pos; i++;
  }

#ifdef G4GEOMETRY_VOXELDEBUG
  if (store->size() < i-1)
    { G4cout << "No volumes deleted. Already deleted by user ?" << G4endl; }
  else
    { G4cout << i-1 << " volumes deleted !" << G4endl; }
#endif

  locked = false;
  store->clear();
}

// ***************************************************************************
// Add volume to container
// ***************************************************************************
//
void G4LogicalVolumeStore::Register(G4LogicalVolume* pVolume)
{
    GetInstance()->push_back(pVolume);
}

// ***************************************************************************
// Remove volume from container
// ***************************************************************************
//
void G4LogicalVolumeStore::DeRegister(G4LogicalVolume* pVolume)
{
  if (!locked)    // Do not de-register if locked !
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

// ***************************************************************************
// Return ptr to Store, setting if necessary
// ***************************************************************************
//
G4LogicalVolumeStore* G4LogicalVolumeStore::GetInstance()
{
    static G4LogicalVolumeStore worldStore;
    if (!fgInstance)
	{
	    fgInstance = &worldStore;
	}
    return fgInstance;
}
