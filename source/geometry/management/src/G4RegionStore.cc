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
// $Id: G4RegionStore.cc,v 1.4 2003-02-07 11:46:37 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4RegionStore
//
// Implementation for singleton container
//
// History:
// 18.09.02 G.Cosmo Initial version
// ********************************************************************

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4GeometryManager.hh"

// ***************************************************************************
// Static class variables
// ***************************************************************************
//
G4RegionStore* G4RegionStore::fgInstance = 0;
G4bool G4RegionStore::locked = false;

// ***************************************************************************
// Protected constructor: Construct underlying container with
// initial size of 20 entries
// ***************************************************************************
//
G4RegionStore::G4RegionStore()
  : G4std::vector<G4Region*>()
{
  reserve(20);
}

// ***************************************************************************
// Destructor
// ***************************************************************************
//
G4RegionStore::~G4RegionStore() 
{
  Clean();
}

// ***************************************************************************
// Delete all elements from the store
// ***************************************************************************
//
void G4RegionStore::Clean()
{
  // Do nothing if geometry is closed
  //
  if (G4GeometryManager::GetInstance()->IsGeometryClosed())
  {
    G4cout << "WARNING - Attempt to delete the region store"
           << " while geometry closed !" << G4endl;
    return;
  }

  // Locks store for deletion of regions. De-registration will be
  // performed at this stage. G4Regions will not de-register themselves.
  //
  locked = true;  

  size_t i=0;
  G4RegionStore* store = GetInstance();
  G4std::vector<G4Region*>::iterator pos;

#ifdef G4GEOMETRY_VOXELDEBUG
  G4cout << "Deleting Solids ... ";
#endif

  for(pos=store->begin(); pos!=store->end(); pos++)
  {
    if (*pos) delete *pos; i++;
  }

#ifdef G4GEOMETRY_VOXELDEBUG
  if (store->size() < i-1)
    { G4cout << "No regions deleted. Already deleted by user ?" << G4endl; }
  else
    { G4cout << i-1 << " regions deleted !" << G4endl; }
#endif

  locked = false;
  store->clear();
}

// ***************************************************************************
// Add Solid to container
// ***************************************************************************
//
void G4RegionStore::Register(G4Region* pRegion)
{
  GetInstance()->push_back(pRegion);
}

// ***************************************************************************
// Remove Solid from container
// ***************************************************************************
//
void G4RegionStore::DeRegister(G4Region* pRegion)
{
  if (!locked)    // Do not de-register if locked !
  {
    for (iterator i=GetInstance()->begin(); i!=GetInstance()->end(); i++)
    {
      if (**i==*pRegion)
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
G4RegionStore* G4RegionStore::GetInstance()
{
  static G4RegionStore worldStore;
  if (!fgInstance)
  {
    fgInstance = &worldStore;
  }
  return fgInstance;
}

// ***************************************************************************
// Loops through all regions to verify if a region has been modified.
// It returns TRUE if just one region is modified.
// ***************************************************************************
//
G4bool G4RegionStore::IsModified() const
{
  for (iterator i=GetInstance()->begin(); i!=GetInstance()->end(); i++)
  {
    if ((*i)->IsModified()) return true;
  }
  return false;
}

// ***************************************************************************
// Loops through all regions to reset flag for modification to FALSE.
// Used by the run manager to notify that the physics table has been updated.
// ***************************************************************************
//
void G4RegionStore::ResetRegionModified()
{
  for (iterator i=GetInstance()->begin(); i!=GetInstance()->end(); i++)
  {
    (*i)->RegionModified(false);
  }
}


// ***************************************************************************
// Forces recomputation of material lists in all regions in the store.
// ***************************************************************************
//
void G4RegionStore::UpdateMaterialList()
{
  for (iterator i=GetInstance()->begin(); i!=GetInstance()->end(); i++)
  {
    (*i)->UpdateMaterialList();
  }
}

// ***************************************************************************
// Returns a region through its name specification.
// ***************************************************************************
//
G4Region* G4RegionStore::GetRegion(const G4String& name, G4bool verbose) const
{
  for (iterator i=GetInstance()->begin(); i!=GetInstance()->end(); i++)
  {
    if ((*i)->GetName() == name)
      return *i;
  }
  if (verbose)
  {
    G4cerr << "ERROR - G4RegionStore::GetRegion()" << G4endl
           << "        Region " << name << " NOT found in store !" << G4endl;
  }
  return 0;
}

// ***************************************************************************
// Returns a region through its name specification, if it exists.
// If it does not exist it will allocate a new region with the given
// name, delegating the ownership to the caller client.
// ***************************************************************************
//
G4Region* G4RegionStore::FindOrCreateRegion(const G4String& name)
{
  G4Region* target = GetRegion(name,false);
  if (!target)
  {
    target = new G4Region(name);
  }
  return target;
}
