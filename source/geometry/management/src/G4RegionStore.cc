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
// $Id: G4RegionStore.cc 103235 2017-03-22 15:53:48Z gcosmo $
//
// G4RegionStore
//
// Implementation for singleton container
//
// History:
// 18.09.02 G.Cosmo Initial version
// --------------------------------------------------------------------

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4GeometryManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4ios.hh"

// ***************************************************************************
// Static class variables
// ***************************************************************************
//
G4RegionStore* G4RegionStore::fgInstance = 0;
G4ThreadLocal G4VStoreNotifier* G4RegionStore::fgNotifier = 0;
G4ThreadLocal G4bool G4RegionStore::locked = false;

// ***************************************************************************
// Protected constructor: Construct underlying container with
// initial size of 20 entries
// ***************************************************************************
//
G4RegionStore::G4RegionStore()
  : std::vector<G4Region*>()
{
  reserve(20);
}

// ***************************************************************************
// Destructor
// ***************************************************************************
//
G4RegionStore::~G4RegionStore() 
{
  Clean();  // Delete all regions in the store
  G4Region::Clean();  // Delete allocated sub-instance data
}

// ***************************************************************************
// Delete all regions from the store except for the world region
// ***************************************************************************
//
void G4RegionStore::Clean()
{
  // Do nothing if geometry is closed
  //
  if (G4GeometryManager::IsGeometryClosed())
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

#ifdef G4GEOMETRY_VOXELDEBUG
  G4cout << "Deleting Regions ... ";
#endif

  for(iterator pos=store->begin(); pos!=store->end(); ++pos)
  {
    if (fgNotifier) { fgNotifier->NotifyDeRegistration(); }
    if (*pos) { delete *pos; }
    i++;
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
// Associate user notifier to the store
// ***************************************************************************
//
void G4RegionStore::SetNotifier(G4VStoreNotifier* pNotifier)
{
  GetInstance();
  fgNotifier = pNotifier;
}

// ***************************************************************************
// Add Region to container
// ***************************************************************************
//
void G4RegionStore::Register(G4Region* pRegion)
{
  GetInstance()->push_back(pRegion);
  if (fgNotifier)  { fgNotifier->NotifyRegistration(); }
}

// ***************************************************************************
// Remove Region from container
// ***************************************************************************
//
void G4RegionStore::DeRegister(G4Region* pRegion)
{
  if (!locked)    // Do not de-register if locked !
  {
    if (fgNotifier)  { fgNotifier->NotifyDeRegistration(); }
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
    if ((*i)->IsModified()) { return true; }
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
void G4RegionStore::UpdateMaterialList(G4VPhysicalVolume* currentWorld)
{
  for (iterator i=GetInstance()->begin(); i!=GetInstance()->end(); i++)
  {
    if((*i)->IsInMassGeometry() || (*i)->IsInParallelGeometry() || currentWorld)
    { (*i)->UpdateMaterialList(); }
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
    if ((*i)->GetName() == name) { return *i; }
  }
  if (verbose)
  {
    std::ostringstream message;
    message << "Region NOT found in store !" << G4endl
            << "        Region " << name << " NOT found in store !" << G4endl
            << "        Returning NULL pointer.";
    G4Exception("G4RegionStore::GetRegion()",
                "GeomMgt1001", JustWarning, message);
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

// **************************************************************************
// Set a world physical volume pointer to a region that belongs to it.
// Scan over all world volumes.
// **************************************************************************
//
void G4RegionStore::SetWorldVolume()
{
  // Reset all pointers first
  //
  for (iterator i=GetInstance()->begin(); i!=GetInstance()->end(); i++)
  { (*i)->SetWorld(0); }

  // Find world volumes
  //
  G4PhysicalVolumeStore* fPhysicalVolumeStore
    = G4PhysicalVolumeStore::GetInstance();
  size_t nPhys = fPhysicalVolumeStore->size();
  for(size_t iPhys=0; iPhys<nPhys; iPhys++)
  {
    G4VPhysicalVolume* fPhys = (*fPhysicalVolumeStore)[iPhys];
    if(fPhys->GetMotherLogical()) { continue; } // not a world volume

    // Now 'fPhys' is a world volume, set it to regions that belong to it.
    //
    for (iterator i=GetInstance()->begin(); i!=GetInstance()->end(); i++)
    { (*i)->SetWorld(fPhys); }
  }
}

