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
// G4RegionStore implementation
//
// 18.09.02, G.Cosmo - Initial version
// --------------------------------------------------------------------

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4GeometryManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4ios.hh"
#include "G4AutoLock.hh"

namespace
{
  G4Mutex mapMutex = G4MUTEX_INITIALIZER;
}

// ***************************************************************************
// Static class variables
// ***************************************************************************
//
G4RegionStore* G4RegionStore::fgInstance = nullptr;
G4ThreadLocal G4VStoreNotifier* G4RegionStore::fgNotifier = nullptr;
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

  G4RegionStore* store = GetInstance();

  for(auto pos=store->cbegin(); pos!=store->cend(); ++pos)
  {
    if (fgNotifier != nullptr) { fgNotifier->NotifyDeRegistration(); }
    delete *pos;
  }

  store->bmap.clear(); store->mvalid = false;
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
// Bring contents of internal map up to date and reset validity flag
// ***************************************************************************
//
void G4RegionStore::UpdateMap()
{
  G4AutoLock l(&mapMutex);  // to avoid thread contention at initialisation
  if (mvalid) return;
  bmap.clear();
  for(auto pos=GetInstance()->cbegin(); pos!=GetInstance()->cend(); ++pos)
  {
    const G4String& reg_name = (*pos)->GetName();
    auto it = bmap.find(reg_name);
    if (it != bmap.cend())
    {
      it->second.push_back(*pos);
    }
    else
    {
      std::vector<G4Region*> reg_vec { *pos };
      bmap.insert(std::make_pair(reg_name, reg_vec));
    }
  }
  mvalid = true;
  l.unlock();
}

// ***************************************************************************
// Add Region to container
// ***************************************************************************
//
void G4RegionStore::Register(G4Region* pRegion)
{
  G4RegionStore* store = GetInstance();
  store->push_back(pRegion);
  const G4String& reg_name = pRegion->GetName();
  auto it = store->bmap.find(reg_name);
  if (it != store->bmap.cend())
  {
    it->second.push_back(pRegion);
  }
  else
  {
    std::vector<G4Region*> reg_vec { pRegion };
    store->bmap.insert(std::make_pair(reg_name, reg_vec));
  }
  if (fgNotifier) { fgNotifier->NotifyRegistration(); }
  store->mvalid = true;
}

// ***************************************************************************
// Remove Region from container
// ***************************************************************************
//
void G4RegionStore::DeRegister(G4Region* pRegion)
{
  G4RegionStore* store = GetInstance();
  if (!locked)    // Do not de-register if locked !
  {
    if (fgNotifier != nullptr)  { fgNotifier->NotifyDeRegistration(); }
    for (auto i=store->cbegin(); i!=store->cend(); ++i)
    {
      if (**i==*pRegion)
      {
        store->erase(i);
        break;
      }
    }
    const G4String& reg_name = pRegion->GetName();
    auto it = store->bmap.find(reg_name);
    if (it != store->bmap.cend())
    {
      if (it->second.size() > 1)
      {
        for (auto i=it->second.cbegin(); i!=it->second.cend(); ++i)
        {
          if (**i==*pRegion)
          {
            it->second.erase(i);
            break;
          }
        }
      }
      else
      {
        store->bmap.erase(it);
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
  if (fgInstance == nullptr)
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
  for (auto i=GetInstance()->cbegin(); i!=GetInstance()->cend(); ++i)
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
  for (auto i=GetInstance()->cbegin(); i!=GetInstance()->cend(); ++i)
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
  for (auto i=GetInstance()->cbegin(); i!=GetInstance()->cend(); ++i)
  {
    if((*i)->IsInMassGeometry() || (*i)->IsInParallelGeometry()
                                || (currentWorld != nullptr))
    { (*i)->UpdateMaterialList(); }
  }
}

// ***************************************************************************
// Returns a region through its name specification.
// ***************************************************************************
//
G4Region* G4RegionStore::GetRegion(const G4String& name, G4bool verbose) const
{
  G4RegionStore* store = GetInstance();
  if (!store->mvalid)  { store->UpdateMap(); }
  auto pos = store->bmap.find(name);
  if(pos != store->bmap.cend())
  {
    if ((verbose) && (pos->second.size()>1))
    {
      std::ostringstream message;
      message << "There exists more than ONE region in store named: "
              << name << "!" << G4endl
              << "Returning the first found.";
      G4Exception("G4RegionStore::GetSolid()",
                  "GeomMgt1001", JustWarning, message);
    }
    return pos->second[0];
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
  return nullptr;
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
  if (target == nullptr)
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
  for (auto i=GetInstance()->cbegin(); i!=GetInstance()->cend(); ++i)
  { (*i)->SetWorld(nullptr); }

  // Find world volumes
  //
  G4PhysicalVolumeStore* fPhysicalVolumeStore
    = G4PhysicalVolumeStore::GetInstance();
  std::size_t nPhys = fPhysicalVolumeStore->size();
  for(std::size_t iPhys=0; iPhys<nPhys; ++iPhys)
  {
    G4VPhysicalVolume* fPhys = (*fPhysicalVolumeStore)[iPhys];
    if(fPhys->GetMotherLogical() != nullptr) { continue; } // not a world volume

    // Now 'fPhys' is a world volume, set it to regions that belong to it.
    //
    for (auto i=GetInstance()->cbegin(); i!=GetInstance()->cend(); ++i)
    { (*i)->SetWorld(fPhys); }
  }
}

