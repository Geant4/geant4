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
// G4LogicalVolumeStore implementation
//
// 10.07.95, P.Kent, G.Cosmo
// --------------------------------------------------------------------

#include "G4Types.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4GeometryManager.hh"

#include "G4AutoLock.hh"

namespace
{
  G4Mutex mapMutex = G4MUTEX_INITIALIZER;
}

// ***************************************************************************
// Static class variables
// ***************************************************************************
//
G4LogicalVolumeStore* G4LogicalVolumeStore::fgInstance = nullptr;
G4ThreadLocal G4VStoreNotifier* G4LogicalVolumeStore::fgNotifier = nullptr;
G4ThreadLocal G4bool G4LogicalVolumeStore::locked = false;

// ***************************************************************************
// Protected constructor: Construct underlying container with
// initial size of 100 entries
// ***************************************************************************
//
G4LogicalVolumeStore::G4LogicalVolumeStore()
 : std::vector<G4LogicalVolume*>()
{
  reserve(100);
}

// ***************************************************************************
// Destructor
// ***************************************************************************
//
G4LogicalVolumeStore::~G4LogicalVolumeStore()
{
  Clean();  // Delete all volumes in the store
  G4LogicalVolume::Clean();  // Delete allocated sub-instance data
}

// ***************************************************************************
// Delete all elements from the store
// ***************************************************************************
//
void G4LogicalVolumeStore::Clean()
{
  // Do nothing if geometry is closed
  //
  if (G4GeometryManager::IsGeometryClosed())
  {
    G4cout << "WARNING - Attempt to delete the logical volume store"
           << " while geometry closed !" << G4endl;
    return;
  }

  // Locks store for deletion of volumes. De-registration will be
  // performed at this stage. G4LogicalVolumes will not de-register themselves.
  //
  locked = true;  

  G4LogicalVolumeStore* store = GetInstance();

  for(auto pos=store->cbegin(); pos!=store->cend(); ++pos)
  {
    if (fgNotifier != nullptr) { fgNotifier->NotifyDeRegistration(); }
    if (*pos != nullptr) { (*pos)->Lock(); delete *pos; }
  }

  store->bmap.clear(); store->mvalid = false;
  locked = false;
  store->clear();
}

// ***************************************************************************
// Associate user notifier to the store
// ***************************************************************************
//
void G4LogicalVolumeStore::SetNotifier(G4VStoreNotifier* pNotifier)
{
  GetInstance();
  fgNotifier = pNotifier;
}

// ***************************************************************************
// Bring contents of internal map up to date and reset validity flag
// ***************************************************************************
//
void G4LogicalVolumeStore::UpdateMap()
{
  G4AutoLock l(&mapMutex);  // to avoid thread contention at initialisation
  if (mvalid) return;
  bmap.clear();
  for(auto pos=GetInstance()->cbegin(); pos!=GetInstance()->cend(); ++pos)
  {
    const G4String& vol_name = (*pos)->GetName();
    auto it = bmap.find(vol_name);
    if (it != bmap.cend())
    {
      it->second.push_back(*pos);
    }
    else
    {
      std::vector<G4LogicalVolume*> vol_vec { *pos };
      bmap.insert(std::make_pair(vol_name, vol_vec));
    }
  }
  mvalid = true;
  l.unlock();
}

// ***************************************************************************
// Add volume to container
// ***************************************************************************
//
void G4LogicalVolumeStore::Register(G4LogicalVolume* pVolume)
{
  G4LogicalVolumeStore* store = GetInstance();
  store->push_back(pVolume);
  const G4String& vol_name = pVolume->GetName();
  auto it = store->bmap.find(vol_name);
  if (it != store->bmap.cend())
  {
    it->second.push_back(pVolume);
  }
  else
  {
    std::vector<G4LogicalVolume*> vol_vec { pVolume };
    store->bmap.insert(std::make_pair(vol_name, vol_vec));
  }
  if (fgNotifier) { fgNotifier->NotifyRegistration(); }
  store->mvalid = true;
}

// ***************************************************************************
// Remove volume from container
// ***************************************************************************
//
void G4LogicalVolumeStore::DeRegister(G4LogicalVolume* pVolume)
{
  G4LogicalVolumeStore* store = GetInstance();
  if (!locked)    // Do not de-register if locked !
  {
    if (fgNotifier != nullptr) { fgNotifier->NotifyDeRegistration(); }
    for (auto i=store->cbegin(); i!=store->cend(); ++i)
    {
      if (**i==*pVolume)
      {
        store->erase(i);
        break;
      }
    }
    const G4String& vol_name = pVolume->GetName();
    auto it = store->bmap.find(vol_name);
    if (it != store->bmap.cend())
    {
      if (it->second.size() > 1)
      {
        for (auto i=it->second.cbegin(); i!=it->second.cend(); ++i)
        {
          if (**i==*pVolume)
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
// Retrieve the first or last volume pointer in the container having that name
// ***************************************************************************
//
G4LogicalVolume*
G4LogicalVolumeStore::GetVolume(const G4String& name, G4bool verbose,
                                G4bool reverseSearch) const
{
  G4LogicalVolumeStore* store = GetInstance();
  if (!store->mvalid)  { store->UpdateMap(); }
  auto pos = store->bmap.find(name);
  if(pos != store->bmap.cend())
  {
    if ((verbose) && (pos->second.size()>1))
    {
      std::ostringstream message;
      message << "There exists more than ONE logical volume in store named: "
              << name << "!" << G4endl
              << "Returning the first found.";
      G4Exception("G4LogicalVolumeStore::GetVolume()",
                  "GeomMgt1001", JustWarning, message);
    }
    if(reverseSearch)
    {
      return pos->second[pos->second.size()-1];
    }
    else
    {
      return pos->second[0];
    }
  }
  if (verbose)
  {
     std::ostringstream message;
     message << "Volume NOT found in store !" << G4endl
             << "        Volume " << name << " NOT found in store !" << G4endl
             << "        Returning NULL pointer.";
     G4Exception("G4LogicalVolumeStore::GetVolume()",
                 "GeomMgt1001", JustWarning, message);
  }
  return nullptr;
}

// ***************************************************************************
// Return ptr to Store, setting if necessary
// ***************************************************************************
//
G4LogicalVolumeStore* G4LogicalVolumeStore::GetInstance()
{
  static G4LogicalVolumeStore worldStore;
  if (fgInstance == nullptr)
  {
    fgInstance = &worldStore;
  }
  return fgInstance;
}
