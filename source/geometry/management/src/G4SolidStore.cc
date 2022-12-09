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
// G4SolidStore implementation
//
// 10.07.95, P.Kent, G.Cosmo
// --------------------------------------------------------------------

#include "globals.hh"
#include "G4SolidStore.hh"
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
G4SolidStore* G4SolidStore::fgInstance = nullptr;
G4ThreadLocal G4VStoreNotifier* G4SolidStore::fgNotifier = nullptr;
G4ThreadLocal G4bool G4SolidStore::locked = false;

// ***************************************************************************
// Protected constructor: Construct underlying container with
// initial size of 100 entries
// ***************************************************************************
//
G4SolidStore::G4SolidStore()
  : std::vector<G4VSolid*>()
{
  reserve(100);
}

// ***************************************************************************
// Destructor
// ***************************************************************************
//
G4SolidStore::~G4SolidStore() 
{
  Clean();
}

// ***************************************************************************
// Delete all elements from the store
// ***************************************************************************
//
void G4SolidStore::Clean()
{
  // Do nothing if geometry is closed
  //
  if (G4GeometryManager::IsGeometryClosed())
  {
    G4cout << "WARNING - Attempt to delete the solid store"
           << " while geometry closed !" << G4endl;
    return;
  }

  // Locks store for deletion of solids. De-registration will be
  // performed at this stage. G4VSolids will not de-register themselves.
  //
  locked = true;  

  G4SolidStore* store = GetInstance();

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
void G4SolidStore::SetNotifier(G4VStoreNotifier* pNotifier)
{
  GetInstance();
  fgNotifier = pNotifier;
}

// ***************************************************************************
// Bring contents of internal map up to date and reset validity flag
// ***************************************************************************
//
void G4SolidStore::UpdateMap()
{
  G4AutoLock l(&mapMutex);  // to avoid thread contention at initialisation
  if (mvalid) return;
  bmap.clear();
  for(auto pos=GetInstance()->cbegin(); pos!=GetInstance()->cend(); ++pos)
  {
    const G4String& sol_name = (*pos)->GetName();
    auto it = bmap.find(sol_name);
    if (it != bmap.cend())
    {
      it->second.push_back(*pos);
    }
    else
    {
      std::vector<G4VSolid*> sol_vec { *pos };
      bmap.insert(std::make_pair(sol_name, sol_vec));
    }
  }
  mvalid = true;
  l.unlock();
}

// ***************************************************************************
// Add Solid to container
// ***************************************************************************
//
void G4SolidStore::Register(G4VSolid* pSolid)
{
  G4SolidStore* store = GetInstance();
  store->push_back(pSolid);
  const G4String& sol_name = pSolid->GetName();
  auto it = store->bmap.find(sol_name);
  if (it != store->bmap.cend())
  {
    it->second.push_back(pSolid);
  }
  else
  {
    std::vector<G4VSolid*> sol_vec { pSolid };
    store->bmap.insert(std::make_pair(sol_name, sol_vec));
  }
  if (fgNotifier) { fgNotifier->NotifyRegistration(); }
  store->mvalid = true;
}

// ***************************************************************************
// Remove Solid from container
// ***************************************************************************
//
void G4SolidStore::DeRegister(G4VSolid* pSolid)
{
  G4SolidStore* store = GetInstance();
  if (!locked)    // Do not de-register if locked !
  {
    if (fgNotifier != nullptr) { fgNotifier->NotifyDeRegistration(); }
    for (auto i=store->crbegin(); i!=store->crend(); ++i)
    {
      if (**i==*pSolid)
      {
        store->erase(std::next(i).base());
        store->mvalid = false;
        break;
      }
    }
    const G4String& sol_name = pSolid->GetName();
    auto it = store->bmap.find(sol_name);
    if (it != store->bmap.cend())
    {
      if (it->second.size() > 1)
      {
        for (auto i=it->second.cbegin(); i!=it->second.cend(); ++i)
        {
          if (**i==*pSolid)
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
// Retrieve the first or last solid pointer in the container having that name
// ***************************************************************************
//
G4VSolid* G4SolidStore::GetSolid(const G4String& name, G4bool verbose,
                                 G4bool reverseSearch) const
{
  G4SolidStore* store = GetInstance();
  if (!store->mvalid)  { store->UpdateMap(); }
  auto pos = store->bmap.find(name);
  if(pos != store->bmap.cend())
  {
    if ((verbose) && (pos->second.size()>1))
    {
      std::ostringstream message;
      message << "There exists more than ONE solid in store named: "
              << name << "!" << G4endl
              << "Returning the first found.";
      G4Exception("G4SolidStore::GetSolid()",
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
    message << "Solid " << name << " not found in store !" << G4endl
            << "Returning NULL pointer.";
    G4Exception("G4SolidStore::GetSolid()",
                "GeomMgt1001", JustWarning, message);
  }
  return nullptr;
}

// ***************************************************************************
// Return ptr to Store, setting if necessary
// ***************************************************************************
//
G4SolidStore* G4SolidStore::GetInstance()
{
  static G4SolidStore worldStore;
  if (fgInstance == nullptr)
  {
    fgInstance = &worldStore;
  }
  return fgInstance;
}
