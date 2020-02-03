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
// G4SolidStore implementation for singleton container
//
// 18.04.01, G.Cosmo - Migrated to STL vector
// 10.07.95, P.Kent - Initial version
// --------------------------------------------------------------------

#include "globals.hh"
#include "G4SolidStore.hh"
#include "G4GeometryManager.hh"

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

  size_t i = 0;
  G4SolidStore* store = GetInstance();

#ifdef G4GEOMETRY_VOXELDEBUG
  G4cout << "Deleting Solids ... ";
#endif

  for(auto pos=store->cbegin(); pos!=store->cend(); ++pos)
  {
    if (fgNotifier != nullptr) { fgNotifier->NotifyDeRegistration(); }
    delete *pos; ++i;
  }

#ifdef G4GEOMETRY_VOXELDEBUG
  if (store->size() < i-1)
    { G4cout << "No solids deleted. Already deleted by user ?" << G4endl; }
  else
    { G4cout << i-1 << " solids deleted !" << G4endl; }
#endif

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
// Add Solid to container
// ***************************************************************************
//
void G4SolidStore::Register(G4VSolid* pSolid)
{
  GetInstance()->push_back(pSolid);
  if (fgNotifier != nullptr) { fgNotifier->NotifyRegistration(); }
}

// ***************************************************************************
// Remove Solid from container
// ***************************************************************************
//
void G4SolidStore::DeRegister(G4VSolid* pSolid)
{
  if (!locked)    // Do not de-register if locked !
  {
    if (fgNotifier != nullptr) { fgNotifier->NotifyDeRegistration(); }
    for (auto i=GetInstance()->crbegin(); i!=GetInstance()->crend(); ++i)
    {
      if (**i==*pSolid)
      {
        GetInstance()->erase(std::next(i).base());
        break;
      }
    }
  }
}

// ***************************************************************************
// Retrieve the first solid pointer in the container having that name
// ***************************************************************************
//
G4VSolid* G4SolidStore::GetSolid(const G4String& name, G4bool verbose) const
{
  for (auto i=GetInstance()->cbegin(); i!=GetInstance()->cend(); ++i)
  {
    if ((*i)->GetName() == name) { return *i; }
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
