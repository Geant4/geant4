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
// G4AssemblyStore
//
// Implementation for singleton container
//
// History:
// 9.10.18 G.Cosmo Initial version
// --------------------------------------------------------------------

#include "G4AssemblyVolume.hh"
#include "G4AssemblyStore.hh"
#include "G4GeometryManager.hh"
#include "G4ios.hh"

// ***************************************************************************
// Static class variables
// ***************************************************************************
//
G4AssemblyStore* G4AssemblyStore::fgInstance = nullptr;
G4ThreadLocal G4VStoreNotifier* G4AssemblyStore::fgNotifier = nullptr;
G4ThreadLocal G4bool G4AssemblyStore::locked = false;

// ***************************************************************************
// Protected constructor: Construct underlying container with
// initial size of 20 entries
// ***************************************************************************
//
G4AssemblyStore::G4AssemblyStore()
  : std::vector<G4AssemblyVolume*>()
{
  reserve(20);
}

// ***************************************************************************
// Destructor
// ***************************************************************************
//
G4AssemblyStore::~G4AssemblyStore() 
{
  Clean();  // Delete all assemblies in the store
}

// ***************************************************************************
// Delete all assemblies from the store
// ***************************************************************************
//
void G4AssemblyStore::Clean()
{
  // Do nothing if geometry is closed
  //
  if (G4GeometryManager::IsGeometryClosed())
  {
    G4cout << "WARNING - Attempt to delete the assembly store"
           << " while geometry closed !" << G4endl;
    return;
  }

  // Locks store for deletion of assemblies. De-registration will be
  // performed at this stage. Assemblies will not de-register themselves.
  //
  locked = true;  

  G4AssemblyStore* store = GetInstance();

  for(auto pos=store->cbegin(); pos!=store->cend(); ++pos)
  {
    if (fgNotifier != nullptr) { fgNotifier->NotifyDeRegistration(); }
    if (*pos) { delete *pos; }
  }

  locked = false;
  store->clear();
}

// ***************************************************************************
// Associate user notifier to the store
// ***************************************************************************
//
void G4AssemblyStore::SetNotifier(G4VStoreNotifier* pNotifier)
{
  GetInstance();
  fgNotifier = pNotifier;
}

// ***************************************************************************
// Add Assembly to container
// ***************************************************************************
//
void G4AssemblyStore::Register(G4AssemblyVolume* pAssembly)
{
  GetInstance()->push_back(pAssembly);
  if (fgNotifier != nullptr)  { fgNotifier->NotifyRegistration(); }
}

// ***************************************************************************
// Remove Assembly from container
// ***************************************************************************
//
void G4AssemblyStore::DeRegister(G4AssemblyVolume* pAssembly)
{
  if (!locked)    // Do not de-register if locked !
  {
    if (fgNotifier != nullptr)  { fgNotifier->NotifyDeRegistration(); }
    for (auto i=GetInstance()->cbegin(); i!=GetInstance()->cend(); ++i)
    {
      if (*i==pAssembly)
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
G4AssemblyStore* G4AssemblyStore::GetInstance()
{
  static G4AssemblyStore assemblyStore;
  if (fgInstance == nullptr)
  {
    fgInstance = &assemblyStore;
  }
  return fgInstance;
}

// ***************************************************************************
// Returns an assembly through its name specification.
// ***************************************************************************
//
G4AssemblyVolume*
G4AssemblyStore::GetAssembly(unsigned int id, G4bool verbose) const
{
  for (auto i=GetInstance()->cbegin(); i!=GetInstance()->cend(); ++i)
  {
    if ((*i)->GetAssemblyID() == id) { return *i; }
  }
  if (verbose)
  {
    std::ostringstream message;
    message << "Assembly NOT found in store !" << G4endl
            << "        Assembly " << id << " NOT found in store !" << G4endl
            << "        Returning NULL pointer.";
    G4Exception("G4AssemblyStore::GetAssembly()",
                "GeomVol1001", JustWarning, message);
  }
  return 0;
}
