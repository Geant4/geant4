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
// G4FieldManagerStore implementation
//
// Author: J.Apostolakis, 07.12.2007 - Adapted from G4LogicalVolumeStore
// --------------------------------------------------------------------

#include "G4Types.hh"
#include "G4FieldManagerStore.hh"
#include "G4ChordFinder.hh" 

// ***************************************************************************
// Static class variables
// ***************************************************************************
//
G4ThreadLocal G4FieldManagerStore* G4FieldManagerStore::fgInstance = nullptr;
G4ThreadLocal G4bool G4FieldManagerStore::locked = false;

// ***************************************************************************
// Protected constructor: Construct underlying container with
// initial size of 100 entries
// ***************************************************************************
//
G4FieldManagerStore::G4FieldManagerStore()
 : std::vector<G4FieldManager*>()
{
  reserve(100);
}

// ***************************************************************************
// Destructor
// ***************************************************************************
//
G4FieldManagerStore::~G4FieldManagerStore()
{
  Clean();
  fgInstance = nullptr;
}

// ***************************************************************************
// Delete all elements from the store
// ***************************************************************************
//
void G4FieldManagerStore::Clean()
{
  // Locks store for deletion of field managers. De-registration will be
  // performed at this stage. G4FieldManagers will not de-register themselves.
  //
  locked = true;  

  G4FieldManagerStore* store = GetInstance();

  for(auto pos=store->cbegin(); pos!=store->cend(); ++pos)
  {
    if (*pos) { delete *pos; }
  }

  locked = false;
  store->clear();
}

// ***************************************************************************
// Add field manager to container
// ***************************************************************************
//
void G4FieldManagerStore::Register(G4FieldManager* pFieldManager)
{
  GetInstance()->push_back(pFieldManager);
}

// ***************************************************************************
// Remove volume from container
// ***************************************************************************
//
void G4FieldManagerStore::DeRegister(G4FieldManager* pFieldMgr)
{
  if (!locked)    // Do not de-register if locked !
  {
    for (auto i=GetInstance()->cbegin(); i!=GetInstance()->cend(); ++i)
    {
      if (*i==pFieldMgr)  //   For LogVol was **i == *pLogVolume ... Reason?
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
G4FieldManagerStore* G4FieldManagerStore::GetInstance()
{
  if (fgInstance == nullptr)
  {
    fgInstance = new G4FieldManagerStore;
  }
  return fgInstance;
}

// ***************************************************************************
// Return ptr to Store
// ***************************************************************************
//
G4FieldManagerStore* G4FieldManagerStore::GetInstanceIfExist()
{
  return fgInstance;
}

// ***************************************************************************
// Globally reset the state
// ***************************************************************************
//
void
G4FieldManagerStore::ClearAllChordFindersState()
{
  G4ChordFinder* pChordFnd;
   
  for (auto i=GetInstance()->cbegin(); i!=GetInstance()->cend(); ++i)
  {
    pChordFnd = (*i)->GetChordFinder();
    if( pChordFnd != nullptr )
    {
      pChordFnd->ResetStepEstimate();
    }
  }
}
