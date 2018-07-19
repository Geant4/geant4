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
// $Id: G4FieldManagerStore.cc 103228 2017-03-22 14:52:32Z gcosmo $
//
// G4FieldManagerStore
//
// Implementation for singleton container
//
// History:
// 07.12.07 J.Apostolakis Adapted from G4LogicalVolumeStore
// --------------------------------------------------------------------

#include "G4Types.hh"
#include "G4FieldManagerStore.hh"
#include "G4ChordFinder.hh" 

// ***************************************************************************
// Static class variables
// ***************************************************************************
//
G4ThreadLocal G4FieldManagerStore* G4FieldManagerStore::fgInstance = 0;
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
  fgInstance = 0;
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

  size_t i=0;
  G4FieldManagerStore* store = GetInstance();

  for(iterator pos=store->begin(); pos!=store->end(); pos++)
  {
    if (*pos) { delete *pos; }
    i++;
  }

#ifdef G4GEOMETRY_DEBUG
  if (store->size() < i-1)
    { G4cout << "No field managers deleted. Already deleted by user ?" << G4endl; }
  else
    { G4cout << i-1 << " field managers deleted !" << G4endl; }
#endif

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
    for (iterator i=GetInstance()->begin(); i!=GetInstance()->end(); i++)
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
  if (!fgInstance)
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
  G4ChordFinder *pChordFnd;
   
  for (iterator i=GetInstance()->begin(); i!=GetInstance()->end(); i++)
  {
    pChordFnd = (*i)->GetChordFinder();
    if( pChordFnd )
    {
      pChordFnd->ResetStepEstimate();
    }
  }
}
