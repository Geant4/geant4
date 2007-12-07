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
// $Id: G4FieldManagerStore.cc,v 1.1 2007-12-07 15:34:10 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4FieldManagerStore
//
// Implementation for singleton container - copied from LogicalVolumeStore
//
// History:
// 07.12.07 J.Apostolakis Adapted for FieldManagerStore
// 10.07.95 P.Kent        Initial version (of LogicalVolumeStore)
// --------------------------------------------------------------------

#include "G4Types.hh"
#include "G4FieldManagerStore.hh"
// #include "G4GeometryManager.hh"

// ***************************************************************************
// Static class variables
// ***************************************************************************
//
G4FieldManagerStore* G4FieldManagerStore::fgInstance = 0;
// G4VStoreNotifier* G4FieldManagerStore::fgNotifier = 0;
G4bool G4FieldManagerStore::locked = false;

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
}

// ***************************************************************************
// Delete all elements from the store
// ***************************************************************************
//
void G4FieldManagerStore::Clean()
{
#if 0
  // Do nothing if geometry is closed
  //
  if (G4GeometryManager::GetInstance()->IsGeometryClosed())
  {
    G4cout << "WARNING - Attempt to delete the field manager store"
           << " while geometry closed !" << G4endl;
    return;
  }
#endif

  // Locks store for deletion of field managers. De-registration will be
  // performed at this stage. G4FieldManagers will not de-register themselves.
  //
  locked = true;  

  size_t i=0;
  G4FieldManagerStore* store = GetInstance();

#ifdef G4DEBUG
  G4cout << "Deleting Field Managers ... ";
#endif

  for(iterator pos=store->begin(); pos!=store->end(); pos++)
  {
    // if (fgNotifier) { fgNotifier->NotifyDeRegistration(); }
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

#if 0 
// ***************************************************************************
// Associate user notifier to the store
// ***************************************************************************
//
void G4FieldManagerStore::SetNotifier(G4VStoreNotifier* pNotifier)
{
  GetInstance();
  fgNotifier = pNotifier;
}
#endif

// ***************************************************************************
// Add field manager to container
// ***************************************************************************
//
void G4FieldManagerStore::Register(G4FieldManager* pFieldManager)
{
  GetInstance()->push_back(pFieldManager);
  // if (fgNotifier) { fgNotifier->NotifyRegistration(); }
}

// ***************************************************************************
// Remove volume from container
// ***************************************************************************
//
void G4FieldManagerStore::DeRegister(G4FieldManager* pFieldMgr)
{
  if (!locked)    // Do not de-register if locked !
  {
    // if (fgNotifier) { fgNotifier->NotifyDeRegistration(); }
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
  static G4FieldManagerStore worldStore;
  if (!fgInstance)
  {
    fgInstance = &worldStore;
  }
  return fgInstance;
}

// iterator<G4FieldManager*>  GetIterator()
// {
//  iterator *itFM; 
