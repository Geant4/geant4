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
// $Id: G4LogicalVolumeStore.cc,v 1.16 2006/06/29 18:33:30 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// G4LogicalVolumeStore
//
// Implementation for singleton container
//
// History:
// 10.07.95 P.Kent Initial version
// --------------------------------------------------------------------

#include "G4Types.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4GeometryManager.hh"

// ***************************************************************************
// Static class variables
// ***************************************************************************
//
G4LogicalVolumeStore* G4LogicalVolumeStore::fgInstance = 0;
G4VStoreNotifier* G4LogicalVolumeStore::fgNotifier = 0;
G4bool G4LogicalVolumeStore::locked = false;

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
  Clean();
}

// ***************************************************************************
// Delete all elements from the store
// ***************************************************************************
//
void G4LogicalVolumeStore::Clean()
{
  // Do nothing if geometry is closed
  //
  if (G4GeometryManager::GetInstance()->IsGeometryClosed())
  {
    G4cout << "WARNING - Attempt to delete the logical volume store"
           << " while geometry closed !" << G4endl;
    return;
  }

  // Locks store for deletion of volumes. De-registration will be
  // performed at this stage. G4LogicalVolumes will not de-register themselves.
  //
  locked = true;  

  size_t i=0;
  G4LogicalVolumeStore* store = GetInstance();
  std::vector<G4LogicalVolume*>::iterator pos;

#ifdef G4GEOMETRY_VOXELDEBUG
  G4cout << "Deleting Logical Volumes ... ";
#endif

  for(pos=store->begin(); pos!=store->end(); pos++)
  {
    if (fgNotifier) fgNotifier->NotifyDeRegistration();
    if (*pos) delete *pos; i++;
  }

#ifdef G4GEOMETRY_VOXELDEBUG
  if (store->size() < i-1)
    { G4cout << "No volumes deleted. Already deleted by user ?" << G4endl; }
  else
    { G4cout << i-1 << " volumes deleted !" << G4endl; }
#endif

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
// Add volume to container
// ***************************************************************************
//
void G4LogicalVolumeStore::Register(G4LogicalVolume* pVolume)
{
  GetInstance()->push_back(pVolume);
  if (fgNotifier) fgNotifier->NotifyRegistration();
}

// ***************************************************************************
// Remove volume from container
// ***************************************************************************
//
void G4LogicalVolumeStore::DeRegister(G4LogicalVolume* pVolume)
{
  if (!locked)    // Do not de-register if locked !
  {
    if (fgNotifier) fgNotifier->NotifyDeRegistration();
    for (iterator i=GetInstance()->begin(); i!=GetInstance()->end(); i++)
    {
      if (**i==*pVolume)
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
G4LogicalVolumeStore* G4LogicalVolumeStore::GetInstance()
{
  static G4LogicalVolumeStore worldStore;
  if (!fgInstance)
  {
    fgInstance = &worldStore;
  }
  return fgInstance;
}
