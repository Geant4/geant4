//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4SolidStore.cc,v 1.9 2002-04-26 16:24:36 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4SolidStore
//
// Implementation for singleton container
//
// History:
// 10.07.95 P.Kent Initial version
// ********************************************************************

#include "globals.hh"
#include "G4SolidStore.hh"
#include "G4GeometryManager.hh"

// ***************************************************************************
// Static class variables
// ***************************************************************************
//
G4SolidStore* G4SolidStore::fgInstance = 0;
G4bool G4SolidStore::locked = false;

// ***************************************************************************
// Protected constructor: Construct underlying container with
// initial size of 100 entries
// ***************************************************************************
//
G4SolidStore::G4SolidStore()
  : G4std::vector<G4VSolid*>()
{
  reserve(100);
}

// ***************************************************************************
// Destructor
// ***************************************************************************
//
G4SolidStore::~G4SolidStore() 
{
  // NOTE: destruction of solids is client responsibility !
  // clear();
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
  if (G4GeometryManager::GetInstance()->IsGeometryClosed())
  {
    G4cout << "WARNING - Attempt to delete the solid store"
           << " while geometry closed !" << G4endl;
    return;
  }

  // Locks store for deletion of solids. De-registration will be
  // performed at this stage. G4VSolids will not de-register themselves.
  //
  locked = true;  

  size_t i=0;
  G4SolidStore* store = GetInstance();
  G4std::vector<G4VSolid*>::iterator pos;

#ifdef G4GEOMETRY_VOXELDEBUG
  G4cout << "Deleting Solids ... ";
#endif

  for(pos=store->begin(); pos!=store->end(); pos++)
  {
    if (*pos) delete *pos; i++;
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
// Add Solid to container
// ***************************************************************************
//
void G4SolidStore::Register(G4VSolid* pSolid)
{
  GetInstance()->push_back(pSolid);
}

// ***************************************************************************
// Remove Solid from container
// ***************************************************************************
//
void G4SolidStore::DeRegister(G4VSolid* pSolid)
{
  if (!locked)    // Do not de-register if locked !
  {
    for (iterator i=GetInstance()->begin(); i!=GetInstance()->end(); i++)
    {
      if (**i==*pSolid)
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
G4SolidStore* G4SolidStore::GetInstance()
{
  static G4SolidStore worldStore;
  if (!fgInstance)
  {
    fgInstance = &worldStore;
  }
  return fgInstance;
}
