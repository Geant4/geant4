// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SolidStore.cc,v 1.6 2001-04-20 20:13:55 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4SolidStore
//
// Implementation for singleton container
//
// History:
// 10.07.95 P.Kent Initial version

#include "G4SolidStore.hh"
#include "globals.hh"

// Protected constructor: Construct underlying container with
// initial size of 100 entries
G4SolidStore::G4SolidStore()
  : G4std::vector<G4VSolid*>()
{
  reserve(100);
}

// Destructor
G4SolidStore::~G4SolidStore() 
{
  while (!empty())
  {
//    delete front();
    erase(begin());
  }
}

// Static class variable
G4SolidStore* G4SolidStore::fgInstance = 0;

// Add Solid to container
void G4SolidStore::Register(G4VSolid* pSolid)
{
    GetInstance()->push_back(pSolid);
}

// Remove Solid from container
void G4SolidStore::DeRegister(G4VSolid* pSolid)
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

// Return ptr to Store, setting if necessary
G4SolidStore* G4SolidStore::GetInstance()
{
    static G4SolidStore worldStore;
    if (!fgInstance)
	{
	    fgInstance = &worldStore;
	}
    return fgInstance;
}
