// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VPersistencyManager.cc,v 1.1.10.1 1999/12/07 20:53:00 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//

#include "G4VPersistencyManager.hh"

G4VPersistencyManager* G4VPersistencyManager::fPersistencyManager = 0;

G4VPersistencyManager* G4VPersistencyManager::GetPersistencyManager()
{
  return fPersistencyManager;
}

G4VPersistencyManager::G4VPersistencyManager()
{
  fPersistencyManager = this;
}

G4VPersistencyManager::~G4VPersistencyManager()
{
  fPersistencyManager = NULL;
}


