// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VPersistencyManager.cc,v 1.2 1999-12-15 14:53:54 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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


