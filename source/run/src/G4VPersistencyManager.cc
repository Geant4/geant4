// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VPersistencyManager.cc,v 1.1 1999-01-07 16:14:19 gunter Exp $
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


