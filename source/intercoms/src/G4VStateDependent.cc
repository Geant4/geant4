// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VStateDependent.cc,v 1.2 1999-12-15 14:50:42 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4VStateDependent ----------------
//             by Gabriele Cosmo, November 1996
// ------------------------------------------------------------

#include "G4VStateDependent.hh"
#include "G4StateManager.hh"

G4VStateDependent::G4VStateDependent() 
{
  G4StateManager * stateManager = G4StateManager::GetStateManager();
  stateManager->RegisterDependent(this);
}

G4VStateDependent::~G4VStateDependent()
{
  G4StateManager * stateManager = G4StateManager::GetStateManager();
  stateManager->DeregisterDependent(this);
}

G4VStateDependent::G4VStateDependent(const G4VStateDependent &right)
{
   *this = right;
}

G4VStateDependent& G4VStateDependent::operator=(const G4VStateDependent &right)
{
   *this = right;
   return *this;
}

G4int G4VStateDependent::operator==(const G4VStateDependent &right) const
{
   return (this == (G4VStateDependent *) &right);
}

G4int G4VStateDependent::operator!=(const G4VStateDependent &right) const
{
   return (this != (G4VStateDependent *) &right);
}
