// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VStateDependent.cc,v 1.2 2000-11-20 17:26:49 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      Geant4 Collaboration
//      ---------------- G4VStateDependent ----------------
//             by Gabriele Cosmo, November 1996
// ------------------------------------------------------------

#include "G4VStateDependent.hh"
#include "G4StateManager.hh"

G4VStateDependent::G4VStateDependent(G4bool bottom) 
{
  G4StateManager * stateManager = G4StateManager::GetStateManager();
  stateManager->RegisterDependent(this,bottom);
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
   if (&right == this) return *this;
   *this = right;
   return *this;
}

G4int G4VStateDependent::operator==(const G4VStateDependent &right) const
{
   return (this == &right);
}

G4int G4VStateDependent::operator!=(const G4VStateDependent &right) const
{
   return (this != &right);
}
