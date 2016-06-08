// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PersistentHitMan.cc,v 1.4 1999/11/28 21:54:16 morita Exp $
// GEANT4 tag $Name: geant4-01-01 $
//
// class G4PersistentHitMan 
//
// History:
// 99.11.16 Y.Morita  Initial version

#include "G4PersistentHitMan.hh"


G4PersistentHitMan* G4PersistentHitMan::f_PersistentHitMan = NULL;

G4PersistentHitMan* G4PersistentHitMan::GetPersistentHitMan()
{
  if(!f_PersistentHitMan)
  {
    f_PersistentHitMan = new G4PersistentHitMan;
  }
  return f_PersistentHitMan;
}

G4PersistentHitMan* G4PersistentHitMan::get_PersistentHitManIfExist()
{ return f_PersistentHitMan; }

G4PersistentHitMan::G4PersistentHitMan()
 : f_CurrentPHCofThisEvent(NULL)
{;}

G4PersistentHitMan::~G4PersistentHitMan()
{;}

