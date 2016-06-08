// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PersistentDigitMan.cc,v 1.4 1999/11/28 21:54:14 morita Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
// class G4PersistentDigitMan 
//
// History:
// 99.11.22 Y.Morita  Initial version

#include "G4PersistentDigitMan.hh"


G4PersistentDigitMan* G4PersistentDigitMan::f_PersistentDigitMan = NULL;

G4PersistentDigitMan* G4PersistentDigitMan::GetPersistentDigitMan()
{
  if(!f_PersistentDigitMan)
  {
    f_PersistentDigitMan = new G4PersistentDigitMan;
  }
  return f_PersistentDigitMan;
}

G4PersistentDigitMan* G4PersistentDigitMan::get_PersistentDigitManIfExist()
{ return f_PersistentDigitMan; }

G4PersistentDigitMan::G4PersistentDigitMan()
 : f_CurrentPDCofThisEvent(NULL)
{;}

G4PersistentDigitMan::~G4PersistentDigitMan()
{;}

