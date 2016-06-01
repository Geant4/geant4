// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HCtable.hh,v 2.1 1998/07/12 02:53:17 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//

#ifndef G4HCtable_H
#define G4HCtable_H 1

#include "globals.hh"
#include <rw/tvordvec.h>

class G4HCtable
{
  public:
    G4HCtable();
    ~G4HCtable();

  public:
    G4int Registor(G4String SDname,G4String HCname);
    G4int GetCollectionID(G4String HCname);

  private:
    RWTValOrderedVector<G4String> SDlist;
    RWTValOrderedVector<G4String> HClist;

  public:
    inline G4int entries() const
    { return HClist.entries(); }

};

#endif

