// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DCtable.hh,v 1.1 1999-01-07 16:14:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4DCtable_H
#define G4DCtable_H 1

#include "globals.hh"
#include <rw/tvordvec.h>

class G4DCtable
{
  public:
    G4DCtable();
    ~G4DCtable();

  public:
    G4int Registor(G4String SDname,G4String DCname);
    G4int GetCollectionID(G4String DCname);

  private:
    RWTValOrderedVector<G4String> DMlist;
    RWTValOrderedVector<G4String> DClist;

  public:
    inline G4int entries() const
    { return DClist.entries(); }

};

#endif

