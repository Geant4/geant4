// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DCtable.hh,v 1.2 1999-11-01 03:20:35 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4DCtable_H
#define G4DCtable_H 1

#include "globals.hh"
#include <rw/tvordvec.h>

// class description:
//
//  This class is used by G4DigiManager for book keeping the
// digitizer modules and digits collections. The order of
// digi collections stored in G4DCofThisEvent is same as the
// order of DClist. 
//  The order may vary from run to run, if the user adds/changes
// some of his/her digitizer modules.
//  In case user wants to make G4Run object persistent, this
// G4DCtable class object should be copied and stored with
// G4Run object.

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

