// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HCtable.hh,v 1.2.2.1.2.1 1999/12/07 20:47:41 gunter Exp $
// GEANT4 tag $Name: geant4-03-00 $
//

#ifndef G4HCtable_H
#define G4HCtable_H 1

#include "globals.hh"
#include "g4rw/tvordvec.h"

// class description:
//
//  This class is used by G4SDManager for book keeping the
// sensitive detector modules and hits collections. The order of
// hits collections stored in G4HCofThisEvent is same as the
// order of HClist. 
//  The order may vary from run to run, if the user adds/changes
// some of his/her sensitive detector modules.
//  In case user wants to make G4Run object persistent, this
// G4HCtable class object should be copied and stored with
// G4Run object.

class G4HCtable
{
  public:
    G4HCtable();
    ~G4HCtable();

  public:
    G4int Registor(G4String SDname,G4String HCname);
    G4int GetCollectionID(G4String HCname);

  private:
    G4RWTValOrderedVector<G4String> SDlist;
    G4RWTValOrderedVector<G4String> HClist;

  public:
    inline G4int entries() const
    { return HClist.entries(); }

};

#endif

