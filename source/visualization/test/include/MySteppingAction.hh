// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MySteppingAction.hh,v 1.3 1999-12-15 14:54:33 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// 16/Apr/1997  J. Allison:  For visualization/test/test19.

#ifndef MYSTEPPINGACTION_HH
#define MYSTEPPINGACTION_HH

#include "G4UserSteppingAction.hh"

#include "Randomize.hh"

class MySteppingAction: public G4UserSteppingAction {
  void UserSteppingAction(const G4Step*);
  private:
    HepJamesRandom theJamesEngine;
    DRand48Engine theDRand48Engine;
};

#endif
