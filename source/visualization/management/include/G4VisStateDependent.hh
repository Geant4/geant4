// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisStateDependent.hh,v 1.2 1999-12-07 11:06:59 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// A "state dependent" service class for G4VisManager.
// John Allison  29th November 1999

// G4VisStateDependent is "state dependent", i.e., it is notified on
// change of state (G4ApplicationState).  This is used to message the
// G4VisManager to draw hits and trajectories in the current scene at the
// end of event, as required.

#ifndef G4VISSTATEDEPENDENT_HH
#define G4VISSTATEDEPENDENT_HH

#include "G4VStateDependent.hh"

class G4VisManager;

class G4VisStateDependent: public G4VStateDependent {
  friend class G4VisManager;
private:
  G4VisStateDependent (G4VisManager *);
  G4bool Notify (G4ApplicationState requestedState);
  G4VisManager* fpVisManager;
};

#endif
