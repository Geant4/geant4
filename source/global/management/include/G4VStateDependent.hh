// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VStateDependent.hh,v 1.2 2000-11-20 17:26:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      For information related to this code contact:
//      Geant4 Collaboration
//
//      ---------------- G4VStateDependent ----------------
//
// Authors: G.Cosmo, M.Asai - November 1996
//
// ------------------------------------------------------------
//
// Class description:
//
// Abstract base class of all classes which need to be notified when
// the state of Geant4 changes. The concrete class object derived from
// this base class will be automatically registered to G4StateManager
// and the virtual method Notify() will be invoked when the state changes.

// ------------------------------------------------------------

#ifndef G4VStateDependent_h
#define G4VStateDependent_h 1

#include "globals.hh"
#include "G4ApplicationState.hh"

class G4VStateDependent
{

public:

  G4VStateDependent(G4bool bottom=false);
  virtual ~G4VStateDependent();
  G4int operator==(const G4VStateDependent &right) const;
  G4int operator!=(const G4VStateDependent &right) const;

public: // with description

  virtual G4bool Notify(G4ApplicationState requestedState) = 0;
    // Pure virtual method which will be invoked by G4StateManager.
    // In case state change must not be allowed by some reason of the
    // concrete class, false should be returned. But this scheme is
    // NOT recommended to use. All command which are state sensitive
    // MUST assign available state(s).

private:

  G4VStateDependent(const G4VStateDependent &right);
  G4VStateDependent& operator=(const G4VStateDependent &right);

};

#endif
