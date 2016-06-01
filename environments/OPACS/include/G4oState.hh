// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4oState.hh,v 2.1 1998/07/12 02:37:02 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//
#ifndef G4oState_h
#define G4oState_h 

#include <G4VStateDependent.hh>
#include <G4UIsession.hh>

class G4oState : public G4VStateDependent, public G4UIsession {
public:
           G4oState          (G4String);
  G4bool   Notify            (G4ApplicationState);
  void     PauseSessionStart (G4String);
private:
  G4String name;
};

#endif

