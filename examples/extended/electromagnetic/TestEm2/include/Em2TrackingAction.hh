// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em2TrackingAction.hh,v 1.2 1999-12-15 14:49:00 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#ifndef Em2TrackingAction_h
#define Em2TrackingAction_h 1

#include "G4UserTrackingAction.hh"

class Em2RunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em2TrackingAction : public G4UserTrackingAction {

  public:  
    Em2TrackingAction(Em2RunAction*);
   ~Em2TrackingAction() {};
   
    void PostUserTrackingAction(const G4Track*);
    
  private:
    Em2RunAction* Em2Run;  
};

#endif
