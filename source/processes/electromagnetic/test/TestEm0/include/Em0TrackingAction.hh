// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em0TrackingAction.hh,v 1.4 1999-12-15 14:51:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#ifndef Em0TrackingAction_h
#define Em0TrackingAction_h 1

#include "G4UserTrackingAction.hh"

class Em0RunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em0TrackingAction : public G4UserTrackingAction {

  public:  
    Em0TrackingAction(Em0RunAction*);
   ~Em0TrackingAction() {};
   
    void PostUserTrackingAction(const G4Track*);
    
  private:
    Em0RunAction* runAction;  
};

#endif
