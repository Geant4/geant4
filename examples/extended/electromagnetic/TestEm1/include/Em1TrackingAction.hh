// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em1TrackingAction.hh,v 1.1 1999-10-11 13:07:42 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#ifndef Em1TrackingAction_h
#define Em1TrackingAction_h 1

#include "G4UserTrackingAction.hh"

class Em1RunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em1TrackingAction : public G4UserTrackingAction {

  public:  
    Em1TrackingAction(Em1RunAction*);
   ~Em1TrackingAction() {};
   
    void PostUserTrackingAction(const G4Track*);
    
  private:
    Em1RunAction* runAction;  
};

#endif
