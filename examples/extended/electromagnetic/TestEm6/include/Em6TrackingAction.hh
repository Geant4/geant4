// Em6TrackingAction.hh

#ifndef Em6TrackingAction_h
#define Em6TrackingAction_h 1

#include "G4UserTrackingAction.hh"

class Em6RunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Em6TrackingAction : public G4UserTrackingAction {

  public:
    Em6TrackingAction(Em6RunAction*);
   ~Em6TrackingAction() {};

    void PreUserTrackingAction(const G4Track* aTrack);
    void PostUserTrackingAction(const G4Track*);

  private:
    Em6RunAction* Em6Run;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
