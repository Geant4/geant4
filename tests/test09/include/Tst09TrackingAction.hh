
#ifndef Tst09TrackingAction_h
#define Tst09TrackingAction_h 1

#include "G4UserTrackingAction.hh"


class Tst09TrackingAction : public G4UserTrackingAction {

  public:
    Tst09TrackingAction(){};
    virtual ~Tst09TrackingAction(){};
   
    virtual void PreUserTrackingAction(const G4Track* aTrack);

};

#endif
