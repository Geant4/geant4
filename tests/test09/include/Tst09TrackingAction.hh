
#ifndef Tst09TrackingAction_h
#define Tst09TrackingAction_h 1

#include "G4UserTrackingAction.hh"


class Tst09TrackingAction : public G4UserTrackingAction {

  public:
    Tst09TrackingAction(){};
    ~Tst09TrackingAction(){};
   
    void PreUserTrackingAction();

};

#endif
