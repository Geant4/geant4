
#ifndef Tst14TrackingAction_h
#define Tst14TrackingAction_h 1

#include "G4UserTrackingAction.hh"


class Tst14TrackingAction : public G4UserTrackingAction {

  public:
    Tst14TrackingAction(){};
    virtual ~Tst14TrackingAction(){};
   
    virtual void PreUserTrackingAction(const G4Track* aTrack);

};

#endif
