
#ifndef ExN04TrackingAction_h
#define ExN04TrackingAction_h 1

#include "G4UserTrackingAction.hh"


class ExN04TrackingAction : public G4UserTrackingAction {

  public:
    ExN04TrackingAction(){};
    virtual ~ExN04TrackingAction(){};
   
    virtual void PreUserTrackingAction(const G4Track*);

};

#endif
