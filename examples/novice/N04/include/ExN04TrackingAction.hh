
#ifndef ExN04TrackingAction_h
#define ExN04TrackingAction_h 1

#include "G4UserTrackingAction.hh"


class ExN04TrackingAction : public G4UserTrackingAction {

  public:
    ExN04TrackingAction(){};
    ~ExN04TrackingAction(){};
   
    void PreUserTrackingAction();

};

#endif
