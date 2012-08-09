#ifndef Tst68TrackingAction_h 
#define Tst68TrackingAction_h 1

#include "G4UserTrackingAction.hh"

class Tst68TrackingAction : public G4UserTrackingAction {

public:

  Tst68TrackingAction();

  ~Tst68TrackingAction();

  void PreUserTrackingAction( const G4Track* );

  void PostUserTrackingAction( const G4Track* );

};

#endif
