#ifndef DetectorSliceTrackingAction_h 
#define DetectorSliceTrackingAction_h 1

#include "G4UserTrackingAction.hh"

class DetectorSliceTrackingAction : public G4UserTrackingAction {

public:

  DetectorSliceTrackingAction();

  ~DetectorSliceTrackingAction();

  void PreUserTrackingAction( const G4Track* );

  void PostUserTrackingAction( const G4Track* );

};

#endif

