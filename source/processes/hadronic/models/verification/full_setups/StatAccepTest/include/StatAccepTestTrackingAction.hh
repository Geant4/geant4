#ifndef StatAccepTestTrackingAction_h 
#define StatAccepTestTrackingAction_h 1

#include "G4UserTrackingAction.hh"

class StatAccepTestTrackingAction : public G4UserTrackingAction {

public:

  StatAccepTestTrackingAction();

  ~StatAccepTestTrackingAction();

  void PreUserTrackingAction( const G4Track* );

  void PostUserTrackingAction( const G4Track* );

};

#endif

