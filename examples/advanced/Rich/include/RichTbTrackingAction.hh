// Rich advanced example for Geant4
// RichTbTrackingAction.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbTrackingAction_h
#define RichTbTrackingAction_h 1

#include "G4UserTrackingAction.hh"
class RichTbTrackingAction : public G4UserTrackingAction {

 public:
    RichTbTrackingAction();
  virtual ~RichTbTrackingAction();
  void PreUserTrackingAction(const G4Track* aTrack){;}
  void PostUserTrackingAction(const G4Track* aTrack){;}


private:
};
#endif
