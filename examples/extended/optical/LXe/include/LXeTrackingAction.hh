
#ifndef LXeTrackingAction_h
#define LXeTrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "globals.hh"

class RecorderBase;

class LXeTrackingAction : public G4UserTrackingAction {

public:  
  LXeTrackingAction(RecorderBase*);
  ~LXeTrackingAction() {};
  
  void PreUserTrackingAction(const G4Track*);
  void PostUserTrackingAction(const G4Track*);
  
private:
  RecorderBase* recorder;
  
};

#endif
