#ifndef LXeSteppingAction_H
#define LXeSteppingACtion_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

class RecorderBase;
class LXeEventAction;
class LXeTrackingAction;
class LXeSteppingMessenger;

class LXeSteppingAction : public G4UserSteppingAction
{
public:
  LXeSteppingAction(RecorderBase*);
  ~LXeSteppingAction();
  virtual void UserSteppingAction(const G4Step*);

  void SetOneStepPrimaries(G4bool b){oneStepPrimaries=b;}
  G4bool GetOneStepPrimaries(){return oneStepPrimaries;}
  
private:
  RecorderBase* recorder;
  G4bool oneStepPrimaries;
  LXeSteppingMessenger* steppingMessenger;
};

#endif
