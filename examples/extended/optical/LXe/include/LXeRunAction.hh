#include "G4UserRunAction.hh"

#ifndef LXeRunAction_h
#define LXeRunAction_h 1

class RecorderBase;

class LXeRunAction : public G4UserRunAction
{
public:
  LXeRunAction(RecorderBase*);
  ~LXeRunAction();
  
  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);

private:
  RecorderBase* recorder;
};

#endif
