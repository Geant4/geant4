
#ifndef LXeEventAction_h
#define LXeEventAction_h 1

#include "LXeEventMessenger.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4Event;
class RecorderBase;

class LXeEventAction : public G4UserEventAction
{
public:
  LXeEventAction(RecorderBase*);
  ~LXeEventAction();
  
public:
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);
  
  void SetSaveThreshold(G4int save);

  void SetEventVerbose(G4int v){verbose=v;}

  void SetPMTThreshold(G4int t){pmtThreshold=t;}

  void SetForceDrawPhotons(G4bool b){forcedrawphotons=b;}
  void SetForceDrawNoPhotons(G4bool b){forcenophotons=b;}

private:
  RecorderBase* recorder;
  LXeEventMessenger* eventMessenger;

  G4int              saveThreshold;

  G4int              scintCollID;
  G4int              pmtCollID;

  G4int              verbose;
  
  G4int              pmtThreshold;
  
  G4bool forcedrawphotons;
  G4bool forcenophotons;

};

#endif

    
