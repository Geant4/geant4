#ifndef G4ParallelSamplerMessenger_hh
#define G4ParallelSamplerMessenger_hh G4ParallelSamplerMessenger_hh

#include "G4UImessenger.hh"
#include "g4std/map"


class G4VSampler;
class G4ImportanceGeometryConstructor;
class G4UIcmdWithAString;
class G4ParallelImportanceScoreSampler;
class G4ParallelScoreSampler;
class G4StandardScorer;

typedef map<G4String, G4VSampler *> G4MapNameSampler;
typedef map<G4String, G4StandardScorer *> G4MapNameScorer;

class G4ParallelSamplerMessenger : public G4UImessenger {
public:
  G4ParallelSamplerMessenger(G4ImportanceGeometryConstructor *igc);
  ~G4ParallelSamplerMessenger(){}
  
  void SetNewValue(G4UIcommand * command, G4String newValue);


private:
  void CheckNewParticle(const G4String &particlename);
  void Error(const G4String &m){
    G4Exception("Error: G4ParallelSamplerMessenger:" + m);
  }
  void ImpCmdAction(const G4String &particlename);
  void ScoreCmdAction(const G4String &particlename);
  void PrintCmdAction(const G4String &particlename);
  

  G4ImportanceGeometryConstructor *fImpGeoConst;

  G4UIcmdWithAString *fImpCmd;
  G4UIcmdWithAString *fScoreCmd;
  G4UIcmdWithAString *fPrintCmd;

  G4ParallelImportanceScoreSampler *fImpSampler;
  G4StandardScorer *fImpScorer;
  G4String fImpPartilce;
  G4MapNameScorer fMapNameScorer;
  G4MapNameSampler fMapNameSampler;

};

#endif


