#ifndef B01MassImportanceScoring_hh
#define B01MassImportanceScoring_hh B01MassImportanceScoring_hh

#include "B01VSimulation.hh"

class B01Run;
class G4MassImportanceScoreSampler;
class G4CellScorerStore;
class G4CellStoreScorer;
class G4IStore;
class B01SlobedConcreteShield;

class B01MassImportanceScoring : public B01VSimulation {
public:
  B01MassImportanceScoring();
  ~B01MassImportanceScoring();
  G4String GetName() const;
  void Construct();
  void Run(G4int nevents);
  void PostRun(G4std::ostream *);
private:
  B01MassImportanceScoring(const B01MassImportanceScoring &);
  B01MassImportanceScoring &operator=(const B01MassImportanceScoring &);

  G4String fName;
  B01SlobedConcreteShield *fGeometry;
  B01Run *fRun;
  G4IStore *fIStore;
  G4CellScorerStore *fCS_store;
  G4CellStoreScorer *fScorer;
  G4MassImportanceScoreSampler *fSampler;
};

#endif
