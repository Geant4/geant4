#ifndef B01ParallelImportanceScoring_hh
#define B01ParallelImportanceScoring_hh B01ParallelImportanceScoring_hh

#include "B01VSimulation.hh"

class B01VGeometry;
class B01Run;
class G4IStore;
class G4ParallelImportanceScoreSampler;
class G4CellScorerStore;
class G4CellStoreScorer;
class B01ParallelGeometry;


class B01ParallelImportanceScoring : public B01VSimulation {
public:
  B01ParallelImportanceScoring();
  ~B01ParallelImportanceScoring();
  G4String GetName() const;
  void Construct();
  void Run(G4int nevents);
  void PostRun(G4std::ostream *);
private:
  B01ParallelImportanceScoring(const B01ParallelImportanceScoring &);
  B01ParallelImportanceScoring &operator=
  (const B01ParallelImportanceScoring &);

  G4String fName;
  B01VGeometry *fMassGeometry;
  B01Run *fRun;
  B01ParallelGeometry *fParallelGeometry;
  G4IStore *fIStore;
  G4CellScorerStore *fCS_store;
  G4CellStoreScorer *fScorer;
  G4ParallelImportanceScoreSampler *fSampler;
};

#endif
