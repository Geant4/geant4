#ifndef B01ParallelScoring_hh
#define B01ParallelScoring_hh B01ParallelScoring_hh

#include "B01VSimulation.hh"

class B01VGeometry;
class B01Run;
class G4ParallelScoreSampler;
class G4CellScorerStore;
class G4CellStoreScorer;
class B01ParallelGeometry;


class B01ParallelScoring : public B01VSimulation {
public:
  B01ParallelScoring();
  ~B01ParallelScoring();
  G4String GetName() const;
  void Construct();
  void Run(G4int nevents);
  void PostRun(G4std::ostream *);
private:
  B01ParallelScoring(const B01ParallelScoring &);
  B01ParallelScoring &operator=
  (const B01ParallelScoring &);

  G4String fName;
  B01VGeometry *fMassGeometry;
  B01Run *fRun;
  B01ParallelGeometry *fParallelGeometry;
  G4CellScorerStore *fCS_store;
  G4CellStoreScorer *fScorer;
  G4ParallelScoreSampler *fSampler;
};

#endif
