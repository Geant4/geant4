#ifndef B01NGapplication_hh
#define B01NGapplication_hh B01NGapplication_hh

#include "B01VSimulation.hh"

class B01VGeometry;
class B01Run;
class G4IStore;
class G4ParallelScoreSampler;
class G4PImportanceWWindowScoreSampler;
class G4CellScorerStore;
class G4CellStoreScorer;
class B01ParallelGeometry;
class G4WeightWindowAlgorithm;

class B01NGapplication : public B01VSimulation {
public:
  B01NGapplication();
  ~B01NGapplication();
  G4String GetName() const;
  void Construct();
  void Run(G4int nevents);
  void PostRun(G4std::ostream *);
private:
  B01NGapplication(const B01NGapplication &);
  B01NGapplication &operator=
  (const B01NGapplication &);

  G4String fName;
  B01VGeometry *fMassGeometry;
  B01Run *fRun;
  B01ParallelGeometry *fParallelGeometry;
  G4IStore *fIStore;
  G4CellScorerStore *fCS_store_n;
  G4CellScorerStore *fCS_store_g; 
  G4CellStoreScorer *fScorer_n;
  G4CellStoreScorer *fScorer_g;  
  G4WeightWindowAlgorithm *fWWalg;
  G4PImportanceWWindowScoreSampler *fSampler_n;
  G4ParallelScoreSampler *fSampler_g;  
};

#endif
