#ifndef B01ParallelImportanceScoring_hh
#define B01ParallelImportanceScoring_hh B01ParallelImportanceScoring_hh

#include "B01VSimulation.hh"

class B01VGeometry;
class G4IStore;
class G4VSampler;
class G4CellScorerStore;
class G4CellStoreScorer;
class B01ParallelGeometry;
class G4CellScorer;

class B01ParallelImportanceScoring : public B01VSimulation {
public:
  B01ParallelImportanceScoring();
  virtual ~B01ParallelImportanceScoring();
  virtual const G4String &GetName() const;

  virtual G4VPhysicalVolume &GetMassGeometry();
  virtual const G4CellScorer *GetG4CellScorer();
  virtual void PrepareSampling();
  virtual void ConfigureSampling();
  virtual void SetWeightRoulette(G4bool wroulette);

  virtual void PostRun(G4std::ostream *);
private:
  B01ParallelImportanceScoring(const B01ParallelImportanceScoring &);
  B01ParallelImportanceScoring &operator=
  (const B01ParallelImportanceScoring &);

  G4String fName;
  G4bool fWeightRoulette;
  B01VGeometry *fMassGeometry;
  B01ParallelGeometry *fParallelGeometry;
  const G4CellScorer *fSpecialCellScorer;
  G4IStore *fIStore;
  G4CellScorerStore *fCS_store;
  G4CellStoreScorer *fScorer;
  G4VSampler *fSampler;
};

#endif
