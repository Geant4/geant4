#ifndef B01MassImportanceScoring_hh
#define B01MassImportanceScoring_hh B01MassImportanceScoring_hh

#include "B01VSimulation.hh"

class G4VSampler;
class G4CellScorerStore;
class G4CellStoreScorer;
class G4IStore;
class B01SlobedConcreteShield;
class G4CellScorer;

class B01MassImportanceScoring : public B01VSimulation {
public:
  B01MassImportanceScoring();
  virtual ~B01MassImportanceScoring();
  virtual const G4String &GetName() const;

  virtual G4VPhysicalVolume &GetMassGeometry();
  virtual const G4CellScorer *GetG4CellScorer();
  virtual void PrepareSampling();
  virtual void ConfigureSampling();
  virtual void SetWeightRoulette(G4bool wroulette);

  virtual void PostRun(G4std::ostream *);
private:
  B01MassImportanceScoring(const B01MassImportanceScoring &);
  B01MassImportanceScoring &operator=(const B01MassImportanceScoring &);

  G4String fName;
  G4bool fWeightRoulette;
  B01SlobedConcreteShield *fGeometry;
  const G4CellScorer *fSpecialCellScorer;
  G4IStore *fIStore;
  G4CellScorerStore *fCS_store;
  G4CellStoreScorer *fScorer;
  G4VSampler *fSampler;
};

#endif
