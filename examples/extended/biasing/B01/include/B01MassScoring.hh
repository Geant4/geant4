#ifndef B01MassScoring_hh
#define B01MassScoring_hh B01MassScoring_hh

#include "B01VSimulation.hh"


class B01SlobedConcreteShield;
class G4VSampler;
class G4CellScorer;
class G4CellScorerStore;
class G4CellStoreScorer;

class B01MassScoring : public B01VSimulation {
public:
  B01MassScoring();
  virtual ~B01MassScoring();
  virtual const G4String &GetName() const;

  virtual G4VPhysicalVolume &GetMassGeometry();
  virtual const G4CellScorer *GetG4CellScorer();
  virtual void PrepareSampling();
  virtual void ConfigureSampling();
  virtual void SetWeightRoulette(G4bool wroulette);

  virtual void PostRun(G4std::ostream *);
private:
  B01MassScoring(const B01MassScoring &);
  B01MassScoring &operator=(const B01MassScoring &);

  G4String fName;
  G4bool fWeightRoulette;
  B01SlobedConcreteShield *fGeometry;
  const G4CellScorer *fSpecialCellScorer;
  G4CellScorerStore *fCS_store;
  G4CellStoreScorer *fScorer;
  G4VSampler *fSampler;
};

#endif
