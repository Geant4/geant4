#ifndef B01ParallelScoring_hh
#define B01ParallelScoring_hh B01ParallelScoring_hh

#include "B01VSimulation.hh"
#include "G4Scorer.hh"
class B01VGeometry;
class G4VSampler;
class B01ParallelGeometry;
class G4CellScorer;


class B01ParallelScoring : public B01VSimulation {
public:
  B01ParallelScoring();
  virtual ~B01ParallelScoring();
  virtual const G4String &GetName() const;

  virtual G4VPhysicalVolume &GetMassGeometry();
  virtual const G4CellScorer *GetG4CellScorer();
  virtual void PrepareSampling();
  virtual void ConfigureSampling();
  virtual void SetWeightRoulette(G4bool wroulette);

  virtual void PostRun(G4std::ostream *);
private:
  B01ParallelScoring(const B01ParallelScoring &);
  B01ParallelScoring &operator=
  (const B01ParallelScoring &);

  G4String fName;
  G4bool fWeightRoulette;
  B01VGeometry *fMassGeometry;
  B01ParallelGeometry *fParallelGeometry;
  const G4CellScorer *fSpecialCellScorer;
  G4Scorer fScorer;
  G4VSampler *fSampler;
};

#endif
