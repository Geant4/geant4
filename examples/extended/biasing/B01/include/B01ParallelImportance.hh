#ifndef B01ParallelImportance_hh
#define B01ParallelImportance_hh B01ParallelImportance_hh

#include "B01VSimulation.hh"

class B01VGeometry;
class G4IStore;
class G4VSampler;
class B01ParallelGeometry;
class G4CellScorer;


class B01ParallelImportance : public B01VSimulation {
public:
  B01ParallelImportance();
  virtual ~B01ParallelImportance();
  virtual const G4String &GetName() const;

  virtual G4VPhysicalVolume &GetMassGeometry();
  virtual const G4CellScorer *GetG4CellScorer();
  virtual void PrepareSampling();
  virtual void ConfigureSampling();
  virtual void SetWeightRoulette(G4bool wroulette);

  virtual void PostRun(G4std::ostream *);
private:
  B01ParallelImportance(const B01ParallelImportance &);
  B01ParallelImportance &operator=
  (const B01ParallelImportance &);

  G4String fName;
  G4bool fWeightRoulette;
  B01VGeometry *fMassGeometry;
  B01ParallelGeometry *fParallelGeometry;
  G4IStore *fIStore;
  G4VSampler *fSampler;
};

#endif
