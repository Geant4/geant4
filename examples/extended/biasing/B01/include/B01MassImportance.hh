#ifndef B01MassImportance_hh
#define B01MassImportance_hh B01MassImportance_hh

#include "B01VSimulation.hh"

class B01SlobedConcreteShield;
class G4CellScorer;
class G4VSampler;
class G4IStore;

class B01MassImportance : public B01VSimulation {
public:
  B01MassImportance();
  virtual ~B01MassImportance();
  virtual const G4String &GetName() const;
  virtual G4VPhysicalVolume &GetMassGeometry();
  virtual const G4CellScorer *GetG4CellScorer();
  virtual void PrepareSampling();
  virtual void ConfigureSampling();
  virtual void SetWeightRoulette(G4bool wroulette);
  virtual void PostRun(G4std::ostream *);
private:
  B01MassImportance(const B01MassImportance &);
  B01MassImportance &operator=(const B01MassImportance &);

  G4String fName;
  G4bool fWeightRoulette;
  B01SlobedConcreteShield *fGeometry;
  G4IStore *fIStore;
  G4VSampler *fSampler;
};

#endif
