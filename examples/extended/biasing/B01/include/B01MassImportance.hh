#ifndef B01MassImportance_hh
#define B01MassImportance_hh B01MassImportance_hh

#include "B01VSimulation.hh"

class B01SlobedConcreteShield;
class B01Run;
class G4MassImportanceSampler;
class G4IStore;

class B01MassImportance : public B01VSimulation {
public:
  B01MassImportance();
  ~B01MassImportance();
  G4String GetName() const;
  void Construct();
  void Run(G4int nevents);
  void PostRun(G4std::ostream *);
private:
  B01MassImportance(const B01MassImportance &);
  B01MassImportance &operator=(const B01MassImportance &);

  G4String fName;
  B01SlobedConcreteShield *fGeometry;
  B01Run *fRun;
  G4IStore *fIStore;
  G4MassImportanceSampler *fSampler;
};

#endif
