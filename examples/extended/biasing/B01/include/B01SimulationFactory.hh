#ifndef B01SimulationFactory_hh
#define B01SimulationFactory_hh B01SimulationFactory_hh

#include "globals.hh"
#include "g4std/vector"
#include "g4std/map"

class B01VSimulation;

struct B01SimStruct{
  B01VSimulation *fSimulation;
  G4bool fOwned;
};


typedef G4std::map<G4String , B01SimStruct> B01MapNameSimulation;
typedef G4std::vector<G4String> B01SimNameVec;


class B01SimulationFactory {
public:
  B01SimulationFactory();
  ~B01SimulationFactory();
  
  B01SimNameVec GetSimulationNames() const;
  G4bool SimulationExists(const G4String &simname) const;
  B01VSimulation *Create(const G4String &simname);

private:
  B01SimulationFactory(const B01SimulationFactory &);

  void AddSimulation(B01VSimulation *sim);

  B01SimulationFactory &operator=(const B01SimulationFactory &);

  B01MapNameSimulation fMapNameSimulation;

};

#endif
