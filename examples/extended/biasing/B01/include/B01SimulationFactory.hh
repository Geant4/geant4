#ifndef B01SimulationFactory_hh
#define B01SimulationFactory_hh B01SimulationFactory_hh

#include "globals.hh"
#include "g4std/vector"
#include "g4std/map"

class B01VSimulation;

struct B01SimStruct{
  B01SimStruct(B01VSimulation *sim = 0) 
    : 
    fSimulation(sim),
    fOwned(true)
{}
  B01VSimulation *fSimulation;
  G4bool fOwned;
};


typedef G4std::map<G4String , B01SimStruct> B01MapNameSimulation;



class B01SimulationFactory {
public:
  B01SimulationFactory();
  ~B01SimulationFactory();
  
  G4String GetSimulationNames() const;
  G4bool SimulationExists(const G4String &simname) const;
  B01VSimulation *Create(const G4String &simname);

private:
  void AddSimulation(B01VSimulation *sim);
  B01MapNameSimulation fMapNameSimulation;

};

#endif
