#include "B01SimulationFactory.hh"
#include "B01MassScoring.hh"
#include "B01MassImportance.hh"
#include "B01MassImportanceScoring.hh"
#include "B01ParallelScoring.hh"
#include "B01ParallelImportance.hh"
#include "B01ParallelImportanceScoring.hh"

B01SimulationFactory::B01SimulationFactory(){
  AddSimulation(new B01MassScoring);
  AddSimulation(new B01MassImportance);
  AddSimulation(new B01MassImportanceScoring);
  AddSimulation(new B01ParallelScoring);
  AddSimulation(new B01ParallelImportance);
  AddSimulation(new B01ParallelImportanceScoring);
}

B01SimulationFactory::~B01SimulationFactory(){
  for (B01MapNameSimulation::iterator it = fMapNameSimulation.begin();
       it != fMapNameSimulation.end(); ++it) {
    B01SimStruct &simstruct = it->second;
    if (simstruct.fOwned) {
      delete simstruct.fSimulation;
    }
  }
}

void B01SimulationFactory::AddSimulation(B01VSimulation *sim) {
  if  (!sim) {
    G4std::G4Exception("ERROR:B01SimulationFactory::AddSimulation: called with pointer to null!");
  }
  B01SimStruct simstruct = {sim, false};
  fMapNameSimulation[sim->GetName()] = simstruct;
}
  
B01SimNameVec B01SimulationFactory::GetSimulationNames() const {

  B01SimNameVec simnamevec;
  B01MapNameSimulation::const_iterator it;
  for (it = fMapNameSimulation.begin();
       it != fMapNameSimulation.end(); ++it) {
    const B01SimStruct &simstruct = it->second;
    G4String simname(simstruct.fSimulation->GetName());
    simnamevec.push_back(simname);
  }
  return simnamevec;
}

G4bool B01SimulationFactory::
SimulationExists(const G4String &simname) const {
  G4bool sim = false;
  B01MapNameSimulation::const_iterator it = 
    fMapNameSimulation.find(simname);
  if (it != fMapNameSimulation.end()) {
    sim =  true;
  }
  return sim;
}

B01VSimulation *B01SimulationFactory::Create(const G4String &simname){

  B01VSimulation *psim = 0;
  if (!SimulationExists(simname)) {
    G4std::G4cout << "B01SimulationFactory::Create: Simulation: " 
	   << simname << ", does not exist" << G4endl;
  }
  else {
    B01SimStruct &simstruct =  fMapNameSimulation[simname];
    simstruct.fOwned = false;
    psim =  simstruct.fSimulation;
  }
  return psim;
}
