#include "B01SimulationFactory.hh"
#include "B01MassScoring.hh"
#include "B01MassImportance.hh"
#include "B01MassImportanceScoring.hh"
#include "B01ParallelScoring.hh"
#include "B01ParallelImportance.hh"
#include "B01ParallelImportanceScoring.hh"
#include "B01NGapplication.hh"

B01SimulationFactory::B01SimulationFactory(){
  AddSimulation(new B01MassScoring);
  AddSimulation(new B01MassImportance);
  AddSimulation(new B01MassImportanceScoring);
  AddSimulation(new B01ParallelScoring);
  AddSimulation(new B01ParallelImportance);
  AddSimulation(new B01ParallelImportanceScoring);
  AddSimulation(new B01NGapplication);
}

B01SimulationFactory::~B01SimulationFactory(){
  for (B01MapNameSimulation::iterator it = fMapNameSimulation.begin();
       it != fMapNameSimulation.end(); it++) {
    B01SimStruct &simstruct = it->second;
    if (simstruct.fOwned) delete simstruct.fSimulation;
  }
}

void B01SimulationFactory::AddSimulation(B01VSimulation *sim) {
  B01SimStruct simstruct(sim);
  fMapNameSimulation[sim->GetName()] = simstruct;
}
  
G4String B01SimulationFactory::GetSimulationNames() const {

  G4String Names;
  std::string lbord(30,' ' );

  B01MapNameSimulation::const_iterator it;
  G4int i;
  for (it = fMapNameSimulation.begin(), i = 0;
       it != fMapNameSimulation.end(); it++, i++) {
    const B01SimStruct &simstruct = it->second;
    G4String text;
    G4String simname(simstruct.fSimulation->GetName());
    Names += lbord + simname + "\n";
  }
  Names += "\0";
  return Names;
}

G4bool B01SimulationFactory::
SimulationExists(const G4String &simname) const {
  B01MapNameSimulation::const_iterator it = 
    fMapNameSimulation.find(simname);
  if (it != fMapNameSimulation.end()) {
    return true;
  }
  return false;
}

B01VSimulation *B01SimulationFactory::Create(const G4String &simname){

  if (!SimulationExists(simname)) {
    G4cout << "B01SimulationFactory::Create: Simulation: " 
	   << simname << ", does not exist" << G4endl;
    return 0;
  }
  B01SimStruct &simstruct =  fMapNameSimulation[simname];
  simstruct.fOwned = false;
  return simstruct.fSimulation;
  
}
