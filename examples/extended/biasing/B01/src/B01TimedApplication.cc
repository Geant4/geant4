#include "B01TimedApplication.hh"
#include "B01VSimulation.hh"


B01TimedApplication::B01TimedApplication(G4int time)
  :fTimedRun(time)
{}
B01TimedApplication::~B01TimedApplication()
{}
void B01TimedApplication::RunSimulation(B01VSimulation *sim){
  sim->PrepareSampling();

  fTimedRun.SetDetector(sim->GetMassGeometry());
  fTimedRun.SetSpecialG4CellScorer(sim->GetG4CellScorer());
  fTimedRun.Initialize();
  
  sim->ConfigureSampling();
  fTimedRun.BeamOn(100000000);
  
  sim->PostRun(&G4std::G4cout);
  
}
