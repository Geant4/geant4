#include "G4ParallelManager.hh"
#include "G4ParallelStepper.hh"
#include "G4ParallelWorld.hh"
#include "G4ParallelNavigator.hh"
#include "G4ParallelTransport.hh"
#include "G4ProcessPlacer.hh"
#include "G4VPhysicalVolume.hh"

G4ParallelManager::G4ParallelManager(G4VPhysicalVolume &worldvolume,
				     const G4String &particlename) :
  fPworld(new G4ParallelWorld(worldvolume)),
  fParticleName(particlename),
  fParallelTransport(0)
{}

  
G4ParallelManager::~G4ParallelManager(){
  delete fPworld;
  if (fParallelTransport) delete fParallelTransport;
}

G4ParallelWorld &G4ParallelManager::GetParallelWorld(){
  return *fPworld;
}

G4String G4ParallelManager::GetParticleName(){
  return fParticleName;
}


G4ParallelTransport *G4ParallelManager::CreateParallelTransport(){
  if (!fParallelTransport) {
    fParallelTransport = new G4ParallelTransport(fPworld->GetGeoDriver(), 
				       fPworld->GetParallelStepper());
  }
  return fParallelTransport;
}

void G4ParallelManager::Initialize(){
  G4ProcessPlacer placer(fParticleName);
  placer.AddProcessAsSecondDoIt(CreateParallelTransport());
}
