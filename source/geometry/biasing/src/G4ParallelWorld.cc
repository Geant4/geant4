#include "G4ParallelWorld.hh"
#include "G4ParallelStepper.hh"
#include "G4ParallelNavigator.hh"
#include "G4VPhysicalVolume.hh"

G4ParallelWorld::G4ParallelWorld(G4VPhysicalVolume &worldvolume):
  fWorldVolume(&worldvolume),
  fPstepper(new G4ParallelStepper),
  fPdriver(new G4ParallelNavigator(*fWorldVolume))
{}

G4ParallelWorld::~G4ParallelWorld(){
  delete fPdriver;
  delete fPstepper;
}

G4ParallelWorld::G4ParallelWorld(const G4ParallelWorld &rhs):
  fWorldVolume(rhs.GetWorldVolume()){
  fPstepper = new G4ParallelStepper;
  fPdriver = new G4ParallelNavigator(*fWorldVolume);
}

G4ParallelWorld &G4ParallelWorld::operator=(const G4ParallelWorld &rhs){
  if (this!=&rhs) {

    fWorldVolume = rhs.GetWorldVolume();
    fPstepper = new G4ParallelStepper;
    fPdriver = new G4ParallelNavigator(*fWorldVolume);
  }
  return *this;
}


G4VPhysicalVolume *
G4ParallelWorld::GetWorldVolume() const {
  return fWorldVolume;
}


G4VParallelStepper &G4ParallelWorld::GetParallelStepper(){
  return *fPstepper;
}

G4VPGeoDriver &G4ParallelWorld::GetGeoDriver(){
  return *fPdriver;
}
