#include "G4ParallelStepper.hh"
#include "G4VPhysicalVolume.hh"

G4ParallelStepper::G4ParallelStepper(): fPStep(0) {}
G4ParallelStepper::~G4ParallelStepper(){
  delete fPStep;
}

G4ParallelStepper::G4ParallelStepper(const G4ParallelStepper &rhs){
  fPStep = new G4PStep(rhs.GetPStep());
}

G4ParallelStepper &G4ParallelStepper::operator=(const G4ParallelStepper &rhs){
  if (this != &rhs) {
    fPStep = new G4PStep(rhs.GetPStep());
  }
  return *this;
}


void G4ParallelStepper::Init(const G4PTouchableKey &aptk){
  if (!fPStep) {
    fPStep = new G4PStep(aptk, aptk);
  }
  else {
    fPStep->fPreTouchableKey = aptk;
    fPStep->fPostTouchableKey = aptk;
  }
  UnSetCrossBoundary();
}

void G4ParallelStepper::Update(const G4PTouchableKey &aptk){
  if (!fPStep) {
    Error("fPStep == 0, Init not called?");
  }
  fPStep->fPreTouchableKey = fPStep->fPostTouchableKey;
  fPStep->fPostTouchableKey = aptk;
  fPStep->fCrossBoundary = true;
}

void G4ParallelStepper::UnSetCrossBoundary(){
  if (!fPStep) {
    Error("fPStep == 0, Init not called?");
  }
  fPStep->fCrossBoundary = false;
}


