#include "B01ParallelScoring.hh"
#include "G4ScoreTable.hh"
#include "B01VGeometry.hh"
#include "B01ConcreteShield.hh"
#include "B01ParallelGeometry.hh"
#include "G4ParallelGeometrySampler.hh"


B01ParallelScoring::B01ParallelScoring()
  :
  fName("ParallelScoring"),
  fWeightRoulette(0),
  fMassGeometry(0),
  fParallelGeometry(0),
  fSpecialCellScorer(0),
  fSampler(0)
{}

B01ParallelScoring::~B01ParallelScoring(){
  if (fSampler) {
    delete fSampler;
  }
  if (fParallelGeometry) {
    delete fParallelGeometry;
  }
  if (fMassGeometry) {
    delete fMassGeometry;
  }
}

const G4String &B01ParallelScoring::GetName() const {
  return fName;
}

G4VPhysicalVolume &B01ParallelScoring::GetMassGeometry(){
  if (!fMassGeometry) {
    fMassGeometry = new B01ConcreteShield;
    if (!fMassGeometry) {
      G4std::G4Exception("B01ParallelScoring::GetMassGeometry: new failed to create B01ConcreteShield!");
    }
  }
  return fMassGeometry->GetWorldVolume();
  
}

const G4CellScorer *B01ParallelScoring::GetG4CellScorer(){
  return 0;  
}

void B01ParallelScoring::PrepareSampling(){
  GetMassGeometry();
  
  fParallelGeometry = new B01ParallelGeometry;
  if (!fParallelGeometry) {
    G4std::G4Exception("B01ParallelScoring::PrepareSampling: new failed to create B01ParallelGeometry!");
  }

}

void B01ParallelScoring::ConfigureSampling(){
  // create scoring sampler to sample neutrons
  // acording to the paralle geometry
  fSampler = new
    G4ParallelGeometrySampler(fParallelGeometry->GetWorldVolume(),
			      "neutron");
  if (!fSampler) {
    G4std::G4Exception("B01ParallelScoring::ConfigureSampling: new failed to create G4ParallelGeometrySampler!");
  }


  fSampler->PrepareScoring(&fScorer); 
  if (fWeightRoulette) {
    fSampler->PrepareWeightRoulett();
  }
  fSampler->Configure();
}

void B01ParallelScoring::SetWeightRoulette(G4bool wroulette){
  fWeightRoulette = wroulette;
}

void B01ParallelScoring::PostRun(G4std::ostream *out) {
  G4ScoreTable sp;
  sp.Print(fScorer.GetMapGeometryCellCellScorer(), out);
}
