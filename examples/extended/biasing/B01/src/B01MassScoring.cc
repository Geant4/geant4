#include "B01MassScoring.hh"
#include "G4RunManager.hh"
#include "G4MassScoreSampler.hh"
#include "G4Scorer.hh"
#include "G4ScoreTable.hh"
#include "B01VGeometry.hh"
#include "B01SlobedConcreteShield.hh"
#include "B01Run.hh"

B01MassScoring::B01MassScoring()
  :
  fName("MassScoring"),
  fGeometry(0),
  fRun(0),
  fScorer(0),
  fSampler(0)
{}

B01MassScoring::~B01MassScoring(){
  if (fSampler) delete fSampler;
  if (fScorer) delete fScorer;
  if (fGeometry) delete fGeometry;
  if (fRun) delete fRun;
}

G4String B01MassScoring::GetName() const {
  return fName;
}
void B01MassScoring::Construct() {
  
  fGeometry = new B01SlobedConcreteShield;
  fRun = new B01Run;
  fRun->SetDetector(fGeometry->GetWorldVolume());
  fRun->Initialize();

  fScorer = new G4Scorer; // a scorer 
  fSampler = new G4MassScoreSampler(*fScorer, "neutron"); 
  
  // to be done after fRun->Initialize()
  fSampler->Initialize(); 


}

void B01MassScoring::Run(G4int nevents) {
  if (!fRun) {
    G4cout << "B01MassScoring::Run: no BooRun constructed yet!" << G4endl;
  }
  else {
    fRun->BeamOn(nevents);
  }
}

void B01MassScoring::PostRun(G4std::ostream *out) {
  G4ScoreTable sp;
  sp.Print(fScorer->GetMapGeometryCellCellScorer(), out);
}

