#include "B01ParallelScoring.hh"
#include "G4RunManager.hh"
#include "G4Scorer.hh"
#include "G4ScoreTable.hh"
#include "B01VGeometry.hh"
#include "B01ConcreteShield.hh"
#include "B01ParallelGeometry.hh"
#include "G4CellScorerStore.hh"
#include "G4GeometryCell.hh"
#include "G4CellStoreScorer.hh"
#include "G4ParallelScoreSampler.hh"
#include "B01Run.hh"
#include "G4Pstring.hh"

B01ParallelScoring::B01ParallelScoring()
  :
  fName("ParallelScoring"),
  fMassGeometry(0),
  fRun(0),
  fParallelGeometry(0),
  fCS_store(0),
  fScorer(0),
  fSampler(0)
{}

B01ParallelScoring::~B01ParallelScoring(){
  if (fSampler) delete fSampler;
  if (fScorer) delete fScorer;
  if (fCS_store) delete fCS_store;
  if (fParallelGeometry) delete fParallelGeometry;
  if (fMassGeometry) delete fMassGeometry;
  if (fRun) delete fRun;
}

G4String B01ParallelScoring::GetName() const {
  return fName;
}

void B01ParallelScoring::Construct() {
  
  fMassGeometry = new B01ConcreteShield;
  fRun = new B01Run;
  fRun->SetDetector(fMassGeometry->GetWorldVolume());
  fRun->Initialize();
  
  fParallelGeometry = new B01ParallelGeometry;

  
  // create an cell scorer stroe and have a cell scorer
  // created for given geometry cells
  fCS_store = new G4CellScorerStore();

  G4int i;
  for (i=1; i <= 18; i++) {
    G4String volname = fParallelGeometry->GetCellName(i);
    
    const G4VPhysicalVolume *pvol = fParallelGeometry->
      GetPhysicalVolumeByName(volname);
    if (pvol) {
      G4GeometryCell gCell(*pvol, 0);
      fCS_store->AddCellScorer(gCell);
    }
  }
  
  // create a scorer
  fScorer = new G4CellStoreScorer(*fCS_store);


  // create scoring sampler to sample neutrons
  // acording to the paralle geometry
  fSampler = new
    G4ParallelScoreSampler(fParallelGeometry->GetWorldVolume(),
			   "neutron",
			   *fScorer);
  fSampler->Initialize();
 
}

void B01ParallelScoring::Run(G4int nevents) {
  if (!fRun) {
    G4cout << "B01ParallelScoring::Run: no BooRun constructed yet!" << G4endl;
  }
  else {
    fRun->BeamOn(nevents);
  }
}

void B01ParallelScoring::PostRun(G4std::ostream *out) {
  G4ScoreTable sp;
  sp.Print(fCS_store->GetMapGeometryCellCellScorer(), out);
}
