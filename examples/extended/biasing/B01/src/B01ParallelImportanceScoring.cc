#include "B01ParallelImportanceScoring.hh"
#include "G4RunManager.hh"
#include "G4Scorer.hh"
#include "G4ScoreTable.hh"
#include "B01VGeometry.hh"
#include "B01ConcreteShield.hh"
#include "B01ParallelGeometry.hh"
#include "G4IStore.hh"
#include "G4CellScorerStore.hh"
#include "G4GeometryCell.hh"
#include "G4CellStoreScorer.hh"
#include "G4ParallelImportanceScoreSampler.hh"
#include "B01Run.hh"
#include "G4Pstring.hh"

B01ParallelImportanceScoring::B01ParallelImportanceScoring()
  :
  fName("ParallelImportanceScoring"),
  fMassGeometry(0),
  fRun(0),
  fParallelGeometry(0),
  fIStore(0),
  fCS_store(0),
  fScorer(0),
  fSampler(0)
{}

B01ParallelImportanceScoring::~B01ParallelImportanceScoring(){
  if (fSampler) delete fSampler;
  if (fScorer) delete fScorer;
  if (fCS_store) delete fCS_store;
  if (fIStore) delete fIStore;
  if (fParallelGeometry) delete fParallelGeometry;
  if (fMassGeometry) delete fMassGeometry;
  if (fRun) delete fRun;
}

G4String B01ParallelImportanceScoring::GetName() const {
  return fName;
}

void B01ParallelImportanceScoring::Construct() {
  
  fMassGeometry = new B01ConcreteShield;
  fRun = new B01Run;
  fRun->SetDetector(fMassGeometry->GetWorldVolume());
  fRun->Initialize();
  
  fParallelGeometry = new B01ParallelGeometry;

  
  // create an importance store and fill it with the importance
  // per cell values
  // and create an cell scorer stroe and have a cell scorer
  // created for given geometry cells
  fIStore = new G4IStore(fParallelGeometry->GetWorldVolume());
  fCS_store = new G4CellScorerStore();

  G4int i;
  for (i=1; i <= 18; i++) {
    G4String volname = fParallelGeometry->GetCellName(i);
    G4double imp = pow(2,i-1);

    const G4VPhysicalVolume *pvol = fParallelGeometry->
      GetPhysicalVolumeByName(volname);
    if (pvol) {
      G4GeometryCell gCell(*pvol, 0);
      fIStore->AddImportanceGeometryCell(imp, gCell);
      fCS_store->AddCellScorer(gCell);
    }
  }
  
  // create a scorer
  fScorer = new G4CellStoreScorer(*fCS_store);


  // create importance sampler to importance sample neutrons
  // in acording to the paralle geometry
  fSampler = new
    G4ParallelImportanceScoreSampler(fParallelGeometry->GetWorldVolume(),
				     *fIStore, 
				     *fScorer, 
				     "neutron");
  fSampler->Initialize();
 
}

void B01ParallelImportanceScoring::Run(G4int nevents) {
  if (!fRun) {
    G4cout << "B01ParallelImportanceScoring::Run: no BooRun constructed yet!" << G4endl;
  }
  else {
    fRun->BeamOn(nevents);
  }
}

void B01ParallelImportanceScoring::PostRun(G4std::ostream *out) {
  G4ScoreTable sp(fIStore);
  sp.Print(fCS_store->GetMapGeometryCellCellScorer(), out);
}
