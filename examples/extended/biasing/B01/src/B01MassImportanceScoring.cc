#include "B01MassImportanceScoring.hh"
#include "G4RunManager.hh"
#include "G4MassImportanceScoreSampler.hh"
#include "G4CellScorerStore.hh"
#include "G4CellStoreScorer.hh"
#include "G4ScoreTable.hh"
#include "G4IStore.hh"
#include "B01VGeometry.hh"
#include "B01SlobedConcreteShield.hh"
#include "B01Run.hh"
#include "G4Pstring.hh"

B01MassImportanceScoring::B01MassImportanceScoring()
  :
  fName("MassImportanceScoring"),
  fGeometry(0),
  fRun(0),
  fIStore(0),
  fCS_store(0),
  fScorer(0),
  fSampler(0)
{}

B01MassImportanceScoring::~B01MassImportanceScoring(){
  if (fSampler) delete fSampler;
  if (fCS_store) delete fCS_store;
  if (fScorer) delete fScorer;
  if (fIStore) delete fIStore;
  if (fGeometry) delete fGeometry;
  if (fRun) delete fRun;
}

G4String B01MassImportanceScoring::GetName() const {
  return fName;
}
void B01MassImportanceScoring::Construct() {
  
  fGeometry = new B01SlobedConcreteShield;
  fRun = new B01Run;
  fRun->SetDetector(fGeometry->GetWorldVolume());
  fRun->Initialize();

  // create an importance store and fill it with the importance
  // per cell values
  // and create an cell scorer stroe and have a cell scorer
  // created for given geometry cells
  fIStore = new G4IStore(fGeometry->GetWorldVolume());
  fCS_store = new G4CellScorerStore();

  G4int i;
  for (i=1; i <= 18; i++) {
    G4String volname = fGeometry->GetCellName(i);
    G4double imp = pow(2,i-1);
    
    const G4VPhysicalVolume *pvol = fGeometry->
      GetPhysicalVolumeByName(volname);
    if (pvol) {
      G4GeometryCell gCell(*pvol, 0);
      fIStore->AddImportanceGeometryCell(imp, gCell);
      fCS_store->AddCellScorer(gCell);
    }
  }
  
  // create a scorer
  fScorer = new G4CellStoreScorer(*fCS_store);

  fSampler = new G4MassImportanceScoreSampler(*fIStore,
					      *fScorer, 
					      "neutron"); 
  
  // to be done after fRun->Initialize()
  fSampler->Initialize(); 


}

void B01MassImportanceScoring::Run(G4int nevents) {
  if (!fRun) {
    G4cout << "B01MassImportanceScoring::Run: no BooRun constructed yet!" << G4endl;
  }
  else {
    fRun->BeamOn(nevents);
  }
}

void B01MassImportanceScoring::PostRun(G4std::ostream *out) {
  G4ScoreTable sp(fIStore);
  sp.Print(fCS_store->GetMapGeometryCellCellScorer(), out);
}
