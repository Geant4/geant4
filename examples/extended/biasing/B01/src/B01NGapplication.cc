#include "B01NGapplication.hh"
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
#include "G4ParallelScoreSampler.hh"
#include "G4WeightWindowAlgorithm.hh"
#include "G4PImportanceWWindowScoreSampler.hh"
#include "B01Run.hh"
#include "G4Pstring.hh"
#include "G4UImanager.hh"

B01NGapplication::B01NGapplication()
  :
  fName("NGapplication"),
  fMassGeometry(0),
  fRun(0),
  fParallelGeometry(0),
  fIStore(0),
  fCS_store_n(0),
  fCS_store_g(0),
  fScorer_n(0),
  fScorer_g(0),
  fWWalg(0),
  fSampler_n(0),
  fSampler_g(0)
{}

B01NGapplication::~B01NGapplication(){
  if (fSampler_n) delete fSampler_n;
  if (fSampler_g) delete fSampler_g;
  if (fWWalg) delete fWWalg;
  if (fScorer_n) delete fScorer_n;
  if (fScorer_g) delete fScorer_g;
  if (fCS_store_n) delete fCS_store_n;
  if (fCS_store_g) delete fCS_store_g;
  if (fIStore) delete fIStore;
  if (fParallelGeometry) delete fParallelGeometry;
  if (fMassGeometry) delete fMassGeometry;
  if (fRun) delete fRun;
}

G4String B01NGapplication::GetName() const {
  return fName;
}

void B01NGapplication::Construct() {
  
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
  fCS_store_n = new G4CellScorerStore();
  fCS_store_g = new G4CellScorerStore();

  G4int i;
  for (i=1; i <= 18; i++) {
    G4String volname = fParallelGeometry->GetCellName(i);
    G4double imp = pow(1.2,i-1);

    const G4VPhysicalVolume *pvol = fParallelGeometry->
      GetPhysicalVolumeByName(volname);
    if (pvol) {
      G4GeometryCell gCell(*pvol, 0);
      fIStore->AddImportanceGeometryCell(imp, gCell);
      fCS_store_n->AddCellScorer(gCell);
      fCS_store_g->AddCellScorer(gCell);     
    }
  }
  G4GeometryCell gCell(fParallelGeometry->GetWorldVolume(), -1);
  fIStore->AddImportanceGeometryCell(1, gCell);

  
  // create a scorer
  fScorer_n = new G4CellStoreScorer(*fCS_store_n);
  fScorer_g = new G4CellStoreScorer(*fCS_store_g);


  // create importance sampler to importance sample neutrons
  // in acording to the paralle geometry

  

  fWWalg = new G4WeightWindowAlgorithm;
  fWWalg->SetUpperLimit(5);
  fWWalg->SetLowerLimit(0.2);
  
  // create the sampler for importance nad weight window sampling 
  // of neutron
  fSampler_n = new G4PImportanceWWindowScoreSampler(fParallelGeometry->
						    GetWorldVolume(),
						    *fIStore, 
						    *fScorer_n, 
						    "neutron",
						    *fWWalg);
  
  fSampler_g = new
    G4ParallelScoreSampler(fParallelGeometry->GetWorldVolume(),
			   "gamma",
			   *fScorer_g);
  fSampler_n->Initialize();
  fSampler_g->Initialize();

  G4UImanager *UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/gun/energy 300.0 MeV"); 
  
}

void B01NGapplication::Run(G4int nevents) {
  if (!fRun) {
    G4cout << "B01NGapplication::Run: no BooRun constructed yet!" << G4endl;
  }
  else {
    fRun->BeamOn(nevents);
  }
}

void B01NGapplication::PostRun(G4std::ostream *out) {
  *out << "Results for neutrons: " << G4endl;
  G4ScoreTable spn(fIStore);
  spn.Print(fCS_store_n->GetMapGeometryCellCellScorer(), out);

  *out << "\n" << "Results for gammas:" << G4endl;
  G4ScoreTable spg;
  spg.Print(fCS_store_g->GetMapGeometryCellCellScorer(), out);
  
}




