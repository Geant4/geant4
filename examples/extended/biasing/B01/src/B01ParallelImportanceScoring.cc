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
#include "G4ParallelGeometrySampler.hh"


B01ParallelImportanceScoring::B01ParallelImportanceScoring()
  :
  fName("ParallelImportanceScoring"),
  fWeightRoulette(0),
  fMassGeometry(0),
  fParallelGeometry(0),
  fSpecialCellScorer(0),
  fIStore(0),
  fCS_store(0),
  fScorer(0),
  fSampler(0)
{}

B01ParallelImportanceScoring::~B01ParallelImportanceScoring(){
  if (fSampler) {
    delete fSampler;
  }
  if (fScorer) {
    delete fScorer;
  }
  if (fCS_store) {
    delete fCS_store;
  }
  if (fIStore) {
    delete fIStore;
  }
  if (fParallelGeometry) {
    delete fParallelGeometry;
  }
  if (fMassGeometry) {
    delete fMassGeometry;
  }
}


const G4String &B01ParallelImportanceScoring::GetName() const {
  return fName;
}

G4VPhysicalVolume &B01ParallelImportanceScoring::GetMassGeometry(){
  if (!fMassGeometry) {
    fMassGeometry = new B01ConcreteShield;
    if (!fMassGeometry) {
      G4std::G4Exception("B01ParallelImportanceScoring::GetMassGeometry: new failes to create B01ConcreteShield!");
    }
  }
  return fMassGeometry->GetWorldVolume();
}

const G4CellScorer *B01ParallelImportanceScoring::GetG4CellScorer(){
  if (!fSpecialCellScorer) {
    G4std::G4Exception("B01ParallelImportanceScoring::GetG4CellScorer: Error");
  }
  return fSpecialCellScorer;
}

void B01ParallelImportanceScoring::PrepareSampling(){
  GetMassGeometry();
  fParallelGeometry = new B01ParallelGeometry;
  if (!fParallelGeometry) {
    G4std::G4Exception("B01ParallelImportanceScoring::PrepareSamplingnew failed to create B01ParallelGeometry!");
  }
  
  // create an importance store and fill it with the importance
  // per cell values
  // and create an cell scorer stroe and have a cell scorer
  // created for given geometry cells
  const G4VPhysicalVolume &pworld = fParallelGeometry->GetWorldVolume();
  fIStore = new G4IStore(pworld);
  if (!fIStore) {
    G4std::G4Exception("B01ParallelImportanceScoring::PrepareSamplingnew failed to create G4IStore!");
  }
  // adding GeometryCell for world volume. ReplicaNumer = -1 !
  G4GeometryCell gWorldCell(pworld, -1);
  fIStore->AddImportanceGeometryCell(1, gWorldCell);

  fCS_store = new G4CellScorerStore();
  if (!fCS_store) {
    G4std::G4Exception("B01ParallelImportanceScoring::PrepareSamplingnew failed to create G4CellScorerStore!");
  }
  fCS_store->AddCellScorer(gWorldCell);

  G4int i = 1;
  for (i=1; i <= 19; ++i) {
    G4String volname = fParallelGeometry->GetCellName(i);
    G4double imp = pow(2,i-1);
    if (i==19) {
      imp = pow(2,17);
    }
    const G4VPhysicalVolume *pvol = fParallelGeometry->
      GetPhysicalVolumeByName(volname);
    if (pvol) {
      G4GeometryCell gCell(*pvol, 0);
      fIStore->AddImportanceGeometryCell(imp, gCell);
      const G4CellScorer *s = fCS_store->AddCellScorer(gCell);
      if (i==18) {
	fSpecialCellScorer = s;
      }
    }
  }
  
  // create a scorer
  fScorer = new G4CellStoreScorer(*fCS_store);
  if (!fScorer) {
    G4std::G4Exception("B01ParallelImportanceScoring::PrepareSamplingnew failed to create G4CellStoreScorer!");
  }
}

void B01ParallelImportanceScoring::ConfigureSampling(){
  fSampler = new
    G4ParallelGeometrySampler(fParallelGeometry->GetWorldVolume(),
			      "neutron");
  if (!fSampler) {
    G4std::G4Exception("B01ParallelImportanceScoring::ConfigureSampling: new failed to create G4ParallelGeometrySampler!");
  }

  fSampler->PrepareImportanceSampling(fIStore);
  fSampler->PrepareScoring(fScorer);
  if (fWeightRoulette) {
    fSampler->PrepareWeightRoulett();
  }
  fSampler->Configure();
}

void B01ParallelImportanceScoring::SetWeightRoulette(G4bool wroulette){
  fWeightRoulette = wroulette;
}


void B01ParallelImportanceScoring::PostRun(G4std::ostream *out) {
  G4ScoreTable sp(fIStore);
  sp.Print(fCS_store->GetMapGeometryCellCellScorer(), out);
}
