#include "B01MassImportanceScoring.hh"
#include "G4MassGeometrySampler.hh"
#include "G4CellScorerStore.hh"
#include "G4CellStoreScorer.hh"
#include "G4ScoreTable.hh"
#include "G4IStore.hh"
#include "B01VGeometry.hh"
#include "B01SlobedConcreteShield.hh"

#include "B01TimedEventAction.hh"



B01MassImportanceScoring::B01MassImportanceScoring()
  :
  fName("MassImportanceScoring"),
  fWeightRoulette(false),
  fGeometry(0),
  fSpecialCellScorer(0),
  fIStore(0),
  fCS_store(0),
  fScorer(0),
  fSampler(0)
{}

B01MassImportanceScoring::~B01MassImportanceScoring(){
  if (fSampler) {
    delete fSampler;
  }
  if (fCS_store) {
    delete fCS_store;
  }
  if (fScorer) {
    delete fScorer;
  }
  if (fIStore) {
    delete fIStore;
  }
  if (fGeometry) {
    delete fGeometry;
  }
}

const G4String &B01MassImportanceScoring::GetName() const {
  return fName;
}

void B01MassImportanceScoring::SetWeightRoulette(G4bool wroulette){
  fWeightRoulette = wroulette;
}


G4VPhysicalVolume &B01MassImportanceScoring::GetMassGeometry(){
  if (!fGeometry) {
    fGeometry = new B01SlobedConcreteShield;
    if (!fGeometry) {
      G4std::G4Exception("B01MassImportanceScoring::GetMassGeometry: new failed to create B01SlobedConcreteShield!");
    }
  }
  return fGeometry->GetWorldVolume();
}

const G4CellScorer *B01MassImportanceScoring::GetG4CellScorer(){
  if (!fSpecialCellScorer) {
    G4std::G4Exception("B01MassImportanceScoring::GetG4CellScorer: Error");
  }
  return fSpecialCellScorer;
}
void B01MassImportanceScoring::PrepareSampling(){
  GetMassGeometry();
  fIStore = new G4IStore(fGeometry->GetWorldVolume());
  if (!fIStore) {
    G4std::G4Exception("B01MassImportanceScoring::PrepareSampling: new failed to create G4IStore!");
  }
  G4GeometryCell gWorldCell(fGeometry->GetWorldVolume(), -1);
  fIStore->AddImportanceGeometryCell(1, gWorldCell);
  
  fCS_store = new G4CellScorerStore();
  if (!fCS_store) {
    G4std::G4Exception("B01MassImportanceScoring::PrepareSampling: new failed to create G4CellScorerStore!");
  }
  fCS_store->AddCellScorer(gWorldCell);


  G4int i = 1;
  for (i=1; i <= 19; i++) {
    G4String volname = fGeometry->GetCellName(i);
    G4double imp = pow(2,i-1);
    if (i==19) {
      imp = pow(2,17);
    }
    
    const G4VPhysicalVolume *pvol = fGeometry->
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
  fScorer = new G4CellStoreScorer(*fCS_store);
  if (!fScorer) {
    G4std::G4Exception("B01MassImportanceScoring::PrepareSampling: new failed to create G4CellStoreScorer!");
  }

}
void B01MassImportanceScoring::ConfigureSampling(){

  fSampler = new G4MassGeometrySampler("neutron"); 
  if  (!fSampler) {
    G4std::G4Exception("B01MassImportanceScoring::ConfigureSampling: new failed to create G4MassGeometrySampler!");
  }
  fSampler->PrepareImportanceSampling(fIStore);
  fSampler->PrepareScoring(fScorer); 
  if (fWeightRoulette) {
    fSampler->PrepareWeightRoulett();
  }
  fSampler->Configure();
}

void B01MassImportanceScoring::PostRun(G4std::ostream *out) {
  G4ScoreTable sp(fIStore);
  sp.Print(fCS_store->GetMapGeometryCellCellScorer(), out);
}
