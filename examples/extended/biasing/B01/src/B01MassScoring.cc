#include "B01MassScoring.hh"
#include "G4MassGeometrySampler.hh"
#include "G4CellScorerStore.hh"
#include "G4CellStoreScorer.hh"
#include "G4ScoreTable.hh"
#include "B01VGeometry.hh"
#include "B01SlobedConcreteShield.hh"


B01MassScoring::B01MassScoring()
  :
  fName("MassScoring"),
  fWeightRoulette(false),
  fGeometry(0),
  fSpecialCellScorer(0),
  fCS_store(0),
  fScorer(0),
  fSampler(0)
{}

B01MassScoring::~B01MassScoring(){
  if (fSampler) {
    delete fSampler;
  }
  if (fGeometry) {
    delete fGeometry;
  }
}

const G4String &B01MassScoring::GetName() const {
  return fName;
}

G4VPhysicalVolume &B01MassScoring::GetMassGeometry(){
  if (!fGeometry) {
    fGeometry = new B01SlobedConcreteShield;
    if (!fGeometry) {
      G4std::G4Exception("B01MassScoring::GetMassGeometry: new failed to create B01SlobedConcreteShield!");
    }
  }
  return fGeometry->GetWorldVolume();
}

const G4CellScorer *B01MassScoring::GetG4CellScorer(){
  if (!fSpecialCellScorer) {
    G4std::G4Exception("B01MassScoring::GetG4CellScorer: Error");
  }
  return fSpecialCellScorer;

}

void B01MassScoring::PrepareSampling(){
  GetMassGeometry();
  
  G4GeometryCell gWorldCell(fGeometry->GetWorldVolume(), -1);
  
  fCS_store = new G4CellScorerStore();
  if (!fCS_store) {
    G4std::G4Exception("B01MassScoring::PrepareSampling: new failed to create G4CellScorerStore!");
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
      const G4CellScorer *s = fCS_store->AddCellScorer(gCell);
      if (i==18) {
	fSpecialCellScorer = s;
      }
    }
  }
  fScorer = new G4CellStoreScorer(*fCS_store);
  if (!fScorer) {
    G4std::G4Exception("B01MassScoring::PrepareSampling: new failed to create G4CellStoreScorer!");
  }
}

void B01MassScoring::ConfigureSampling(){
  fSampler = new G4MassGeometrySampler("neutron"); 
  if (!fSampler) {
    G4std::G4Exception("B01MassScoring::ConfigureSampling new failed to create G4MassGeometrySampler!");
  }
  

  fSampler->PrepareScoring(fScorer); 
  if (fWeightRoulette) {
    fSampler->PrepareWeightRoulett();  
  }
  fSampler->Configure();

}

void B01MassScoring::SetWeightRoulette(G4bool wroulette){
  fWeightRoulette = wroulette;
}


void B01MassScoring::PostRun(G4std::ostream *out) {
  G4ScoreTable sp;
  sp.Print(fCS_store->GetMapGeometryCellCellScorer(), out);
}

