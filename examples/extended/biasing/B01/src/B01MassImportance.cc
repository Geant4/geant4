#include "B01MassImportance.hh"
#include "G4MassGeometrySampler.hh"
#include "G4IStore.hh"
#include "B01VGeometry.hh"
#include "B01SlobedConcreteShield.hh"


B01MassImportance::B01MassImportance()
  :
  fName("MassImportance"),
  fWeightRoulette(false),
  fGeometry(0),
  fIStore(0),
  fSampler(0)
{}

B01MassImportance::~B01MassImportance(){
  if (fSampler) {
    delete fSampler;
  }
  if (fIStore) {
    delete fIStore;
  }
  if (fGeometry) {
    delete fGeometry;
  }
}

const G4String &B01MassImportance::GetName() const {
  return fName;
}

G4VPhysicalVolume &B01MassImportance::GetMassGeometry(){
  if (!fGeometry) {
    fGeometry = new B01SlobedConcreteShield;
    if (!fGeometry) {
      G4std::G4Exception("B01MassImportance::GetMassGeometry: failed to create B01SlobedConcreteShield");
    }
  }
  return fGeometry->GetWorldVolume();
}

const G4CellScorer *B01MassImportance::GetG4CellScorer(){
  return 0;
}

void B01MassImportance::PrepareSampling(){
  GetMassGeometry();

  fIStore = new G4IStore(fGeometry->GetWorldVolume());
  if (!fIStore) {
    G4std::G4Exception(" B01MassImportance::PrepareSampling: new failed to create G4IStore");
  }
  G4GeometryCell gWorldCell(fGeometry->GetWorldVolume(), -1);
  fIStore->AddImportanceGeometryCell(1, gWorldCell);

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
    }
  }
  
  
}

void B01MassImportance::ConfigureSampling(){

  fSampler = new G4MassGeometrySampler("neutron"); 
  if (!fSampler) {
    G4std::G4Exception("B01MassImportance::ConfigureSampling: new failed to create G4MassGeometrySampler!");
  }

  fSampler->PrepareImportanceSampling(fIStore); 
  if (fWeightRoulette) {
    fSampler->PrepareWeightRoulett();
  }
  fSampler->Configure();
  
}

void B01MassImportance::SetWeightRoulette(G4bool wroulette){
  fWeightRoulette = wroulette;
}

void B01MassImportance::PostRun(G4std::ostream *out) {
  *out << "B01MassImportance::PostRun: no results since no scoring applied" << G4endl; 
}
