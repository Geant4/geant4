#include "B01ParallelImportance.hh"
#include "G4RunManager.hh"
#include "B01VGeometry.hh"
#include "B01ConcreteShield.hh"
#include "B01ParallelGeometry.hh"
#include "G4IStore.hh"
#include "G4GeometryCell.hh"
#include "G4ParallelGeometrySampler.hh"


B01ParallelImportance::B01ParallelImportance()
  :
  fName("ParallelImportance"),
  fWeightRoulette(0),
  fMassGeometry(0),
  fParallelGeometry(0),
  fIStore(0),
  fSampler(0)
{}

B01ParallelImportance::~B01ParallelImportance(){
  if (fSampler) {
    delete fSampler;
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

const G4String &B01ParallelImportance::GetName() const {
  return fName;
}

G4VPhysicalVolume &B01ParallelImportance::GetMassGeometry(){
  if (!fMassGeometry) {
    fMassGeometry = new B01ConcreteShield;
    if (!fMassGeometry) {
      G4std::G4Exception("B01ParallelImportance::GetMassGeometry: new failed to create B01ConcreteShield!");
    }
  }
  return fMassGeometry->GetWorldVolume();
}

const G4CellScorer *B01ParallelImportance::GetG4CellScorer(){
  return 0;  
}

void B01ParallelImportance::PrepareSampling(){
  GetMassGeometry();
  fParallelGeometry = new B01ParallelGeometry;
  if (!fParallelGeometry) {
    G4std::G4Exception("B01ParallelImportance::PrepareSampling: new failed to create B01ParallelGeometry!");
  }
  
  // create an importance store and fill it with the importance
  // per cell values
  const G4VPhysicalVolume &pworld = fParallelGeometry->GetWorldVolume();
  fIStore = new G4IStore(pworld);
  if (!fIStore) {
    G4std::G4Exception("B01ParallelImportance::PrepareSampling: new failed to create  G4IStore!");
  }
  // adding GeometryCell for world volume. ReplicaNumer = -1 !
  G4GeometryCell gWorldCell(pworld, -1);
  fIStore->AddImportanceGeometryCell(1, gWorldCell);

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
    }
  }

}

void B01ParallelImportance::ConfigureSampling(){
  fSampler = new
    G4ParallelGeometrySampler(fParallelGeometry->GetWorldVolume(),
			      "neutron");
  if (!fSampler) {
    G4std::G4Exception("B01ParallelImportance::ConfigureSampling: new failed to create G4ParallelGeometrySampler!");
  }
  fSampler->PrepareImportanceSampling(fIStore); 
  if (fWeightRoulette) {
    fSampler->PrepareWeightRoulett();
  }
  fSampler->Configure();

}

void B01ParallelImportance::SetWeightRoulette(G4bool wroulette){
  fWeightRoulette = wroulette;
}

void B01ParallelImportance::PostRun(G4std::ostream *out) {
  G4std::G4cout << "B01ParallelImportance::PostRun  no result since no scoring applied" << G4endl;
}


