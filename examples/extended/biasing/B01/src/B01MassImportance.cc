#include "B01MassImportance.hh"
#include "G4RunManager.hh"
#include "G4MassImportanceSampler.hh"
#include "G4IStore.hh"
#include "B01VGeometry.hh"
#include "B01SlobedConcreteShield.hh"
#include "B01Run.hh"
#include "G4Pstring.hh"

B01MassImportance::B01MassImportance()
  :
  fName("MassImportance"),
  fGeometry(0),
  fRun(0),
  fIStore(0),
  fSampler(0)
{}

B01MassImportance::~B01MassImportance(){
  if (fSampler) delete fSampler;
  if (fIStore) delete fIStore;
  if (fGeometry) delete fGeometry;
  if (fRun) delete fRun;
}

G4String B01MassImportance::GetName() const {
  return fName;
}
void B01MassImportance::Construct() {
  
  fGeometry = new B01SlobedConcreteShield;
  fRun = new B01Run;
  fRun->SetDetector(fGeometry->GetWorldVolume());
  fRun->Initialize();

  // create an importance store and fill it with the importance
  // per cell values

  fIStore = new G4IStore(fGeometry->GetWorldVolume());

  G4int i;
  for (i=1; i <= 18; i++) {
    G4String volname = fGeometry->GetCellName(i);
    G4double imp = pow(2,i-1);
    
    const G4VPhysicalVolume *pvol = fGeometry->
      GetPhysicalVolumeByName(volname);
    if (pvol) {
      G4GeometryCell gCell(*pvol, 0);
      fIStore->AddImportanceGeometryCell(imp, gCell);
    }
  }
  
  // create a scorer

  fSampler = new G4MassImportanceSampler(*fIStore,
					 "neutron"); 
  
  // to be done after fRun->Initialize()
  fSampler->Initialize(); 


}

void B01MassImportance::Run(G4int nevents) {
  if (!fRun) {
    G4cout << "B01MassImportance::Run: no BooRun constructed yet!" << G4endl;
  }
  else {
    fRun->BeamOn(nevents);
  }
}

void B01MassImportance::PostRun(G4std::ostream *out) {
  *out << "B01MassImportance::PostRun: no results since no scoring applied" << G4endl; 
}
