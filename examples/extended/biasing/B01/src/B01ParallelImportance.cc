#include "B01ParallelImportance.hh"
#include "G4RunManager.hh"
#include "B01VGeometry.hh"
#include "B01ConcreteShield.hh"
#include "B01ParallelGeometry.hh"
#include "G4IStore.hh"
#include "G4GeometryCell.hh"
#include "G4ParallelImportanceSampler.hh"
#include "B01Run.hh"
#include "G4Pstring.hh"

B01ParallelImportance::B01ParallelImportance()
  :
  fName("ParallelImportance"),
  fMassGeometry(0),
  fRun(0),
  fParallelGeometry(0),
  fIStore(0),
  fSampler(0)
{}

B01ParallelImportance::~B01ParallelImportance(){
  if (fSampler) delete fSampler;
  if (fIStore) delete fIStore;
  if (fParallelGeometry) delete fParallelGeometry;
  if (fMassGeometry) delete fMassGeometry;
  if (fRun) delete fRun;
}

G4String B01ParallelImportance::GetName() const {
  return fName;
}

void B01ParallelImportance::Construct() {
  
  fMassGeometry = new B01ConcreteShield;
  fRun = new B01Run;
  fRun->SetDetector(fMassGeometry->GetWorldVolume());
  fRun->Initialize();
  
  fParallelGeometry = new B01ParallelGeometry;

  
  // create an importance store and fill it with the importance
  // per cell values
  fIStore = new G4IStore(fParallelGeometry->GetWorldVolume());

  G4int i;
  for (i=1; i <= 18; i++) {
    G4String volname = fParallelGeometry->GetCellName(i);
    G4double imp = pow(2,i-1);


    const G4VPhysicalVolume *pvol = fParallelGeometry->
      GetPhysicalVolumeByName(volname);
    if (pvol) {
      G4GeometryCell gCell(*pvol, 0);
      fIStore->AddImportanceGeometryCell(imp, gCell);
    }
  }
  
  // create importance sampler to importance sample neutrons
  // in acording to the paralle geometry
  fSampler = new
    G4ParallelImportanceSampler(fParallelGeometry->GetWorldVolume(),
				*fIStore, 
				"neutron");
  fSampler->Initialize();
 
}

void B01ParallelImportance::Run(G4int nevents) {
  if (!fRun) {
    G4cout << "B01ParallelImportance::Run: no BooRun constructed yet!" << G4endl;
  }
  else {
    fRun->BeamOn(nevents);
  }
}

void B01ParallelImportance::PostRun(G4std::ostream *out) {
  G4cout << "B01ParallelImportance::PostRun  no result since no scoring applied" << G4endl;
}
