#include "G4DiscreteScatteringProcess.hh"
#include "G4DiscreteScatteringModel.hh"
#include "G4Electron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DiscreteScatteringProcess::G4DiscreteScatteringProcess(G4int iNumAngles)
: G4VEmProcess("DiscreteScattering", fElectromagnetic), 
  fIsInitialised(false), fNumAngles(iNumAngles)
{
  SetBuildTableFlag(true);
  SetStartFromNullFlag(false);
  SetMinKinEnergy(2*CLHEP::keV);
  SetMaxKinEnergy(100*CLHEP::MeV);
  SetIntegral(false);
  SetSplineFlag(false);
  SetProcessSubType(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DiscreteScatteringProcess::~G4DiscreteScatteringProcess(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool G4DiscreteScatteringProcess::IsApplicable(const G4ParticleDefinition& p){
  return (&p == G4Electron::Definition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DiscreteScatteringProcess::PrintInfo(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4DiscreteScatteringProcess::InitialiseProcess(const G4ParticleDefinition*)
{
  if(fIsInitialised) { return; }
  fIsInitialised = true;
  G4VEmModel* model = new G4DiscreteScatteringModel(fNumAngles);
  AddEmModel(1, model);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
