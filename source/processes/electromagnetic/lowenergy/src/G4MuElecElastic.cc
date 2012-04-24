//
// G4MuElecElastic.cc, 2011/08/29 A.Valentin, M. Raine
//
// Based on the following publications
//	    - Geant4 physics processes for microdosimetry simulation:
//	    very low energy electromagnetic models for electrons in Si,
//	    to be published in TNS
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "G4MuElecElastic.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4MuElecElastic::G4MuElecElastic(const G4String& processName,
  G4ProcessType type):G4VEmProcess (processName, type),
    isInitialised(false)
{
  SetProcessSubType(51);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4MuElecElastic::~G4MuElecElastic()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4MuElecElastic::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4Electron::Electron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuElecElastic::InitialiseProcess(const G4ParticleDefinition*)
{
  if(!isInitialised) 
  {
    isInitialised = true;
    SetBuildTableFlag(false);
    if(!Model()) SetModel(new G4MuElecElasticModel);
    Model()->SetLowEnergyLimit(5*eV);
    Model()->SetHighEnergyLimit(50*keV);
    AddEmModel(1, Model());
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuElecElastic::PrintInfo()
{
  G4cout
    << " Total cross sections computed from " << Model()->GetName() << " model"
    << G4endl;
}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
