
//
// GEANT4 Class file
//
// File name:     G4AdjointhMultipleScattering
//
// Author:        Desorgher Laurent
//
// Creation date: 03.06.2009 cloned from G4hMultipleScattering by U.Laszlo with slight modification for adjoint_ion. 
//
// -----------------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4AdjointhMultipleScattering.hh"
#include "G4UrbanMscModel.hh"
#include "G4UrbanMscModel90.hh"
#include "G4MscStepLimitType.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4AdjointhMultipleScattering::G4AdjointhMultipleScattering(const G4String& processName)
  : G4VMultipleScattering(processName)
{
  isInitialized = false;  
  isIon         = false;
  SetStepLimitType(fMinimal);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4AdjointhMultipleScattering::~G4AdjointhMultipleScattering()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4AdjointhMultipleScattering::IsApplicable (const G4ParticleDefinition& p)
{
  return (p.GetPDGCharge() != 0.0 && !p.IsShortLived());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AdjointhMultipleScattering::InitialiseProcess(const G4ParticleDefinition* p)
{
  // Modification of parameters between runs
  if(isInitialized) {
    if (p->GetParticleType() != "adjoint_nucleus" && p->GetPDGMass() < GeV) {
      mscUrban->SetStepLimitType(StepLimitType());
      mscUrban->SetLateralDisplasmentFlag(LateralDisplasmentFlag());
      mscUrban->SetSkin(Skin());
      mscUrban->SetRangeFactor(RangeFactor());
      mscUrban->SetGeomFactor(GeomFactor());
    }
    return;
  }

  // defaults for ions, which cannot be overwritten
  if (p->GetParticleType() == "adjoint_nucleus" || p->GetPDGMass() > GeV) {
    SetStepLimitType(fMinimal);
    SetLateralDisplasmentFlag(false);
    SetBuildLambdaTable(false);
    if(p->GetParticleType() == "adjoint_nucleus") isIon = true;
  }

  // initialisation of parameters
  G4String part_name = p->GetParticleName();
  mscUrban = new G4UrbanMscModel90();

  mscUrban->SetStepLimitType(StepLimitType());
  mscUrban->SetLateralDisplasmentFlag(LateralDisplasmentFlag());
  mscUrban->SetSkin(Skin());
  mscUrban->SetRangeFactor(RangeFactor());
  mscUrban->SetGeomFactor(GeomFactor());

  AddEmModel(1,mscUrban);
  isInitialized = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AdjointhMultipleScattering::PrintInfo()
{
  G4cout << "      RangeFactor= " << RangeFactor()
	 << ", step limit type: " << StepLimitType()
         << ", lateralDisplacement: " << LateralDisplasmentFlag()
	 << ", skin= " << Skin()  
    //	 << ", geomFactor= " << GeomFactor()  
	 << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*G4double G4AdjointhMultipleScattering::AlongStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double,
                             G4double currentMinimalStep,
                             G4double& currentSafety,
                             G4GPILSelection* selection)
{
  // get Step limit proposed by the process
  valueGPILSelectionMSC = NotCandidateForSelection;

  G4double escaled = track.GetKineticEnergy();
  if(isIon) escaled *= track.GetDynamicParticle()->GetMass()/proton_mass_c2;

  G4double steplength = GetMscContinuousStepLimit(track,
						  escaled,
						  currentMinimalStep,
						  currentSafety);
  // G4cout << "StepLimit= " << steplength << G4endl;
  // set return value for G4GPILSelection
  *selection = valueGPILSelectionMSC;
  return  steplength;
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

