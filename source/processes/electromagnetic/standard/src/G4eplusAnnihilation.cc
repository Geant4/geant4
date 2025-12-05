//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4eplusAnnihilation
//
// Author:        Vladimir Ivanchenko on base of Michel Maire code
//
// Creation date: 02.08.2004
//
// Modified by Michel Maire, Vladimir Ivanchenko and Daren Sawkey
//
// Introduced Quantum Entanglement  April 2021  John Allison
// This is activated by /process/em/QuantumEntanglement
// For e+e- -> gamma gamma, the gammas are "tagged" here
// and must be "analysed" in a Compton scattering process - see, for
// example, G4LivermorePolarizedComptonModel. Otherwise entanglement
// has no effect even if activated.

//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eplusAnnihilation.hh"
#include "G4PhysicalConstants.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4eeToTwoGammaModel.hh"
#include "G4EmBiasingManager.hh"
#include "G4EntanglementAuxInfo.hh"
#include "G4eplusAnnihilationEntanglementClipBoard.hh"
#include "G4SimplePositronAtRestModel.hh"
#include "G4AllisonPositronAtRestModel.hh"
#include "G4OrePowellAtRestModel.hh"
#include "G4PolarizedOrePowellAtRestModel.hh"
#include "G4EmParameters.hh"
#include "G4PhysicsModelCatalog.hh"
#include "G4Log.hh"

namespace
{
  constexpr G4double fTauPara = 0.12*CLHEP::ns;
  constexpr G4double fTauOrto = 142.*CLHEP::ns;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eplusAnnihilation::G4eplusAnnihilation(const G4String& name)
  : G4VEmProcess(name)
{
  SetCrossSectionType(fEmDecreasing);
  SetBuildTableFlag(false);
  SetStartFromNullFlag(false);
  SetSecondaryParticle(G4Gamma::Gamma());
  SetProcessSubType(fAnnihilation);
  enableAtRestDoIt = true;
  mainSecondaries = 2;
  fEntanglementModelID = G4PhysicsModelCatalog::GetModelID("model_GammaGammaEntanglement");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eplusAnnihilation::~G4eplusAnnihilation()
{
  delete f2GammaAtRestModel;
  delete f3GammaAtRestModel;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4eplusAnnihilation::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eplusAnnihilation::AtRestGetPhysicalInteractionLength(
                              const G4Track&, G4ForceCondition* condition)
{
  *condition = NotForced;
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eplusAnnihilation::InitialiseProcess(const G4ParticleDefinition*)
{
  if (!isInitialised) {
    isInitialised = true;
    if(nullptr == EmModel(0)) { SetEmModel(new G4eeToTwoGammaModel()); }
    EmModel(0)->SetLowEnergyLimit(MinKinEnergy());
    EmModel(0)->SetHighEnergyLimit(MaxKinEnergy());
    AddEmModel(1, EmModel(0));
  }
  auto param = G4EmParameters::Instance();

  // AtRest models should be chosen only once
  if (nullptr == f2GammaAtRestModel) {
    auto type = param->PositronAtRestModelType();
    if (type == fAllisonPositronium) {
      f2GammaAtRestModel = new G4AllisonPositronAtRestModel();
    } else if (type == fOrePowell) {
      f2GammaAtRestModel = new G4AllisonPositronAtRestModel();
      f3GammaAtRestModel = new G4OrePowellAtRestModel();
    } else if (type == fOrePowellPolar) {
      f2GammaAtRestModel = new G4AllisonPositronAtRestModel();
      f3GammaAtRestModel = new G4PolarizedOrePowellAtRestModel();
    } else {
      f2GammaAtRestModel = new G4SimplePositronAtRestModel();
    }
  }
  // Check that entanglement is switched on
  // It may be set by the UI command "/process/em/QuantumEntanglement true".
  fEntangled = param->QuantumEntanglement();
  fApplyCuts = param->ApplyCuts();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eplusAnnihilation::StreamProcessInfo(std::ostream&) const
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4eplusAnnihilation::AtRestDoIt(const G4Track& track,
						   const G4Step& step)
{
  // positron at rest should be killed
  fParticleChange.InitializeForPostStep(track);
  fParticleChange.SetProposedKineticEnergy(0.);
  fParticleChange.ProposeTrackStatus(fStopAndKill);

  auto couple = step.GetPreStepPoint()->GetMaterialCutsCouple();
  DefineMaterial(couple);
  
  G4double gammaCut = GetGammaEnergyCut();

  // apply cuts 
  if (fApplyCuts && gammaCut > CLHEP::electron_mass_c2) {
    fParticleChange.ProposeLocalEnergyDeposit(2*CLHEP::electron_mass_c2); 
    return &fParticleChange;
  }

  // sample secondaries
  secParticles.clear();
  G4double edep = 0.0;
  G4double ltime;
  if (nullptr != f3GammaAtRestModel &&
      G4UniformRand() < currentMaterial->GetIonisation()->GetOrtoPositroniumFraction()) {
    f3GammaAtRestModel->SampleSecondaries(secParticles, edep, couple->GetMaterial());
    ltime = fTauOrto;
  } else {
    f2GammaAtRestModel->SampleSecondaries(secParticles, edep, couple->GetMaterial());
    ltime = fTauPara;
  }

  // define new weight for primary and secondaries
  G4double weight = fParticleChange.GetParentWeight();				  
  std::size_t num0 = secParticles.size();
  
  // Russian roulette
  if (nullptr != biasManager) {
    G4int idx = couple->GetIndex();
    if (biasManager->SecondaryBiasingRegion(idx) &&
	!biasManager->GetDirectionalSplitting()) {
      G4VEmModel* mod = EmModel(0);
      G4double eloss = 0.0;
      weight *= biasManager->ApplySecondaryBiasing(secParticles, track, mod,
						   &fParticleChange, eloss, 
                                                   idx, gammaCut);
      edep += eloss;
    }
  }

  // save secondaries
  std::size_t num = secParticles.size();

  // Prepare a shared pointer only for two first gamma. If it is used, the
  // shared pointer is copied into the tracks through G4EntanglementAuxInfo.
  // This ensures the clip board lasts until both tracks are destroyed.
  // It is assumed that 2 first secondary particles are the most energetic gamma
  std::shared_ptr<G4eplusAnnihilationEntanglementClipBoard> clipBoard;
  if (fEntangled && num >= 2) {
    clipBoard = std::make_shared<G4eplusAnnihilationEntanglementClipBoard>();
    clipBoard->SetParentParticleDefinition(track.GetDefinition());
  }

  if (num > 0) {
    const G4double time = track.GetGlobalTime() - ltime*G4Log(G4UniformRand());
    const G4ThreeVector& pos = track.GetPosition();
    auto touch = track.GetTouchableHandle();
    for (std::size_t i=0; i<num; ++i) {
      G4DynamicParticle* dp = secParticles[i];
      G4Track* t = new G4Track(dp, time, pos);
      t->SetTouchableHandle(touch);
      if (fEntangled && i < 2) {
	// entangledgammagamma is only true when there are only two gammas
	// (See code above where entangledgammagamma is calculated.)
	if (nullptr != clipBoard) { 
	  if (i == 0) { // First gamma
	    clipBoard->SetTrackA(t);
	  } else if (i == 1) {  // Second gamma
	    clipBoard->SetTrackB(t);
	  }
	}
	t->SetAuxiliaryTrackInformation
	  (fEntanglementModelID, new G4EntanglementAuxInfo(clipBoard));
      }
      if (nullptr != biasManager) {
	t->SetWeight(weight * biasManager->GetWeight((G4int)i));
      } else {
	t->SetWeight(weight);
      }
      pParticleChange->AddSecondary(t);

      // define type of secondary
      if (i < num0) {
	t->SetCreatorModelID(secID);
      }
      else {
	t->SetCreatorModelID(biasID);
      }
    }
  }
  fParticleChange.ProposeLocalEnergyDeposit(edep);
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eplusAnnihilation::ProcessDescription(std::ostream& out) const
{
  out << "  Positron annihilation";
  G4VEmProcess::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
