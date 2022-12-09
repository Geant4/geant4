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
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4eeToTwoGammaModel.hh"
#include "G4EmBiasingManager.hh"
#include "G4EntanglementAuxInfo.hh"
#include "G4eplusAnnihilationEntanglementClipBoard.hh"
#include "G4EmParameters.hh"
#include "G4PhysicsModelCatalog.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eplusAnnihilation::G4eplusAnnihilation(const G4String& name)
  : G4VEmProcess(name)
{
  theGamma = G4Gamma::Gamma();
  theElectron = G4Electron::Electron();
  SetCrossSectionType(fEmDecreasing);
  SetBuildTableFlag(false);
  SetStartFromNullFlag(false);
  SetSecondaryParticle(theGamma);
  SetProcessSubType(fAnnihilation);
  enableAtRestDoIt = true;
  mainSecondaries = 2;
  fEntanglementModelID = G4PhysicsModelCatalog::GetModelID("model_GammaGammaEntanglement");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eplusAnnihilation::~G4eplusAnnihilation() = default;

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
  if(!isInitialised) {
    isInitialised = true;
    if(nullptr == EmModel(0)) { SetEmModel(new G4eeToTwoGammaModel()); }
    EmModel(0)->SetLowEnergyLimit(MinKinEnergy());
    EmModel(0)->SetHighEnergyLimit(MaxKinEnergy());
    AddEmModel(1, EmModel(0));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eplusAnnihilation::StreamProcessInfo(std::ostream&) const
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4eplusAnnihilation::AtRestDoIt(const G4Track& track,
						   const G4Step& step)
// Performs the e+ e- annihilation when both particles are assumed at rest.
{
  fParticleChange.InitializeForPostStep(track);

  DefineMaterial(track.GetMaterialCutsCouple());
  G4int idx = (G4int)CurrentMaterialCutsCoupleIndex();
  G4double ene(0.0);
  G4VEmModel* model = SelectModel(ene, idx);

  // define new weight for primary and secondaries
  G4double weight = fParticleChange.GetParentWeight();

  // sample secondaries
  secParticles.clear();
  G4double gammaCut = GetGammaEnergyCut();
  model->SampleSecondaries(&secParticles, MaterialCutsCouple(), 
			   track.GetDynamicParticle(), gammaCut);

  G4int num0 = (G4int)secParticles.size();

  // splitting or Russian roulette
  if(biasManager) {
    if(biasManager->SecondaryBiasingRegion(idx)) {
      G4double eloss = 0.0;
      weight *= biasManager->ApplySecondaryBiasing(
	secParticles, track, model, &fParticleChange, eloss, 
        idx, gammaCut, step.GetPostStepPoint()->GetSafety());
      if(eloss > 0.0) {
        eloss += fParticleChange.GetLocalEnergyDeposit();
        fParticleChange.ProposeLocalEnergyDeposit(eloss);
      }
    }
  }

  // save secondaries
  G4int num = (G4int)secParticles.size();

  // Check that entanglement is switched on... (the following flag is
  // set by /process/em/QuantumEntanglement).
  G4bool entangled = G4EmParameters::Instance()->QuantumEntanglement();
  // ...and that we have two gammas with both gammas' energies above
  // gammaCut (entanglement is only programmed for e+ e- -> gamma gamma).
  G4bool entangledgammagamma = false;
  if (entangled) {
    if (num == 2) {
      entangledgammagamma = true;
      for (const auto* p: secParticles) {
        if (p->GetDefinition() != theGamma ||
            p->GetKineticEnergy() < gammaCut) {
          entangledgammagamma = false;
        }
      }
    }
  }

  // Prepare a shared pointer for psossible use below. If it is used, the
  // shared pointer is copied into the tracks through G4EntanglementAuxInfo.
  // This ensures the clip board lasts until both tracks are destroyed.
  std::shared_ptr<G4eplusAnnihilationEntanglementClipBoard> clipBoard;
  if (entangledgammagamma) {
    clipBoard = std::make_shared<G4eplusAnnihilationEntanglementClipBoard>();
    clipBoard->SetParentParticleDefinition(track.GetDefinition());
  }

  if(num > 0) {

    fParticleChange.SetNumberOfSecondaries(num);
    G4double edep = fParticleChange.GetLocalEnergyDeposit();
    G4double time = track.GetGlobalTime();
     
    for (G4int i=0; i<num; ++i) {
      if (secParticles[i]) {
        G4DynamicParticle* dp = secParticles[i];
        const G4ParticleDefinition* p = dp->GetParticleDefinition();
        G4double e = dp->GetKineticEnergy();
        G4bool good = true;
        if(ApplyCuts()) {
          if (p == theGamma) {
            if (e < gammaCut) { good = false; }
          } else if (p == theElectron) {
            if (e < GetElectronEnergyCut()) { good = false; }
          }
          // added secondary if it is good
        }
        if (good) { 
          G4Track* t = new G4Track(dp, time, track.GetPosition());
          t->SetTouchableHandle(track.GetTouchableHandle());
          if (entangledgammagamma) {
	    // entangledgammagamma is only true when there are only two gammas
	    // (See code above where entangledgammagamma is calculated.)
            if (i == 0) { // First gamma
              clipBoard->SetTrackA(t);
            } else if (i == 1) {  // Second gamma
              clipBoard->SetTrackB(t);
            }
            t->SetAuxiliaryTrackInformation
            (fEntanglementModelID,new G4EntanglementAuxInfo(clipBoard));
          }
          if (biasManager) {
            t->SetWeight(weight * biasManager->GetWeight(i));
          } else {
            t->SetWeight(weight);
          }
          pParticleChange->AddSecondary(t);

          // define type of secondary
          if(i < mainSecondaries) { t->SetCreatorModelID(secID); }
          else if(i < num0) {
            if(p == theGamma) { 
              t->SetCreatorModelID(fluoID);
            } else {
              t->SetCreatorModelID(augerID);
	    }
	  } else {
            t->SetCreatorModelID(biasID);
          }
          /* 
          G4cout << "Secondary(post step) has weight " << t->GetWeight() 
                 << ", Ekin= " << t->GetKineticEnergy()/MeV << " MeV "
                 << GetProcessName() << " fluoID= " << fluoID
                 << " augerID= " << augerID <<G4endl;
          */
        } else {
          delete dp;
          edep += e;
        }
      } 
    }
    fParticleChange.ProposeLocalEnergyDeposit(edep);
  }
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eplusAnnihilation::ProcessDescription(std::ostream& out) const
{
  out << "  Positron annihilation";
  G4VEmProcess::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
