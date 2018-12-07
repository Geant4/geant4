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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4eplusAnnihilation::G4eplusAnnihilation(const G4String& name)
  : G4VEmProcess(name), isInitialised(false)
{
  theGamma = G4Gamma::Gamma();
  SetIntegral(true);
  SetBuildTableFlag(false);
  SetStartFromNullFlag(false);
  SetSecondaryParticle(theGamma);
  SetProcessSubType(fAnnihilation);
  enableAtRestDoIt = true;
  mainSecondaries = 2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eplusAnnihilation::~G4eplusAnnihilation()
{}

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
    if(!EmModel(0)) { SetEmModel(new G4eeToTwoGammaModel()); }
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
  size_t idx = CurrentMaterialCutsCoupleIndex();
  G4double ene(0.0);
  G4VEmModel* model = SelectModel(ene, idx);

  // define new weight for primary and secondaries
  G4double weight = fParticleChange.GetParentWeight();

  // sample secondaries
  secParticles.clear();
  G4double gammaCut = GetGammaEnergyCut();
  model->SampleSecondaries(&secParticles, MaterialCutsCouple(), 
			   track.GetDynamicParticle(), gammaCut);
 
  G4int num0 = secParticles.size();

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
  G4int num = secParticles.size();
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
          if (biasManager) {
            t->SetWeight(biasManager->GetWeight(i));
          } else {
            t->SetWeight(weight);
          }
          pParticleChange->AddSecondary(t);

          // define type of secondary
          if(i < mainSecondaries) { t->SetCreatorModelIndex(secID); }
          else if(i < num0) {
            if(p == theGamma) { 
              t->SetCreatorModelIndex(fluoID);
            } else {
              t->SetCreatorModelIndex(augerID);
	    }
	  } else {
            t->SetCreatorModelIndex(biasID);
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
