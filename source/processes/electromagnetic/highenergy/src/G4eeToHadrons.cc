//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4eeToHadrons.cc,v 1.2 2004/12/01 19:39:15 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4eeToHadrons
//
// Author:        Vladimir Ivanchenko on base of Michel Maire code
//
// Creation date: 02.08.2004
//
// Modifications:
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
//

//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eeToHadrons.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Gamma.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4eeToTwoPiModel.hh"
#include "G4eeToHadronsModel.hh"
#include "G4eeCrossSections.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4eeToHadrons::G4eeToHadrons(const G4String& name)
  : G4VEmProcess(name),
    lambdaFactor(0.1),
    csFactor(1.0),
    nModels(0),
    isInitialised(false)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eeToHadrons::~G4eeToHadrons()
{
  G4int n = models.size();
  if(n>1) {
    for(G4int i=1; i<n; i++) {
      delete models[i];
    }
  }
  delete cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToHadrons::InitialiseProcess(const G4ParticleDefinition*)
{
  if(!isInitialised) {
    isInitialised = true;

    G4int ver = 0;

    SetBuildTableFlag(false);
    SetIntegral(false);

    SetSecondaryParticle(G4Gamma::Gamma());
    SetParticle(G4Positron::Positron());

    currentCouple = 0;
    preStepMFP = DBL_MAX;
    thKineticEnergy = DBL_MAX;
    maxKineticEnergy = MaxKinEnergy();
    //maxKineticEnergy = 10.*TeV;
    const G4DataVector v;

    cross = new G4eeCrossSections();

    G4eeToHadronsModel* model =
      new G4eeToHadronsModel(new G4eeToTwoPiModel(cross), ver);
    models.push_back(model);
    model->SetHighEnergyLimit(maxKineticEnergy);
    model->Initialise(0, v);
    G4double emin = model->LowEnergyLimit();
    if(emin < thKineticEnergy) thKineticEnergy = emin;
    ekinMin.push_back(emin);
    ekinMax.push_back(model->HighEnergyLimit());
    ekinPeak.push_back(model->PeakEnergy());
    cumSum.push_back(0.0);
    AddEmModel(1, model);
    nModels = 1;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToHadrons::PrintInfoDefinition()
{
  G4VEmProcess::PrintInfoDefinition();

  G4cout << "      e+ annihilation into hadrons active above " 
         << thKineticEnergy/GeV << " GeV"
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4eeToHadrons::LambdaPhysicsVector(const G4MaterialCutsCouple*)
{
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeToHadrons::ComputeMeanFreePath(G4double kineticEnergy,
                                      const G4MaterialCutsCouple* couple)
{
  if( (couple != currentCouple) ||
      (abovePeak && kineticEnergy < mfpKineticEnergy) ) {

    ResetNumberOfInteractionLengthLeft();
    currentCouple = couple;
    mfpKineticEnergy = 0.0;
    for(G4int i=0; i<nModels; i++) {
      G4double peak = ekinPeak[i]; 
      if(peak <= kineticEnergy && peak > mfpKineticEnergy) 
        mfpKineticEnergy = peak;
    }  
    
    G4double cross = CrossSection(kineticEnergy, couple);
    if(mfpKineticEnergy == 0.0) {
       mfpKineticEnergy = kineticEnergy;
       abovePeak = false; 
    } else {
      G4double e = kineticEnergy*lambdaFactor;
      if(e > mfpKineticEnergy) mfpKineticEnergy = e;
      abovePeak = true; 
      G4double cross1 = CrossSection(mfpKineticEnergy, couple);
      if(cross1 > cross) cross = cross1;
    }
    preStepCS  = cross;
    preStepMFP = DBL_MAX;
    if(cross > 0.0) preStepMFP = 1.0/cross;
  }
  //  G4cout << "MFP: e= " << kineticEnergy << " mfpE= " << mfpKineticEnergy
  //       << " abovePeak= " << abovePeak << G4endl;
  return preStepMFP;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4DynamicParticle*>* G4eeToHadrons::GenerateSecondaries(
      const G4DynamicParticle* dp)
{
  std::vector<G4DynamicParticle*>* newp = 0;
  G4double e = dp->GetKineticEnergy();
  G4double cross = CrossSection(e, currentCouple);
  G4double q = preStepCS*G4UniformRand();
  if(q > cross) {
    ResetNumberOfInteractionLengthLeft();
  } else {
    for(G4int i=0; i<nModels; i++) {
      if(q < cumSum[i]) {
        newp = (models[i])->SampleSecondaries(currentCouple,dp,0.0,0.0);
        fParticleChange.ProposeTrackStatus(fStopAndKill);
        break;
      }
    }
  }
  return newp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToHadrons::SetCrossSecFactor(G4double fac)
{
  if(fac > 1.0) {
    csFactor = fac;
    G4cout << "### G4eeToHadrons: The cross section for G4eeToHadrons is  "
           << "increased by the Factor= " << csFactor << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
