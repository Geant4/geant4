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
// $Id: G4eeToHadronsModel.cc,v 1.1 2004-08-19 16:30:06 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4eeToHadronsModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 12.08.2003
//
// Modifications:
//
//
// -------------------------------------------------------------------
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eeToHadronsModel.hh"
#include "Randomize.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4PionPlus.hh"
#include "Randomize.hh"
#include "G4Vee2hadrons.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLogVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eeToHadronsModel::G4eeToHadronsModel(const G4Vee2hadrons* m, const G4String& nam)
  : G4VEmModel(nam),
  model(m),
  crossPerElectron(0),
  crossBornPerElectron(0),
  isInitialised(false),
  nbins(100)
{
  highKinEnergy = 10000.*TeV;
  lowKinEnergy  = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eeToHadronsModel::~G4eeToHadronsModel()
{
  delete model;
  delete crossPerElectron;
  delete crossBornPerElectron;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeToHadronsModel::HighEnergyLimit(const G4ParticleDefinition*)
{
  return highKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeToHadronsModel::LowEnergyLimit(const G4ParticleDefinition*)
{
  return lowKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToHadronsModel::SetHighEnergyLimit(G4double e)
{
  highKinEnergy = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToHadronsModel::SetLowEnergyLimit(G4double e)
{
  lowKinEnergy = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4eeToHadronsModel::IsInCharge(const G4ParticleDefinition* p)
{
  return (p == G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToHadronsModel::Initialise(const G4ParticleDefinition*,
                                    const G4DataVector&)
{
  if(isInitialised) return;
  isInitialised  = true;

  emin  = model->ThresholdEnergy();
  emax = 2.0*electron_mass_c2*sqrt(1.0 + 0.5*highKinEnergy/electron_mass_c2);
  if(emin > emax) emin = emax;

  lowKinEnergy  = 0.5*emin*emin/electron_mass_c2 - 2.0*electron_mass_c2;

  epeak = std::min(model->PeakEnergy(), emax);
  peakKinEnergy  = 0.5*epeak*epeak/electron_mass_c2 - 2.0*electron_mass_c2;

  if(lowKinEnergy < peakKinEnergy) {
    crossBornPerElectron = model->PhysicsVector(emin, emax); 
    crossPerElectron     = crossBornPerElectron;
    G4int nbins = crossPerElectron->GetVectorLength();
    for(G4int i=0; i<nbins; i++) {
      G4double e  = crossPerElectron->GetLowEdgeEnergy(i);
      G4double cs = model->ComputeCrossSection(e);     
      crossBornPerElectron->PutValue(i, cs);
    }
    ComputeCMCrossSectionPerElectron();         
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeToHadronsModel::ComputeDEDX(const G4MaterialCutsCouple*,
                                        const G4ParticleDefinition*,
                                              G4double,
                                              G4double)
{
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeToHadronsModel::CrossSection(const G4MaterialCutsCouple* couple,
                                          const G4ParticleDefinition*,
                                                G4double kineticEnergy,
                                                G4double,
                                                G4double)
{
  G4double cross = 0.0;
  if(crossPerElectron) {
    G4bool b;   
    G4double e  = 2.0*electron_mass_c2*sqrt(1.0 + 0.5*kineticEnergy/electron_mass_c2);
    cross = (couple->GetMaterial()->GetElectronDensity())*(crossPerElectron->GetValue(e, b));
  }
  //  G4cout << "e= " << kineticEnergy << " cross= " << cross << G4endl;
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DynamicParticle* G4eeToHadronsModel::SampleSecondary(
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle*,
                                   G4double,
                                   G4double)
{
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4DynamicParticle*>* G4eeToHadronsModel::SampleSecondaries(
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle* dParticle,
                                   G4double,
                                   G4double)
{
  std::vector<G4DynamicParticle*>* newp = 0;
  if(crossPerElectron) {
    G4double t = dParticle->GetKineticEnergy();
    G4double e = 2.0*electron_mass_c2*sqrt(1.0 + 0.5*t/electron_mass_c2);
    G4LorentzVector inlv = dParticle->Get4Momentum();
    G4ThreeVector inBoost = inlv.boostVector();
    if(e > emin) {
      G4DynamicParticle* gamma = GenerateCMPhoton(e);
      G4LorentzVector gLv = gamma->Get4Momentum();
      G4LorentzVector lv(0.0,0.0,0.0,e);
      lv -= gLv;
      G4double m = lv.m();
      G4ThreeVector boost = lv.boostVector();
      const G4ThreeVector dir = gamma->GetMomentumDirection();
      newp = model->SampleSecondaries(m, dir);
      if(newp) {
        G4int np = newp->size();
        for(G4int j=0; j<np; j++) {
          G4DynamicParticle* dp = (*newp)[j];
          G4LorentzVector v = dp->Get4Momentum();
          v.boost(boost);
          v.boost(inBoost);
          dp->Set4Momentum(v);
	}
      } else {
        newp = new std::vector<G4DynamicParticle*>;
      }
      gLv.boost(inBoost);
      gamma->Set4Momentum(gLv);
      newp->push_back(gamma);
    }
  }
  return newp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToHadronsModel::ComputeCMCrossSectionPerElectron()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DynamicParticle* G4eeToHadronsModel::GenerateCMPhoton(G4double)
{
  G4DynamicParticle* gamma = 0;
  return gamma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

