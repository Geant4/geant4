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
// $Id: G4eeToHadronsMultiModel.cc,v 1.9 2010-10-26 14:15:40 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4eeToHadronsMultiModel
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 02.08.2004
//
// Modifications:
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
//

//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eeToHadronsMultiModel.hh"
#include "G4eeToTwoPiModel.hh"
#include "G4eeTo3PiModel.hh"
#include "G4eeToPGammaModel.hh"
#include "G4ee2KNeutralModel.hh"
#include "G4ee2KChargedModel.hh"
#include "G4eeCrossSections.hh"
#include "G4Vee2hadrons.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4eeToHadronsMultiModel::G4eeToHadronsMultiModel(G4int ver, const G4String& name)
  : G4VEmModel(name),
    csFactor(1.0),
    nModels(0),
    verbose(ver),
    isInitialised(false)
{
  thKineticEnergy  = DBL_MAX;
  maxKineticEnergy = 1.2*GeV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eeToHadronsMultiModel::~G4eeToHadronsMultiModel()
{
  G4int n = models.size();
  if(n>0) {
    for(G4int i=0; i<n; i++) {
      delete models[i];
    }
  }
  delete cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToHadronsMultiModel::Initialise(const G4ParticleDefinition*, 
					 const G4DataVector&)
{
  if(!isInitialised) {
    isInitialised = true;

    cross = new G4eeCrossSections();

    G4eeToTwoPiModel* m2pi = new G4eeToTwoPiModel(cross);
    m2pi->SetHighEnergy(maxKineticEnergy);
    AddEEModel(m2pi);

    G4eeTo3PiModel* m3pi1 = new G4eeTo3PiModel(cross);
    m3pi1->SetHighEnergy(0.95*GeV);
    AddEEModel(m3pi1);

    G4eeTo3PiModel* m3pi2 = new G4eeTo3PiModel(cross);
    m3pi2->SetLowEnergy(0.95*GeV);
    m3pi2->SetHighEnergy(maxKineticEnergy);
    AddEEModel(m3pi2);

    G4ee2KChargedModel* m2kc = new G4ee2KChargedModel(cross);
    m2kc->SetHighEnergy(maxKineticEnergy);
    AddEEModel(m2kc);

    G4ee2KNeutralModel* m2kn = new G4ee2KNeutralModel(cross);
    m2kn->SetHighEnergy(maxKineticEnergy);
    AddEEModel(m2kn);

    G4eeToPGammaModel* mpg1 = new G4eeToPGammaModel(cross,"pi0");
    mpg1->SetLowEnergy(0.7*GeV);
    mpg1->SetHighEnergy(maxKineticEnergy);
    AddEEModel(mpg1);

    G4eeToPGammaModel* mpg2 = new G4eeToPGammaModel(cross,"eta");
    mpg2->SetLowEnergy(0.7*GeV);
    mpg2->SetHighEnergy(maxKineticEnergy);
    AddEEModel(mpg2);

    nModels = models.size();

    fParticleChange = GetParticleChangeForGamma();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToHadronsMultiModel::AddEEModel(G4Vee2hadrons* mod)
{
  G4eeToHadronsModel* model = new G4eeToHadronsModel(mod, verbose);
  model->SetLowEnergyLimit(LowEnergyLimit());
  model->SetHighEnergyLimit(HighEnergyLimit());
  models.push_back(model);
  G4double elow = mod->ThresholdEnergy();
  ekinMin.push_back(elow);
  if(thKineticEnergy > elow) thKineticEnergy = elow;
  ekinMax.push_back(mod->HighEnergy());
  ekinPeak.push_back(mod->PeakEnergy());
  cumSum.push_back(0.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeToHadronsMultiModel::CrossSectionPerVolume(
				      const G4Material* mat,
				      const G4ParticleDefinition* p,
				      G4double kineticEnergy,
				      G4double, G4double)
{
  return mat->GetElectronDensity()*
    ComputeCrossSectionPerElectron(p, kineticEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeToHadronsMultiModel::ComputeCrossSectionPerAtom(
                                      const G4ParticleDefinition* p,
				      G4double kineticEnergy,
				      G4double Z, G4double,
				      G4double, G4double)
{
  return Z*ComputeCrossSectionPerElectron(p, kineticEnergy);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToHadronsMultiModel::SampleSecondaries(std::vector<G4DynamicParticle*>* newp,
						const G4MaterialCutsCouple* couple,
						const G4DynamicParticle* dp,
						G4double, G4double)
{
  G4double kinEnergy = dp->GetKineticEnergy();
  if (kinEnergy > thKineticEnergy) {
    G4double q = cumSum[nModels-1]*G4UniformRand();
    for(G4int i=0; i<nModels; i++) {
      if(q <= cumSum[i]) {
        (models[i])->SampleSecondaries(newp, couple,dp);
	if(newp->size() > 0) fParticleChange->ProposeTrackStatus(fStopAndKill);
	break;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToHadronsMultiModel::PrintInfo()
{
  if(verbose > 0) {
    G4double e1 = 0.5*thKineticEnergy*thKineticEnergy/electron_mass_c2 
      - 2.0*electron_mass_c2; 
    G4double e2 = 0.5*maxKineticEnergy*maxKineticEnergy/electron_mass_c2 
      - 2.0*electron_mass_c2; 
    G4cout << "      e+ annihilation into hadrons active from "
           << e1/GeV << " GeV to " << e2/GeV << " GeV"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToHadronsMultiModel::SetCrossSecFactor(G4double fac)
{
  if(fac > 1.0) {
    csFactor = fac;
    if(verbose > 0)
      G4cout << "### G4eeToHadronsMultiModel: The cross section for G4eeToHadronsMultiModel "
             << " is increased by the Factor= " << csFactor << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
