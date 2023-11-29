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
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4eeToTwoPiModel.hh"
#include "G4eeTo3PiModel.hh"
#include "G4eeToPGammaModel.hh"
#include "G4ee2KNeutralModel.hh"
#include "G4ee2KChargedModel.hh"
#include "G4eeCrossSections.hh"
#include "G4Vee2hadrons.hh"
#include "G4Positron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4eeToHadronsMultiModel::G4eeToHadronsMultiModel(G4int ver, 
  const G4String& mname) : G4VEmModel(mname), verbose(ver)
{
  maxKineticEnergy = 4.521*CLHEP::GeV;  //crresponding to 10TeV in lab
  delta = 1.0*CLHEP::MeV;  //for bin width
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eeToHadronsMultiModel::~G4eeToHadronsMultiModel()
{
  delete cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToHadronsMultiModel::Initialise(const G4ParticleDefinition*, 
					 const G4DataVector& cuts) 
{
  if(!isInitialised) {
    isInitialised = true;

    //G4cout<<"###Initialise in HadronMultiModel###"<<G4endl;

    cross = new G4eeCrossSections();

    G4eeToTwoPiModel* m2pi = 
      new G4eeToTwoPiModel(cross,maxKineticEnergy,delta);
    AddEEModel(m2pi,cuts);

    G4eeTo3PiModel* m3pi = 
      new G4eeTo3PiModel(cross,maxKineticEnergy,delta);
    AddEEModel(m3pi,cuts);

    G4ee2KChargedModel* m2kc = 
      new G4ee2KChargedModel(cross,maxKineticEnergy,delta);
    AddEEModel(m2kc,cuts);

    G4ee2KNeutralModel* m2kn = 
      new G4ee2KNeutralModel(cross,maxKineticEnergy,delta);
    AddEEModel(m2kn,cuts);

    G4eeToPGammaModel* mpg1 = 
      new G4eeToPGammaModel(cross,"pi0",maxKineticEnergy,delta);
    AddEEModel(mpg1,cuts);

    G4eeToPGammaModel* mpg2 = 
      new G4eeToPGammaModel(cross,"eta",maxKineticEnergy,delta);
    AddEEModel(mpg2,cuts);

    nModels = (G4int)models.size();

    fParticleChange = GetParticleChangeForGamma();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToHadronsMultiModel::AddEEModel(G4Vee2hadrons* mod,
					 const G4DataVector& cuts)
{
  G4eeToHadronsModel* model = new G4eeToHadronsModel(mod, verbose);
  models.push_back(model);
  G4double elow = mod->LowEnergy();   
  ekinMin.push_back(elow);
  if(thKineticEnergy > elow) { thKineticEnergy = elow; }
  ekinMax.push_back(mod->HighEnergy()); 
  ekinPeak.push_back(mod->PeakEnergy());
  cumSum.push_back(0.0);

  const G4ParticleDefinition* positron = G4Positron::Positron();
  model->Initialise(positron,cuts);
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

G4double G4eeToHadronsMultiModel::ComputeCrossSectionPerElectron(const G4ParticleDefinition*,
								 G4double kineticEnergy,
								 G4double, G4double)
{
  G4double res = 0.0;

  G4double energy = LabToCM(kineticEnergy);

  if (energy > thKineticEnergy) {
    for(G4int i=0; i<nModels; i++) {
      if(energy >= ekinMin[i] && energy <= ekinMax[i]){
        res += (models[i])->ComputeCrossSectionPerElectron(0,energy);
      }
      cumSum[i] = res;
    }
  }
  return res*csFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToHadronsMultiModel::SampleSecondaries(
           std::vector<G4DynamicParticle*>* newp,
	   const G4MaterialCutsCouple* couple,
	   const G4DynamicParticle* dp,
	   G4double, G4double)
{
  G4double kinEnergy = dp->GetKineticEnergy();
  G4double energy = LabToCM(kinEnergy);
  if (energy > thKineticEnergy) {
    G4double q = cumSum[nModels-1]*G4UniformRand();
    for(G4int i=0; i<nModels; ++i) {
      if(q <= cumSum[i]) {
        (models[i])->SampleSecondaries(newp, couple,dp);
	if(newp->size() > 0) {
	  fParticleChange->ProposeTrackStatus(fStopAndKill);
	}
	break;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToHadronsMultiModel::ModelDescription(std::ostream& outFile) const
{
  if(verbose > 0) {
    G4double e1 = 0.5*thKineticEnergy*thKineticEnergy/electron_mass_c2 
      - 2.0*electron_mass_c2; 
    G4double e2 = 0.5*maxKineticEnergy*maxKineticEnergy/electron_mass_c2 
      - 2.0*electron_mass_c2; 
    outFile << "      e+ annihilation into hadrons active from "
	    << e1/GeV << " GeV to " << e2/GeV << " GeV" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToHadronsMultiModel::SetCrossSecFactor(G4double fac)
{
  if(fac > 1.0) {
    csFactor = fac;
    if(verbose > 0) {
      G4cout << "### G4eeToHadronsMultiModel: The cross section for "
	     << "G4eeToHadronsMultiModel is increased by " 
	     << csFactor << " times" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
