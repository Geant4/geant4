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
// $Id: G4eplusAnnihilation.cc,v 1.22 2005/05/12 11:06:43 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $
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
// Modifications:
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivanchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivanchenko)
// 03-05-05 suppress Integral option (mma)
// 04-05-05, Make class to be default (V.Ivanchenko)
//

//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eplusAnnihilation.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Gamma.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4eeToTwoGammaModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4eplusAnnihilation::G4eplusAnnihilation(const G4String& name)
  : G4VEmProcess(name), isInitialised(false)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eplusAnnihilation::~G4eplusAnnihilation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eplusAnnihilation::InitialiseProcess(const G4ParticleDefinition*)
{
  if(!isInitialised) {
    isInitialised = true;
    //    SetVerboseLevel(3);
    SetBuildTableFlag(true);
    SetStartFromNullFlag(false);
    SetSecondaryParticle(G4Gamma::Gamma());
    G4double emin = 0.1*keV;
    G4double emax = 100.*TeV;
    SetLambdaBinning(120);
    SetMinKinEnergy(emin);
    SetMaxKinEnergy(emax);
    G4VEmModel* em = new G4eeToTwoGammaModel();
    em->SetLowEnergyLimit(emin);
    em->SetHighEnergyLimit(emax);
    AddEmModel(1, em);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eplusAnnihilation::PrintInfo()
{
  G4cout << "      Heilter model of formula of annihilation into 2 photons"
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4eplusAnnihilation::AtRestDoIt(const G4Track& aTrack,
                                                     const G4Step& )
//
// Performs the e+ e- annihilation when both particles are assumed at rest.
// It generates two back to back photons with energy = electron_mass.
// The angular distribution is isotropic.
// GEANT4 internal units
//
// Note : Effects due to binding of atomic electrons are negliged.
{
  fParticleChange.InitializeForPostStep(aTrack);

  // Below gamma production threshold
  if (GetGammaEnergyCut() > electron_mass_c2) {
    fParticleChange.ProposeLocalEnergyDeposit(2.0*electron_mass_c2);
    
  } else {   // Real gamma production 
    fParticleChange.SetNumberOfSecondaries(2);

    G4double cosTeta = 2.*G4UniformRand()-1. , sinTeta = sqrt(1.-cosTeta*cosTeta);
    G4double phi     = twopi * G4UniformRand();
    G4ThreeVector direction (sinTeta*cos(phi), sinTeta*sin(phi), cosTeta);
    fParticleChange.AddSecondary( new G4DynamicParticle (G4Gamma::Gamma(),
                                            direction, electron_mass_c2) );
    fParticleChange.AddSecondary( new G4DynamicParticle (G4Gamma::Gamma(),
                                           -direction, electron_mass_c2) );
  }
  // Kill the incident positron
  //
  fParticleChange.ProposeTrackStatus(fStopAndKill);
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
