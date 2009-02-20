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
// $Id: G4eplusAnnihilation.cc,v 1.30 2009-02-20 12:06:37 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
// 04-05-05 Make class to be default (V.Ivanchenko)
// 25-01-06 remove cut dependance in AtRestDoIt (mma)
// 09-08-06 add SetModel(G4VEmModel*) (mma)
// 12-09-06, move SetModel(G4VEmModel*) in G4VEmProcess (mma)
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
{
  SetProcessSubType(fAnnihilation);
  enableAtRestDoIt = true;
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
    SetBuildTableFlag(true);
    SetStartFromNullFlag(false);
    SetSecondaryParticle(G4Gamma::Gamma());
    if(!Model()) SetModel(new G4eeToTwoGammaModel());
    Model()->SetLowEnergyLimit(MinKinEnergy());
    Model()->SetHighEnergyLimit(MaxKinEnergy());
    AddEmModel(1, Model());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eplusAnnihilation::PrintInfo()
{} 

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

  fParticleChange.SetNumberOfSecondaries(2);

  G4double cosTeta = 2.*G4UniformRand()-1. , sinTeta = sqrt(1.-cosTeta*cosTeta);
  G4double phi     = twopi * G4UniformRand();
  G4ThreeVector direction (sinTeta*cos(phi), sinTeta*sin(phi), cosTeta);
  fParticleChange.AddSecondary( new G4DynamicParticle (G4Gamma::Gamma(),
                                            direction, electron_mass_c2) );
  fParticleChange.AddSecondary( new G4DynamicParticle (G4Gamma::Gamma(),
                                           -direction, electron_mass_c2) );
  // Kill the incident positron
  //
  fParticleChange.ProposeTrackStatus(fStopAndKill);
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
