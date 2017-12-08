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
// $Id: G4eplusAnnihilation.cc 107058 2017-11-01 14:54:12Z gcosmo $
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
// 30-05-12 propagate parent weight to secondaries (D. Sawkey)
//

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

void G4eplusAnnihilation::StreamProcessInfo(std::ostream&,
                                            G4String) const
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
  CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();

  G4double cosTeta = 2.*rndmEngine->flat()-1.; 
  G4double sinTeta = sqrt((1.-cosTeta)*(1.0 + cosTeta));
  G4double phi     = twopi * rndmEngine->flat();
  G4ThreeVector dir(sinTeta*cos(phi), sinTeta*sin(phi), cosTeta);
  phi = twopi * rndmEngine->flat();
  G4double cosphi = cos(phi);
  G4double sinphi = sin(phi);
  G4ThreeVector pol(cosphi, sinphi, 0.0);
  pol.rotateUz(dir);

  // e+ parameters
  G4double weight = aTrack.GetWeight();
  G4double time   = aTrack.GetGlobalTime();

  // add gammas
  fParticleChange.SetNumberOfSecondaries(2);
  G4DynamicParticle* dp = 
    new G4DynamicParticle(theGamma, dir, electron_mass_c2);
  dp->SetPolarization(pol.x(),pol.y(),pol.z());
  G4Track* track = new G4Track(dp, time, aTrack.GetPosition());
  track->SetTouchableHandle(aTrack.GetTouchableHandle());
  track->SetWeight(weight); 
  pParticleChange->AddSecondary(track);

  dp = new G4DynamicParticle(theGamma,-dir, electron_mass_c2);
  pol.set(-sinphi, cosphi, 0.0);
  pol.rotateUz(dir);
  dp->SetPolarization(pol.x(),pol.y(),pol.z());
  track = new G4Track(dp, time, aTrack.GetPosition());
  track->SetTouchableHandle(aTrack.GetTouchableHandle());
  track->SetWeight(weight); 
  pParticleChange->AddSecondary(track);

  // Kill the incident positron
  //
  fParticleChange.ProposeTrackStatus(fStopAndKill);
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eplusAnnihilation::ProcessDescription(std::ostream& out) const
{
  out << "<strong>Positron annihilation</strong>";
  G4VEmProcess::ProcessDescription(out);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
