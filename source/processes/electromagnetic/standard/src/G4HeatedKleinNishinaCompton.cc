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
// $Id: G4HeatedKleinNishinaCompton.cc 91726 2015-08-03 15:41:36Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4HeatedKleinNishinaCompton
//
// Author:        Vladimir Grichine on base of M. Maire and V. Ivanchenko code
//
// Creation date: 15.03.2009
//
// Modifications: 07.07.2014 V.Ivanchenko make direct inheritence from 
//                           G4KleinNishinaCompton 
// 
//
// Class Description:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4HeatedKleinNishinaCompton.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4RandomDirection.hh"
#include "Randomize.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4ParticleChangeForGamma.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4HeatedKleinNishinaCompton::G4HeatedKleinNishinaCompton(
  const G4ParticleDefinition* p, const G4String& nam)
  : G4KleinNishinaCompton(p, nam)
{
  fTemperature = 1.0*keV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4HeatedKleinNishinaCompton::~G4HeatedKleinNishinaCompton()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4HeatedKleinNishinaCompton::SampleSecondaries(
                         std::vector<G4DynamicParticle*>* fvect,
			 const G4MaterialCutsCouple*,
			 const G4DynamicParticle* aDynamicGamma,
			 G4double, G4double)
{
  // do nothing below the threshold
  if(aDynamicGamma->GetKineticEnergy() <= LowEnergyLimit()) { return; }

  // The scattered gamma energy is sampled according to Klein - Nishina formula.
  // The random number techniques of Butcher & Messel are used 
  // (Nuc Phys 20(1960),15).
  // Note : Effects due to binding of atomic electrons are negliged.

  // We start to prepare a heated electron from Maxwell distribution. 
  // Then we try to boost to the electron rest frame and make scattering.
  // The final step is to recover new gamma 4momentum in the lab frame

  G4double eMomentumC2  = G4RandGamma::shoot(1.5, 1.);
  eMomentumC2          *= 2*electron_mass_c2*fTemperature; // electron (pc)^2
  G4ThreeVector eMomDir = G4RandomDirection();
  eMomDir              *= std::sqrt(eMomentumC2);
  G4double eEnergy      = std::sqrt(eMomentumC2+electron_mass_c2*electron_mass_c2);
  G4LorentzVector electron4v = G4LorentzVector(eMomDir,eEnergy);
  G4ThreeVector bst = electron4v.boostVector();

  G4LorentzVector gamma4v = aDynamicGamma->Get4Momentum();
  gamma4v.boost(-bst);

  G4ThreeVector gammaMomV = gamma4v.vect();
  G4double  gamEnergy0    = gammaMomV.mag();
 
 
  // G4double gamEnergy0 = aDynamicGamma->GetKineticEnergy();
  G4double E0_m = gamEnergy0 / electron_mass_c2 ;

  // G4ThreeVector gamDirection0 = /aDynamicGamma->GetMomentumDirection();

  G4ThreeVector gamDirection0 = gammaMomV/gamEnergy0;
  
  // sample the energy rate of the scattered gamma in the electron rest frame
  //

  G4double epsilon, epsilonsq, onecost, sint2, greject ;

  G4double eps0       = 1./(1. + 2.*E0_m);
  G4double epsilon0sq = eps0*eps0;
  G4double alpha1     = - G4Log(eps0);
  G4double alpha2     = 0.5*(1.- epsilon0sq);

  G4int nloop = 0;
  do 
  {
    ++nloop;
    // false interaction if too many iterations
    if(nloop > 1000) { return; }

    if ( alpha1/(alpha1+alpha2) > G4UniformRand() ) 
    {
      epsilon   = G4Exp(-alpha1*G4UniformRand());   // eps0**r
      epsilonsq = epsilon*epsilon; 

    } 
    else 
    {
      epsilonsq = epsilon0sq + (1.- epsilon0sq)*G4UniformRand();
      epsilon   = sqrt(epsilonsq);
    };

    onecost = (1.- epsilon)/(epsilon*E0_m);
    sint2   = onecost*(2.-onecost);
    greject = 1. - epsilon*sint2/(1.+ epsilonsq);

    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while (greject < G4UniformRand());
 
  //
  // scattered gamma angles. ( Z - axis along the parent gamma)
  //

  G4double cosTeta = 1. - onecost; 
  G4double sinTeta = sqrt (sint2);
  G4double Phi     = twopi * G4UniformRand();
  G4double dirx    = sinTeta*cos(Phi), diry = sinTeta*sin(Phi), dirz = cosTeta;

  //
  // update G4VParticleChange for the scattered gamma
  //
   
  G4ThreeVector gamDirection1 ( dirx,diry,dirz );
  gamDirection1.rotateUz(gamDirection0);
  G4double gamEnergy1  = epsilon*gamEnergy0;
  gamDirection1       *= gamEnergy1;

  G4LorentzVector gamma4vfinal = G4LorentzVector(gamDirection1,gamEnergy1);

  
  // kinematic of the scattered electron
  //

  G4double eKinEnergy = gamEnergy0 - gamEnergy1;
  G4ThreeVector eDirection = gamEnergy0*gamDirection0 - gamEnergy1*gamDirection1;
  eDirection = eDirection.unit();
  G4double eFinalMom = std::sqrt(eKinEnergy*(eKinEnergy+2*electron_mass_c2));
  eDirection *= eFinalMom;
  G4LorentzVector e4vfinal = G4LorentzVector(eDirection,gamEnergy1+electron_mass_c2);
  
  gamma4vfinal.boost(bst);
  e4vfinal.boost(bst);

  gamDirection1 = gamma4vfinal.vect();
  gamEnergy1 = gamDirection1.mag(); 
  gamDirection1 /= gamEnergy1;

  G4double edep = 0.0;
  if(gamEnergy1 > lowestSecondaryEnergy) {
    fParticleChange->ProposeMomentumDirection(gamDirection1);
    fParticleChange->SetProposedKineticEnergy(gamEnergy1);
  } else { 
    fParticleChange->ProposeTrackStatus(fStopAndKill);
    fParticleChange->SetProposedKineticEnergy(0.0);
    edep = gamEnergy1;
  }

  //
  // kinematic of the scattered electron
  //
  eKinEnergy = e4vfinal.t()-electron_mass_c2;

  if(eKinEnergy > lowestSecondaryEnergy) {

    eDirection = e4vfinal.vect().unit();

    // create G4DynamicParticle object for the electron.
    G4DynamicParticle* dp = 
      new G4DynamicParticle(theElectron,eDirection,eKinEnergy);
    fvect->push_back(dp);
  } else {
    edep += eKinEnergy;  
  }
  // energy balance
  if(edep > 0.0) { 
    fParticleChange->ProposeLocalEnergyDeposit(edep);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


