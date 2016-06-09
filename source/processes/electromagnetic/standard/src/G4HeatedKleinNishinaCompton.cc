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
// $Id: G4HeatedKleinNishinaCompton.cc,v 1.5 2009-04-12 17:09:57 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
// Modifications:
// 
//
// Class Description:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include <CLHEP/Random/RandGamma.h>
#include "globals.hh"
#include "G4RandomDirection.hh"
#include "Randomize.hh"

#include "G4HeatedKleinNishinaCompton.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4DataVector.hh"
#include "G4ParticleChangeForGamma.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4HeatedKleinNishinaCompton::G4HeatedKleinNishinaCompton(const G4ParticleDefinition*,
                                             const G4String& nam)
  : G4VEmModel(nam)
{
  theGamma = G4Gamma::Gamma();
  theElectron = G4Electron::Electron();
  lowestGammaEnergy = 1.0*eV;
  fTemperature = 1.0*keV;
  fParticleChange = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4HeatedKleinNishinaCompton::~G4HeatedKleinNishinaCompton()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4HeatedKleinNishinaCompton::Initialise(const G4ParticleDefinition*,
                                       const G4DataVector&)
{
  if(!fParticleChange) fParticleChange = GetParticleChangeForGamma();
}

////////////////////////////////////////////////////////////////////////////
//
//

G4double G4HeatedKleinNishinaCompton::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                             G4double GammaEnergy,
                                             G4double Z, G4double,
                                             G4double, G4double)
{
  G4double CrossSection = 0.0 ;
  if ( Z < 0.9999 )                 return CrossSection;
  if ( GammaEnergy < 0.01*keV      ) return CrossSection;
  //  if ( GammaEnergy > (100.*GeV/Z) ) return CrossSection;

  static const G4double a = 20.0 , b = 230.0 , c = 440.0;
  
  static const G4double
    d1= 2.7965e-1*barn, d2=-1.8300e-1*barn, d3= 6.7527   *barn, d4=-1.9798e+1*barn,
    e1= 1.9756e-5*barn, e2=-1.0205e-2*barn, e3=-7.3913e-2*barn, e4= 2.7079e-2*barn,
    f1=-3.9178e-7*barn, f2= 6.8241e-5*barn, f3= 6.0480e-5*barn, f4= 3.0274e-4*barn;
     
  G4double p1Z = Z*(d1 + e1*Z + f1*Z*Z), p2Z = Z*(d2 + e2*Z + f2*Z*Z),
           p3Z = Z*(d3 + e3*Z + f3*Z*Z), p4Z = Z*(d4 + e4*Z + f4*Z*Z);

  G4double T0  = 15.0*keV; 
  if (Z < 1.5) T0 = 40.0*keV; 

  G4double X   = max(GammaEnergy, T0) / electron_mass_c2;
  CrossSection = p1Z*std::log(1.+2.*X)/X
               + (p2Z + p3Z*X + p4Z*X*X)/(1. + a*X + b*X*X + c*X*X*X);
		
  //  modification for low energy. (special case for Hydrogen)
  if (GammaEnergy < T0) {
    G4double dT0 = 1.*keV;
    X = (T0+dT0) / electron_mass_c2 ;
    G4double sigma = p1Z*log(1.+2*X)/X
                    + (p2Z + p3Z*X + p4Z*X*X)/(1. + a*X + b*X*X + c*X*X*X);
    G4double   c1 = -T0*(sigma-CrossSection)/(CrossSection*dT0);             
    G4double   c2 = 0.150; 
    if (Z > 1.5) c2 = 0.375-0.0556*log(Z);
    G4double    y = log(GammaEnergy/T0);
    CrossSection *= exp(-y*(c1+c2*y));          
  }
  //  G4cout << "e= " << GammaEnergy << " Z= " << Z << " cross= " << CrossSection << G4endl;
  return CrossSection;
}

//////////////////////////////////////////////////////////////////////////
//
//

void G4HeatedKleinNishinaCompton::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
					      const G4MaterialCutsCouple*,
					      const G4DynamicParticle* aDynamicGamma,
					      G4double,
					      G4double)
{
  // The scattered gamma energy is sampled according to Klein - Nishina formula.
  // The random number techniques of Butcher & Messel are used 
  // (Nuc Phys 20(1960),15).
  // Note : Effects due to binding of atomic electrons are negliged.

  // We start to prepare a heated electron from Maxwell distribution. 
  // Then we try to boost to the electron rest frame and make scattering.
  // The final step is to recover new gamma 4momentum in the lab frame

  G4double eMomentumC2   = CLHEP::RandGamma::shoot(1.5,1.);
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

  G4double epsilon0   = 1./(1. + 2.*E0_m);
  G4double epsilon0sq = epsilon0*epsilon0;
  G4double alpha1     = - log(epsilon0);
  G4double alpha2     = 0.5*(1.- epsilon0sq);

  do 
  {
    if ( alpha1/(alpha1+alpha2) > G4UniformRand() ) 
    {
      epsilon   = exp(-alpha1*G4UniformRand());   // epsilon0**r
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




  fParticleChange->SetProposedKineticEnergy(gamEnergy1);

  if( gamEnergy1 > lowestGammaEnergy ) 
  {
    gamDirection1 /= gamEnergy1;
    fParticleChange->ProposeMomentumDirection(gamDirection1);
  } 
  else 
  { 
    fParticleChange->ProposeTrackStatus(fStopAndKill);
    gamEnergy1 += fParticleChange->GetLocalEnergyDeposit();
    fParticleChange->ProposeLocalEnergyDeposit(gamEnergy1);
  }

  eKinEnergy = e4vfinal.t()-electron_mass_c2;

  if( eKinEnergy > DBL_MIN ) 
  {
    // create G4DynamicParticle object for the electron.
    eDirection = e4vfinal.vect();
    G4double eFinMomMag = eDirection.mag();
    eDirection /= eFinMomMag;
    G4DynamicParticle* dp = new G4DynamicParticle(theElectron,eDirection,eKinEnergy);
    fvect->push_back(dp);
  }
}

//////////////////////////////////////////////////////////////////////////


