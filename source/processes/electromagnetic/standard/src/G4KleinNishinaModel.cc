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
// $Id: G4KleinNishinaModel.cc,v 1.5 2011-01-02 19:58:54 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4KleinNishinaModel
//
// Author:        Vladimir Ivanchenko on base of G4KleinNishinaCompton
//
// Creation date: 13.06.2010
//
// Modifications:
//
// Class Description:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4KleinNishinaModel.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4DataVector.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4LossTableManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4KleinNishinaModel::G4KleinNishinaModel(const G4String& nam)
  : G4VEmModel(nam),isInitialized(false)
{
  theGamma = G4Gamma::Gamma();
  theElectron = G4Electron::Electron();
  lowestGammaEnergy = 1.0*eV;
  fProbabilities.resize(9,0.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4KleinNishinaModel::~G4KleinNishinaModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4KleinNishinaModel::Initialise(const G4ParticleDefinition* p,
				     const G4DataVector& cuts)
{
  fAtomDeexcitation = G4LossTableManager::Instance()->AtomDeexcitation();
  InitialiseElementSelectors(p, cuts);

  if (isInitialized) { return; }
  fParticleChange = GetParticleChangeForGamma();
  isInitialized = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4KleinNishinaModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
						G4double GammaEnergy,
						G4double Z, G4double,
						G4double, G4double)
{
  G4double CrossSection = 0.0 ;
  if ( Z < 0.9999 || GammaEnergy < 0.1*keV) { return CrossSection; }

  static const G4double a = 20.0 , b = 230.0 , c = 440.0;
  
  static const G4double
    d1= 2.7965e-1*barn, d2=-1.8300e-1*barn, d3= 6.7527   *barn, d4=-1.9798e+1*barn,
    e1= 1.9756e-5*barn, e2=-1.0205e-2*barn, e3=-7.3913e-2*barn, e4= 2.7079e-2*barn,
    f1=-3.9178e-7*barn, f2= 6.8241e-5*barn, f3= 6.0480e-5*barn, f4= 3.0274e-4*barn;
     
  G4double p1Z = Z*(d1 + e1*Z + f1*Z*Z), p2Z = Z*(d2 + e2*Z + f2*Z*Z),
           p3Z = Z*(d3 + e3*Z + f3*Z*Z), p4Z = Z*(d4 + e4*Z + f4*Z*Z);

  G4double T0  = 15.0*keV; 
  if (Z < 1.5) { T0 = 40.0*keV; } 

  G4double X   = max(GammaEnergy, T0) / electron_mass_c2;
  CrossSection = p1Z*std::log(1.+2.*X)/X
               + (p2Z + p3Z*X + p4Z*X*X)/(1. + a*X + b*X*X + c*X*X*X);
		
  //  modification for low energy. (special case for Hydrogen)
  if (GammaEnergy < T0) {
    G4double dT0 = keV;
    X = (T0+dT0) / electron_mass_c2 ;
    G4double sigma = p1Z*log(1.+2*X)/X
                    + (p2Z + p3Z*X + p4Z*X*X)/(1. + a*X + b*X*X + c*X*X*X);
    G4double   c1 = -T0*(sigma-CrossSection)/(CrossSection*dT0);             
    G4double   c2 = 0.150; 
    if (Z > 1.5) { c2 = 0.375-0.0556*log(Z); }
    G4double    y = log(GammaEnergy/T0);
    CrossSection *= exp(-y*(c1+c2*y));          
  }
  //  G4cout << "e= " << GammaEnergy << " Z= " << Z 
  //  << " cross= " << CrossSection << G4endl;
  return CrossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4KleinNishinaModel::SampleSecondaries(
			     std::vector<G4DynamicParticle*>* fvect,
			     const G4MaterialCutsCouple* couple,
			     const G4DynamicParticle* aDynamicGamma,
			     G4double,
			     G4double)
{
  // primary gamma
  G4double energy = aDynamicGamma->GetKineticEnergy();
  G4ThreeVector direction = aDynamicGamma->GetMomentumDirection();

  // select atom
  const G4Element* elm = SelectRandomAtom(couple, theGamma, energy);

  // select shell first
  G4int Z = (G4int)elm->GetZ();
  G4int nShells = elm->GetNbOfAtomicShells();
  if(nShells > (G4int)fProbabilities.size()) { fProbabilities.resize(nShells); }
  G4double totprob = 0.0;
  G4int i = 0;
  for(; i<nShells; ++i) {
    G4double prob = 0.0;
    if(energy > elm->GetAtomicShell(i)) { 
      prob = (G4double)elm->GetNbOfShellElectrons(i);
    }
    totprob += prob;
    fProbabilities[i] = totprob; 
  }
  if(totprob == 0.0) { return; }

  G4LorentzVector lv1, lv2, lv3;
  G4double eKinEnergy = 0.0;
  G4double gamEnergy1 = 0.0;

  // Loop on sampling
  do {
    G4double xprob = totprob*G4UniformRand();

    // select shell
    for(i=0; i<nShells; ++i) { if(xprob <= fProbabilities[i]) {break;} }
    if( i == nShells ) { return; }
   
    // rest frame of the electron
    G4double bindingEnergy = elm->GetAtomicShell(i);
    G4double eEnergy = bindingEnergy + electron_mass_c2;
    G4double eTotMomentum = sqrt(bindingEnergy*(bindingEnergy + 2*electron_mass_c2));
    G4ThreeVector eDir = G4RandomDirection();
    lv1.set(0.0,0.0,energy,energy);
    lv2.set(eTotMomentum*eDir.x(),eTotMomentum*eDir.y(),
	    eTotMomentum*eDir.z(),eEnergy);
    G4ThreeVector bst = lv2.boostVector();
    lv1.boost(-bst);

    // In the rest frame of the electron
    // The scattered gamma energy is sampled according to Klein - Nishina formula.
    // The random number techniques of Butcher & Messel are used 
    // (Nuc Phys 20(1960),15).
 
    G4double gamEnergy0 = lv1.e();
    G4double E0_m = gamEnergy0 / electron_mass_c2 ;

    //
    // sample the energy rate of the scattered gamma 
    //

    G4double epsilon, epsilonsq, onecost, sint2, greject ;

    G4double epsilon0   = 1./(1. + 2.*E0_m);
    G4double epsilon0sq = epsilon0*epsilon0;
    G4double alpha1     = - log(epsilon0);
    G4double alpha2     = 0.5*(1.- epsilon0sq);

    do {
      if ( alpha1/(alpha1+alpha2) > G4UniformRand() ) {
	epsilon   = exp(-alpha1*G4UniformRand());   // epsilon0**r
	epsilonsq = epsilon*epsilon; 

      } else {
	epsilonsq = epsilon0sq + (1.- epsilon0sq)*G4UniformRand();
	epsilon   = sqrt(epsilonsq);
      };

      onecost = (1.- epsilon)/(epsilon*E0_m);
      sint2   = onecost*(2.-onecost);
      greject = 1. - epsilon*sint2/(1.+ epsilonsq);

    } while (greject < G4UniformRand());
    gamEnergy1 = epsilon*gamEnergy0;
 
    //
    // scattered gamma angles. ( Z - axis along the parent gamma)
    //

    G4double cosTeta = 1. - onecost; 
    G4double sinTeta = sqrt (sint2);
    G4double Phi  = twopi * G4UniformRand();
    G4double dirx = sinTeta*cos(Phi);
    G4double diry = sinTeta*sin(Phi); 

    // before scattering
    lv2.set(0.0,0.0,0.0,electron_mass_c2);
    lv2 += lv1;
 
    // after scattering
    lv1.set(gamEnergy1*dirx,gamEnergy1*diry,gamEnergy1*cosTeta,gamEnergy1);
    lv2 -= lv1;
    lv2.boost(bst);
    lv1.boost(bst);
    eKinEnergy = lv2.e() - electron_mass_c2 - 2*bindingEnergy;
  } while ( eKinEnergy < 0.0 );

  // gamma kinematics
  gamEnergy1 = lv1.e();

  //
  // update G4VParticleChange for the scattered gamma
  //
   
  if(gamEnergy1 > lowestGammaEnergy) {
    G4ThreeVector gamDirection1 = lv1.vect().unit();
    gamDirection1.rotateUz(direction);
    fParticleChange->ProposeMomentumDirection(gamDirection1);
  } else { 
    fParticleChange->ProposeTrackStatus(fStopAndKill);
    gamEnergy1 = 0.0;
  }
  fParticleChange->SetProposedKineticEnergy(gamEnergy1);
  //
  // kinematic of the scattered electron
  //
  G4ThreeVector eDirection = (lv2.vect()).unit();
  eDirection.rotateUz(direction);
  G4DynamicParticle* dp = new G4DynamicParticle(theElectron,eDirection,eKinEnergy);
  fvect->push_back(dp);

  G4double edep = energy - gamEnergy1 - eKinEnergy;
  
  // sample deexcitation
  //
  if(fAtomDeexcitation) {
    G4int index = couple->GetIndex();
    if(fAtomDeexcitation->CheckDeexcitationActiveRegion(index)) {
      G4AtomicShellEnumerator as = G4AtomicShellEnumerator(i);
      const G4AtomicShell* shell = fAtomDeexcitation->GetAtomicShell(Z, as);    
      size_t nbefore = fvect->size();
      fAtomDeexcitation->GenerateParticles(fvect, shell, Z, index);
      size_t nafter = fvect->size();
      if(nafter > nbefore) {
	for (size_t j=nbefore; j<nafter; ++j) {
	  edep -= ((*fvect)[j])->GetKineticEnergy();
	} 
      }
    }
  }
  // energy balance
  if(edep < 0.0) { edep = 0.0; }
  fParticleChange->ProposeLocalEnergyDeposit(edep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

