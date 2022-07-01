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
// File name:     G4KleinNishinaCompton
//
// Author:        Vladimir Ivanchenko on base of Michel Maire code
//
// Creation date: 15.03.2005
//
// Modifications:
// 18-04-05 Use G4ParticleChangeForGamma (V.Ivantchenko)
// 27-03-06 Remove upper limit of cross section (V.Ivantchenko)
//
// Class Description:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4KleinNishinaCompton.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4DataVector.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4KleinNishinaCompton::G4KleinNishinaCompton(const G4ParticleDefinition*,
                                             const G4String& nam)
  : G4VEmModel(nam)
{
  theGamma = G4Gamma::Gamma();
  theElectron = G4Electron::Electron();
  lowestSecondaryEnergy = 100.0*eV;
  fParticleChange = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4KleinNishinaCompton::~G4KleinNishinaCompton() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4KleinNishinaCompton::Initialise(const G4ParticleDefinition* p,
                                       const G4DataVector& cuts)
{
  if(IsMaster()) { InitialiseElementSelectors(p, cuts); }
  if(nullptr == fParticleChange) { 
    fParticleChange = GetParticleChangeForGamma(); 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4KleinNishinaCompton::InitialiseLocal(const G4ParticleDefinition*,
                                            G4VEmModel* masterModel)
{
  SetElementSelectors(masterModel->GetElementSelectors());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4KleinNishinaCompton::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                             G4double GammaEnergy,
                                             G4double Z, G4double,
                                             G4double, G4double)
{
  G4double xSection = 0.0 ;
  if (GammaEnergy <= LowEnergyLimit()) { return xSection; }

  static const G4double a = 20.0 , b = 230.0 , c = 440.0;

  static const G4double
  d1= 2.7965e-1*CLHEP::barn, d2=-1.8300e-1*CLHEP::barn, 
  d3= 6.7527   *CLHEP::barn, d4=-1.9798e+1*CLHEP::barn,
  e1= 1.9756e-5*CLHEP::barn, e2=-1.0205e-2*CLHEP::barn, 
  e3=-7.3913e-2*CLHEP::barn, e4= 2.7079e-2*CLHEP::barn,
  f1=-3.9178e-7*CLHEP::barn, f2= 6.8241e-5*CLHEP::barn, 
  f3= 6.0480e-5*CLHEP::barn, f4= 3.0274e-4*CLHEP::barn;
       
  G4double p1Z = Z*(d1 + e1*Z + f1*Z*Z), p2Z = Z*(d2 + e2*Z + f2*Z*Z),
           p3Z = Z*(d3 + e3*Z + f3*Z*Z), p4Z = Z*(d4 + e4*Z + f4*Z*Z);

  G4double T0  = 15.0*keV; 
  if (Z < 1.5) { T0 = 40.0*keV; }

  G4double X   = max(GammaEnergy, T0) / electron_mass_c2;
  xSection = p1Z*G4Log(1.+2.*X)/X
               + (p2Z + p3Z*X + p4Z*X*X)/(1. + a*X + b*X*X + c*X*X*X);
                
  //  modification for low energy. (special case for Hydrogen)
  if (GammaEnergy < T0) {
    static const G4double dT0 = keV;
    X = (T0+dT0) / electron_mass_c2 ;
    G4double sigma = p1Z*G4Log(1.+2*X)/X
                    + (p2Z + p3Z*X + p4Z*X*X)/(1. + a*X + b*X*X + c*X*X*X);
    G4double   c1 = -T0*(sigma-xSection)/(xSection*dT0);             
    G4double   c2 = 0.150; 
    if (Z > 1.5) { c2 = 0.375-0.0556*G4Log(Z); }
    G4double    y = G4Log(GammaEnergy/T0);
    xSection *= G4Exp(-y*(c1+c2*y));          
  }
  // G4cout<<"e= "<< GammaEnergy<<" Z= "<<Z<<" cross= " << xSection << G4endl;
  return xSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4KleinNishinaCompton::SampleSecondaries(
                            std::vector<G4DynamicParticle*>* fvect,
                            const G4MaterialCutsCouple*,
                            const G4DynamicParticle* aDynamicGamma,
                            G4double,
                            G4double)
{
  // The scattered gamma energy is sampled according to Klein - Nishina formula.
  // The random number techniques of Butcher & Messel are used 
  // (Nuc Phys 20(1960),15).
  // Note : Effects due to binding of atomic electrons are negliged.
 
  G4double gamEnergy0 = aDynamicGamma->GetKineticEnergy();

  // do nothing below the threshold
  if(gamEnergy0 <= LowEnergyLimit()) { return; }

  G4double E0_m = gamEnergy0 / electron_mass_c2 ;

  G4ThreeVector gamDirection0 = aDynamicGamma->GetMomentumDirection();

  //
  // sample the energy rate of the scattered gamma 
  //

  G4double epsilon, epsilonsq, onecost, sint2, greject ;

  G4double eps0       = 1./(1. + 2.*E0_m);
  G4double epsilon0sq = eps0*eps0;
  G4double alpha1     = - G4Log(eps0);
  G4double alpha2     = alpha1 + 0.5*(1.- epsilon0sq);

  CLHEP::HepRandomEngine* rndmEngineMod = G4Random::getTheEngine();
  G4double rndm[3];

  static const G4int nlooplim = 1000;
  G4int nloop = 0;
  do {
    ++nloop;
    // false interaction if too many iterations
    if(nloop > nlooplim) { return; }

    // 3 random numbers to sample scattering
    rndmEngineMod->flatArray(3, rndm);

    if ( alpha1 > alpha2*rndm[0] ) {
      epsilon   = G4Exp(-alpha1*rndm[1]);   // eps0**r
      epsilonsq = epsilon*epsilon; 

    } else {
      epsilonsq = epsilon0sq + (1.- epsilon0sq)*rndm[1];
      epsilon   = sqrt(epsilonsq);
    };

    onecost = (1.- epsilon)/(epsilon*E0_m);
    sint2   = onecost*(2.-onecost);
    greject = 1. - epsilon*sint2/(1.+ epsilonsq);

    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while (greject < rndm[2]);
 
  //
  // scattered gamma angles. ( Z - axis along the parent gamma)
  //

  if(sint2 < 0.0) { sint2 = 0.0; }
  G4double cosTeta = 1. - onecost; 
  G4double sinTeta = sqrt (sint2);
  G4double Phi     = twopi * rndmEngineMod->flat();

  //
  // update G4VParticleChange for the scattered gamma
  //
   
  G4ThreeVector gamDirection1(sinTeta*cos(Phi), sinTeta*sin(Phi), cosTeta);
  gamDirection1.rotateUz(gamDirection0);
  G4double gamEnergy1 = epsilon*gamEnergy0;
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

  G4double eKinEnergy = gamEnergy0 - gamEnergy1;

  if(eKinEnergy > lowestSecondaryEnergy) {
    G4ThreeVector eDirection = gamEnergy0*gamDirection0 - gamEnergy1*gamDirection1;
    eDirection = eDirection.unit();

    // create G4DynamicParticle object for the electron.
    auto dp = new G4DynamicParticle(theElectron,eDirection,eKinEnergy);
    fvect->push_back(dp);
  } else {
    edep += eKinEnergy;  
  }
  // energy balance
  if(edep > 0.0) { 
    fParticleChange->ProposeLocalEnergyDeposit(edep);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


