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
// $Id: G4BetheHeitlerModel.cc,v 1.15 2010-10-25 19:02:32 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4BetheHeitlerModel
//
// Author:        Vladimir Ivanchenko on base of Michel Maire code
//
// Creation date: 15.03.2005
//
// Modifications:
// 18-04-05 Use G4ParticleChangeForGamma (V.Ivantchenko)
// 24-06-05 Increase number of bins to 200 (V.Ivantchenko)
// 16-11-05 replace shootBit() by G4UniformRand()  mma
// 04-12-05 SetProposedKineticEnergy(0.) for the killed photon (mma)
// 20-02-07 SelectRandomElement is called for any initial gamma energy 
//          in order to have selected element for polarized model (VI)
// 25-10-10 Removed unused table, added element selector (VI) 
//
// Class Description:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4BetheHeitlerModel.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4ParticleChangeForGamma.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4BetheHeitlerModel::G4BetheHeitlerModel(const G4ParticleDefinition*,
					 const G4String& nam)
  : G4VEmModel(nam)
{
  fParticleChange = 0;
  theGamma    = G4Gamma::Gamma();
  thePositron = G4Positron::Positron();
  theElectron = G4Electron::Electron();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4BetheHeitlerModel::~G4BetheHeitlerModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4BetheHeitlerModel::Initialise(const G4ParticleDefinition* p,
				     const G4DataVector& cuts)
{
  if(!fParticleChange) { fParticleChange = GetParticleChangeForGamma(); }
  InitialiseElementSelectors(p, cuts);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4BetheHeitlerModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
						G4double GammaEnergy, G4double Z,
						G4double, G4double, G4double)
// Calculates the microscopic cross section in GEANT4 internal units.
// A parametrized formula from L. Urban is used to estimate
// the total cross section.
// It gives a good description of the data from 1.5 MeV to 100 GeV.
// below 1.5 MeV: sigma=sigma(1.5MeV)*(GammaEnergy-2electronmass)
//                                   *(GammaEnergy-2electronmass) 
{
  static const G4double GammaEnergyLimit = 1.5*MeV;
  G4double CrossSection = 0.0 ;
  if ( Z < 0.9 || GammaEnergy <= 2.0*electron_mass_c2 ) { return CrossSection; }

  static const G4double
    a0= 8.7842e+2*microbarn, a1=-1.9625e+3*microbarn, a2= 1.2949e+3*microbarn,
    a3=-2.0028e+2*microbarn, a4= 1.2575e+1*microbarn, a5=-2.8333e-1*microbarn;

  static const G4double
    b0=-1.0342e+1*microbarn, b1= 1.7692e+1*microbarn, b2=-8.2381   *microbarn,
    b3= 1.3063   *microbarn, b4=-9.0815e-2*microbarn, b5= 2.3586e-3*microbarn;

  static const G4double
    c0=-4.5263e+2*microbarn, c1= 1.1161e+3*microbarn, c2=-8.6749e+2*microbarn,
    c3= 2.1773e+2*microbarn, c4=-2.0467e+1*microbarn, c5= 6.5372e-1*microbarn;

  G4double GammaEnergySave = GammaEnergy;
  if (GammaEnergy < GammaEnergyLimit) { GammaEnergy = GammaEnergyLimit; }

  G4double X=log(GammaEnergy/electron_mass_c2), X2=X*X, X3=X2*X, X4=X3*X, X5=X4*X;

  G4double F1 = a0 + a1*X + a2*X2 + a3*X3 + a4*X4 + a5*X5,
           F2 = b0 + b1*X + b2*X2 + b3*X3 + b4*X4 + b5*X5,
           F3 = c0 + c1*X + c2*X2 + c3*X3 + c4*X4 + c5*X5;     

  CrossSection = (Z + 1.)*(F1*Z + F2*Z*Z + F3);

  if (GammaEnergySave < GammaEnergyLimit) {

    X = (GammaEnergySave  - 2.*electron_mass_c2)
      / (GammaEnergyLimit - 2.*electron_mass_c2);
    CrossSection *= X*X;
  }

  if (CrossSection < 0.) { CrossSection = 0.; }
  return CrossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4BetheHeitlerModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
					    const G4MaterialCutsCouple* couple,
					    const G4DynamicParticle* aDynamicGamma,
					    G4double,
					    G4double)
// The secondaries e+e- energies are sampled using the Bethe - Heitler
// cross sections with Coulomb correction.
// A modified version of the random number techniques of Butcher & Messel
// is used (Nuc Phys 20(1960),15).
//
// GEANT4 internal units.
//
// Note 1 : Effects due to the breakdown of the Born approximation at
//          low energy are ignored.
// Note 2 : The differential cross section implicitly takes account of 
//          pair creation in both nuclear and atomic electron fields.
//          However triplet prodution is not generated.
{
  const G4Material* aMaterial = couple->GetMaterial();

  G4double GammaEnergy = aDynamicGamma->GetKineticEnergy();
  G4ParticleMomentum GammaDirection = aDynamicGamma->GetMomentumDirection();

  G4double epsil ;
  G4double epsil0 = electron_mass_c2/GammaEnergy ;
  if(epsil0 > 1.0) { return; }

  // do it fast if GammaEnergy < 2. MeV
  static const G4double Egsmall=2.*MeV;

  // select randomly one element constituing the material
  const G4Element* anElement = SelectRandomAtom(aMaterial, theGamma, GammaEnergy);

  if (GammaEnergy < Egsmall) {

    epsil = epsil0 + (0.5-epsil0)*G4UniformRand();

  } else {
    // now comes the case with GammaEnergy >= 2. MeV

    // Extract Coulomb factor for this Element
    G4double FZ = 8.*(anElement->GetIonisation()->GetlogZ3());
    if (GammaEnergy > 50.*MeV) { FZ += 8.*(anElement->GetfCoulomb()); }

    // limits of the screening variable
    G4double screenfac = 136.*epsil0/(anElement->GetIonisation()->GetZ3());
    G4double screenmax = exp ((42.24 - FZ)/8.368) - 0.952 ;
    G4double screenmin = min(4.*screenfac,screenmax);

    // limits of the energy sampling
    G4double epsil1 = 0.5 - 0.5*sqrt(1. - screenmin/screenmax) ;
    G4double epsilmin = max(epsil0,epsil1) , epsilrange = 0.5 - epsilmin;

    //
    // sample the energy rate of the created electron (or positron)
    //
    //G4double epsil, screenvar, greject ;
    G4double  screenvar, greject ;

    G4double F10 = ScreenFunction1(screenmin) - FZ;
    G4double F20 = ScreenFunction2(screenmin) - FZ;
    G4double NormF1 = max(F10*epsilrange*epsilrange,0.); 
    G4double NormF2 = max(1.5*F20,0.);

    do {
      if ( NormF1/(NormF1+NormF2) > G4UniformRand() ) {
	epsil = 0.5 - epsilrange*pow(G4UniformRand(), 0.333333);
	screenvar = screenfac/(epsil*(1-epsil));
	greject = (ScreenFunction1(screenvar) - FZ)/F10;
              
      } else { 
	epsil = epsilmin + epsilrange*G4UniformRand();
	screenvar = screenfac/(epsil*(1-epsil));
	greject = (ScreenFunction2(screenvar) - FZ)/F20;
      }

    } while( greject < G4UniformRand() );

  }   //  end of epsil sampling
   
  //
  // fixe charges randomly
  //

  G4double ElectTotEnergy, PositTotEnergy;
  if (G4UniformRand() > 0.5) {

    ElectTotEnergy = (1.-epsil)*GammaEnergy;
    PositTotEnergy = epsil*GammaEnergy;
     
  } else {
    
    PositTotEnergy = (1.-epsil)*GammaEnergy;
    ElectTotEnergy = epsil*GammaEnergy;
  }

  //
  // scattered electron (positron) angles. ( Z - axis along the parent photon)
  //
  //  universal distribution suggested by L. Urban 
  // (Geant3 manual (1993) Phys211),
  //  derived from Tsai distribution (Rev Mod Phys 49,421(1977))

  G4double u;
  const G4double a1 = 0.625 , a2 = 3.*a1 , d = 27. ;

  if (9./(9.+d) >G4UniformRand()) u= - log(G4UniformRand()*G4UniformRand())/a1;
  else                            u= - log(G4UniformRand()*G4UniformRand())/a2;

  G4double TetEl = u*electron_mass_c2/ElectTotEnergy;
  G4double TetPo = u*electron_mass_c2/PositTotEnergy;
  G4double Phi  = twopi * G4UniformRand();
  G4double dxEl= sin(TetEl)*cos(Phi),dyEl= sin(TetEl)*sin(Phi),dzEl=cos(TetEl);
  G4double dxPo=-sin(TetPo)*cos(Phi),dyPo=-sin(TetPo)*sin(Phi),dzPo=cos(TetPo);
   
  //
  // kinematic of the created pair
  //
  // the electron and positron are assumed to have a symetric
  // angular distribution with respect to the Z axis along the parent photon.

  G4double ElectKineEnergy = max(0.,ElectTotEnergy - electron_mass_c2);

  G4ThreeVector ElectDirection (dxEl, dyEl, dzEl);
  ElectDirection.rotateUz(GammaDirection);   

  // create G4DynamicParticle object for the particle1  
  G4DynamicParticle* aParticle1= new G4DynamicParticle(
		     theElectron,ElectDirection,ElectKineEnergy);
  
  // the e+ is always created (even with Ekine=0) for further annihilation.

  G4double PositKineEnergy = max(0.,PositTotEnergy - electron_mass_c2);

  G4ThreeVector PositDirection (dxPo, dyPo, dzPo);
  PositDirection.rotateUz(GammaDirection);   

  // create G4DynamicParticle object for the particle2 
  G4DynamicParticle* aParticle2= new G4DynamicParticle(
                      thePositron,PositDirection,PositKineEnergy);

  // Fill output vector
  fvect->push_back(aParticle1);
  fvect->push_back(aParticle2);

  // kill incident photon
  fParticleChange->SetProposedKineticEnergy(0.);
  fParticleChange->ProposeTrackStatus(fStopAndKill);   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
