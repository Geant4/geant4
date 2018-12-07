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
/// \file medical/fanoCavity/src/MyKleinNishinaCompton.cc
/// \brief Implementation of the MyKleinNishinaCompton class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "MyKleinNishinaCompton.hh"
#include "MyKleinNishinaMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4DataVector.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

MyKleinNishinaCompton::MyKleinNishinaCompton(DetectorConstruction* det,
                                             const G4ParticleDefinition*,
                                             const G4String& nam)
  :G4KleinNishinaCompton(0,nam), fDetector(det), fMessenger(0)
{
  fCrossSectionFactor = 1.;
  fMessenger = new MyKleinNishinaMessenger(this);    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MyKleinNishinaCompton::~MyKleinNishinaCompton()
{  
  delete fMessenger;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double MyKleinNishinaCompton::CrossSectionPerVolume(
                                       const G4Material* mat,
                                       const G4ParticleDefinition* part,
                                             G4double GammaEnergy,
                                             G4double, G4double)
{
  G4double xsection = G4VEmModel::CrossSectionPerVolume(mat,part,GammaEnergy);

  return xsection*fCrossSectionFactor;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MyKleinNishinaCompton::SampleSecondaries(
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
  G4double E0_m = gamEnergy0 / electron_mass_c2 ;

  G4ThreeVector gamDirection0 = aDynamicGamma->GetMomentumDirection();

  //
  // sample the energy rate of the scattered gamma 
  //

  G4double epsilon, epsilonsq, onecost, sint2, greject ;

  G4double eps0   = 1./(1. + 2.*E0_m);
  G4double eps0sq = eps0*eps0;
  G4double alpha1     = - log(eps0);
  G4double alpha2     = 0.5*(1.- eps0sq);

  do {
    if ( alpha1/(alpha1+alpha2) > G4UniformRand() ) {
      epsilon   = exp(-alpha1*G4UniformRand());   // eps0**r
      epsilonsq = epsilon*epsilon; 

    } else {
      epsilonsq = eps0sq + (1.- eps0sq)*G4UniformRand();
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
  G4double dirx = sinTeta*cos(Phi), diry = sinTeta*sin(Phi), dirz = cosTeta;

  //
  // update G4VParticleChange for the scattered gamma
  //
  // beam regeneration trick : restore incident beam
   
  G4ThreeVector gamDirection1 ( dirx,diry,dirz );
  gamDirection1.rotateUz(gamDirection0);
  G4double gamEnergy1 = epsilon*gamEnergy0;
  fParticleChange->SetProposedKineticEnergy(gamEnergy0);
  fParticleChange->ProposeMomentumDirection(gamDirection0);

  //
  // kinematic of the scattered electron
  //

  G4double eKinEnergy = gamEnergy0 - gamEnergy1;

  if(eKinEnergy > DBL_MIN) {
    G4ThreeVector eDirection 
                   = gamEnergy0*gamDirection0 - gamEnergy1*gamDirection1;
    eDirection = eDirection.unit();

    // create G4DynamicParticle object for the electron.
    G4DynamicParticle* dp 
                   = new G4DynamicParticle(theElectron,eDirection,eKinEnergy);
    fvect->push_back(dp);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

