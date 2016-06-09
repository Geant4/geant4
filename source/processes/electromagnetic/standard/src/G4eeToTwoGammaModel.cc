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
// $Id: G4eeToTwoGammaModel.cc,v 1.8 2005/04/29 18:02:35 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:   G4eeToTwoGammaModel
//
// Author:        Vladimir Ivanchenko on base of Michel Maire code
//
// Creation date: 02.08.2004
//
// Modifications:
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 18-04-05 Compute CrossSectionPerVolume (V.Ivantchenko)
//
//
// Class Description:
//
// Implementation of e+ annihilation into 2 gamma
//
// The secondaries Gamma energies are sampled using the Heitler cross section.
//
// A modified version of the random number techniques of Butcher & Messel
// is used (Nuc Phys 20(1960),15).
//
// GEANT4 internal units.
//
// Note 1: The initial electron is assumed free and at rest.
//
// Note 2: The annihilation processes producing one or more than two photons are
//         ignored, as negligible compared to the two photons process.



//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eeToTwoGammaModel.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4eeToTwoGammaModel::G4eeToTwoGammaModel(const G4ParticleDefinition*,
                                         const G4String& nam)
  : G4VEmModel(nam),
  pi_rcl2(pi*classic_electr_radius*classic_electr_radius)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eeToTwoGammaModel::~G4eeToTwoGammaModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToTwoGammaModel::Initialise(const G4ParticleDefinition*,
                                     const G4DataVector&)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eeToTwoGammaModel::CrossSectionPerVolume(const G4Material* material,
						    const G4ParticleDefinition*,
						    G4double kineticEnergy,
						    G4double,
						    G4double)
{
  // Calculates the cross section per atom of annihilation into two photons
  // from the Heilter formula.
  G4double eDensity = material->GetElectronDensity();

  G4double tau   = kineticEnergy/electron_mass_c2;
  G4double gam   = tau + 1.0;
  G4double gamma2= gam*gam;
  G4double bg2   = tau * (tau+2.0);
  G4double bg    = sqrt(bg2);

  G4double cross = pi_rcl2*eDensity*((gamma2+4*gam+1.)*log(gam+bg) - (gam+3.)*bg)
                 / (bg2*(gam+1.));
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

vector<G4DynamicParticle*>* G4eeToTwoGammaModel::SampleSecondaries(
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle* dp,
                                   G4double,
                                   G4double)
{
  G4double PositKinEnergy = dp->GetKineticEnergy();
  G4ThreeVector PositDirection = dp->GetMomentumDirection();

  G4double tau     = PositKinEnergy/electron_mass_c2;
  G4double gam     = tau + 1.0;
  G4double tau2    = tau + 2.0;
  G4double sqgrate = sqrt(tau/tau2)*0.5;
  G4double sqg2m1  = sqrt(tau*tau2);

  // limits of the energy sampling
  G4double epsilmin = 0.5 - sqgrate;
  G4double epsilmax = 0.5 + sqgrate;
  G4double epsilqot = epsilmax/epsilmin;

  //
  // sample the energy rate of the created gammas
  //
  G4double epsil, greject;

  do {
     epsil = epsilmin*pow(epsilqot,G4UniformRand());
     greject = 1. - epsil + (2.*gam*epsil-1.)/(epsil*tau2*tau2);
  } while( greject < G4UniformRand() );

  //
  // scattered Gamma angles. ( Z - axis along the parent positron)
  //

  G4double cost = (epsil*tau2-1.)/(epsil*sqg2m1);
  G4double sint = sqrt((1.+cost)*(1.-cost));
  G4double phi  = twopi * G4UniformRand();

  G4double dirx = sint*cos(phi) , diry = sint*sin(phi) , dirz = cost;

  //
  // kinematic of the created pair
  //

  G4double TotalAvailableEnergy = PositKinEnergy + 2.0*electron_mass_c2;
  G4double Phot1Energy = epsil*TotalAvailableEnergy;

  vector<G4DynamicParticle*>* vdp = new vector<G4DynamicParticle*>;

  G4ThreeVector Phot1Direction (dirx, diry, dirz);
  Phot1Direction.rotateUz(PositDirection);
  G4DynamicParticle* aParticle1 = new G4DynamicParticle (G4Gamma::Gamma(),
                                                 Phot1Direction, Phot1Energy);
  vdp->push_back(aParticle1);

  G4double Phot2Energy =(1.-epsil)*TotalAvailableEnergy;
  G4double Eratio= Phot1Energy/Phot2Energy;
  G4double PositP= sqrt(PositKinEnergy*(PositKinEnergy+2.*electron_mass_c2));
  G4ThreeVector Phot2Direction (-dirx*Eratio, -diry*Eratio,
                                    (PositP-dirz*Phot1Energy)/Phot2Energy);
  Phot2Direction.unit();
  Phot2Direction.rotateUz(PositDirection);
  // create G4DynamicParticle object for the particle2
  G4DynamicParticle* aParticle2= new G4DynamicParticle (G4Gamma::Gamma(),
                                                 Phot2Direction, Phot2Energy);
  vdp->push_back(aParticle2);
  return vdp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
