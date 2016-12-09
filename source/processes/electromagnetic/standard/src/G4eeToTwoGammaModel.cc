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
// $Id: G4eeToTwoGammaModel.cc 101198 2016-11-09 09:34:52Z gcosmo $
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
// 08-04-05 Major optimisation of internal interfaces (V.Ivanchenko)
// 18-04-05 Compute CrossSectionPerVolume (V.Ivanchenko)
// 06-02-06 ComputeCrossSectionPerElectron, ComputeCrossSectionPerAtom (mma)
// 29-06-06 Fix problem for zero energy incident positron (V.Ivanchenko) 
// 20-10-06 Add theGamma as a member (V.Ivanchenko)
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
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4TrackStatus.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4eeToTwoGammaModel::G4eeToTwoGammaModel(const G4ParticleDefinition*,
                                         const G4String& nam)
  : G4VEmModel(nam),
    pi_rcl2(pi*classic_electr_radius*classic_electr_radius),
    isInitialised(false)
{
  theGamma = G4Gamma::Gamma();
  fParticleChange = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eeToTwoGammaModel::~G4eeToTwoGammaModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eeToTwoGammaModel::Initialise(const G4ParticleDefinition*,
                                     const G4DataVector&)
{
  if(isInitialised) { return; }
  fParticleChange = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4eeToTwoGammaModel::ComputeCrossSectionPerElectron(
                                       const G4ParticleDefinition*,
                                       G4double kineticEnergy,
				       G4double, G4double)
{
  // Calculates the cross section per electron of annihilation into two photons
  // from the Heilter formula.

  G4double ekin  = std::max(eV,kineticEnergy);   

  G4double tau   = ekin/electron_mass_c2;
  G4double gam   = tau + 1.0;
  G4double gamma2= gam*gam;
  G4double bg2   = tau * (tau+2.0);
  G4double bg    = sqrt(bg2);

  G4double cross = pi_rcl2*((gamma2+4*gam+1.)*G4Log(gam+bg) - (gam+3.)*bg)
                 / (bg2*(gam+1.));
  return cross;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4eeToTwoGammaModel::ComputeCrossSectionPerAtom(
                                    const G4ParticleDefinition* p,
                                    G4double kineticEnergy, G4double Z,
				    G4double, G4double, G4double)
{
  // Calculates the cross section per atom of annihilation into two photons
  
  G4double cross = Z*ComputeCrossSectionPerElectron(p,kineticEnergy);
  return cross;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4eeToTwoGammaModel::CrossSectionPerVolume(
					const G4Material* material,
					const G4ParticleDefinition* p,
					      G4double kineticEnergy,
					      G4double, G4double)
{
  // Calculates the cross section per volume of annihilation into two photons
  
  G4double eDensity = material->GetElectronDensity();
  G4double cross = eDensity*ComputeCrossSectionPerElectron(p,kineticEnergy);
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Polarisation of gamma according to M.H.L.Pryce and J.C.Ward, 
// Nature 4065 (1947) 435.

void G4eeToTwoGammaModel::SampleSecondaries(vector<G4DynamicParticle*>* vdp,
					    const G4MaterialCutsCouple*,
					    const G4DynamicParticle* dp,
					    G4double,
					    G4double)
{
  G4double PositKinEnergy = dp->GetKineticEnergy();
  G4DynamicParticle *aGamma1, *aGamma2;

  CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();
   
  // Case at rest
  if(PositKinEnergy == 0.0) {
    G4double cost = 2.*rndmEngine->flat()-1.;
    G4double sint = sqrt((1. - cost)*(1. + cost));
    G4double phi  = twopi * rndmEngine->flat();
    G4ThreeVector dir(sint*cos(phi), sint*sin(phi), cost);
    phi = twopi * rndmEngine->flat();
    G4double cosphi = cos(phi);
    G4double sinphi = sin(phi);
    G4ThreeVector pol(cosphi, sinphi, 0.0);
    pol.rotateUz(dir);
    aGamma1 = new G4DynamicParticle(theGamma, dir, electron_mass_c2);
    aGamma1->SetPolarization(pol.x(),pol.y(),pol.z());
    aGamma2 = new G4DynamicParticle(theGamma,-dir, electron_mass_c2);
    pol.set(-sinphi, cosphi, 0.0);
    pol.rotateUz(dir);
    aGamma2->SetPolarization(pol.x(),pol.y(),pol.z());

  } else {

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
      epsil = epsilmin*G4Exp(G4Log(epsilqot)*rndmEngine->flat());
      greject = 1. - epsil + (2.*gam*epsil-1.)/(epsil*tau2*tau2);
      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
    } while( greject < rndmEngine->flat());

    //
    // scattered Gamma angles. ( Z - axis along the parent positron)
    //

    G4double cost = (epsil*tau2-1.)/(epsil*sqg2m1);
    if(std::abs(cost) > 1.0) {
      G4cout << "### G4eeToTwoGammaModel WARNING cost= " << cost
	     << " positron Ekin(MeV)= " << PositKinEnergy
	     << " gamma epsil= " << epsil
	     << G4endl;
      if(cost > 1.0) cost = 1.0;
      else cost = -1.0; 
    }
    G4double sint = sqrt((1.+cost)*(1.-cost));
    G4double phi  = twopi * rndmEngine->flat();

    //
    // kinematic of the created pair
    //

    G4double TotalAvailableEnergy = PositKinEnergy + 2.0*electron_mass_c2;
    G4double Phot1Energy = epsil*TotalAvailableEnergy;

    G4ThreeVector Phot1Direction(sint*cos(phi), sint*sin(phi), cost);
    Phot1Direction.rotateUz(PositDirection);
    aGamma1 = new G4DynamicParticle (theGamma,Phot1Direction, Phot1Energy);
    phi = twopi * rndmEngine->flat();
    G4double cosphi = cos(phi);
    G4double sinphi = sin(phi);
    G4ThreeVector pol(cosphi, sinphi, 0.0);
    pol.rotateUz(Phot1Direction);
    aGamma1->SetPolarization(pol.x(),pol.y(),pol.z());

    G4double Phot2Energy =(1.-epsil)*TotalAvailableEnergy;
    G4double PositP= sqrt(PositKinEnergy*(PositKinEnergy+2.*electron_mass_c2));
    G4ThreeVector dir = PositDirection*PositP - Phot1Direction*Phot1Energy;
    G4ThreeVector Phot2Direction = dir.unit();

    // create G4DynamicParticle object for the particle2
    aGamma2 = new G4DynamicParticle (theGamma,Phot2Direction, Phot2Energy);

    //!!! likely problematic direction to be checked
    pol.set(-sinphi, cosphi, 0.0);
    pol.rotateUz(Phot1Direction);
    cost = pol*Phot2Direction;
    pol -= cost*Phot2Direction;
    pol = pol.unit();
    aGamma2->SetPolarization(pol.x(),pol.y(),pol.z());
  }
  /*
    G4cout << "Annihilation in fly: e0= " << PositKinEnergy
    << " m= " << electron_mass_c2
    << " e1= " << Phot1Energy 
    << " e2= " << Phot2Energy << " dir= " <<  dir 
    << " -> " << Phot1Direction << " " 
    << Phot2Direction << G4endl;
  */
 
  vdp->push_back(aGamma1);
  vdp->push_back(aGamma2);
  fParticleChange->SetProposedKineticEnergy(0.);
  fParticleChange->ProposeTrackStatus(fStopAndKill);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
