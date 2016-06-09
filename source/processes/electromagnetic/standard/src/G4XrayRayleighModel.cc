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
//
// Author: Vladimir Grichine
//
// History:
// 
// 25.052011   first implementation

#include "G4XrayRayleighModel.hh"

////////////////////////////////////////////////////////////////////////////////////

using namespace std;

const G4double G4XrayRayleighModel::fCofA = 2.*pi2*Bohr_radius*Bohr_radius;

const G4double G4XrayRayleighModel::fCofR = 8.*pi*classic_electr_radius*classic_electr_radius/3.;

//////////////////////////////////////////////////////////////////////////////////.

G4XrayRayleighModel::G4XrayRayleighModel(const G4ParticleDefinition*,
						   const G4String& nam)
  :G4VEmModel(nam),isInitialised(false)
{
  fParticleChange = 0;
  lowEnergyLimit  = 250*eV; 
  highEnergyLimit = 1.*MeV;
  
  //  SetLowEnergyLimit(lowEnergyLimit);
  SetHighEnergyLimit(highEnergyLimit);
  //
  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if(verboseLevel > 0) 
  {
    G4cout << "Xray Rayleigh is constructed " << G4endl
	   << "Energy range: "
	   << lowEnergyLimit / eV << " eV - "
	   << highEnergyLimit / MeV << " MeV"
	   << G4endl;
  }
}

//////////////////////////////////////////////////////////////////////////////////

G4XrayRayleighModel::~G4XrayRayleighModel()
{  

}

/////////////////////////////////////////////////////////////////////////////////////

void G4XrayRayleighModel::Initialise(const G4ParticleDefinition* particle,
					  const G4DataVector& cuts)
{
  if (verboseLevel > 3) 
  {
    G4cout << "Calling G4XrayRayleighModel::Initialise()" << G4endl;
  }

  InitialiseElementSelectors(particle,cuts);


  if(isInitialised) return; 
  fParticleChange = GetParticleChangeForGamma();
  isInitialised = true;

}

/////////////////////////////////////////////////////////////////////////////////

G4double G4XrayRayleighModel::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                             G4double gammaEnergy,
                                             G4double Z, G4double,
                                             G4double, G4double)
{
  if (verboseLevel > 3) 
  {
    G4cout << "Calling CrossSectionPerAtom() of G4XrayRayleighModel" << G4endl;
  }
  if (gammaEnergy < lowEnergyLimit || gammaEnergy > highEnergyLimit) 
  {
    return 0.0;
  }
  const G4double p0 =  1.78076e+00;  // 1.77457e+00;
  const G4double p1 = -6.0911e-02;    // -1.78171e-02;
  // const G4double p2 =  4.60444e-04;
  G4double    alpha = p0 + p1*std::log(Z); //  + p2*Z*Z;

  G4double k    = gammaEnergy/hbarc;
  G4double k2   = std::pow(k, alpha); // k*k; 1.69 for Z=6, 1.5 for Z=82
  G4double z2p3 = std::pow( Z, 0.66667); // 2/3~0.66667
  G4double a    = fCofA*z2p3*k2; 
  G4double b    = 1. + 2.*a;
  G4double b2   = b*b;
  G4double b3   = b*b2;
  G4double xsc  = fCofR*Z*Z/b3;
           xsc *= a*a + (1. + a)*(1. + a);  
  return   xsc;   

}

///////////////////////////////////////////////////////////////////////////////////

void G4XrayRayleighModel::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
					    const G4MaterialCutsCouple*, // couple,
					      const G4DynamicParticle* aDynamicGamma,
					      G4double,
					      G4double)
{
  if ( verboseLevel > 3)
  {
    G4cout << "Calling SampleSecondaries() of G4XrayRayleighModel" << G4endl;
  }
  G4double photonEnergy0 = aDynamicGamma->GetKineticEnergy();

  G4ParticleMomentum photonDirection0 = aDynamicGamma->GetMomentumDirection();


  // Sample the angle of the scattered photon
  // according to 1 + cosTheta*cosTheta distribution

  G4double cosTheta, sinTheta;
  G4double c, delta, cofA, signc = 1., a, power = 1./3.;

  c = 4. - 8.*G4UniformRand();
  a = c;
 
  if( c < 0. )
  {
    signc = -1.;
    a     = -c;
  }
  delta  = std::sqrt(a*a+4.);
  delta += a;
  delta *= 0.5; 
  cofA = -signc*std::pow(delta, power);
  cosTheta = cofA - 1./cofA;

  if( cosTheta >  1.) cosTheta =  1.;
  if( cosTheta < -1.) cosTheta = -1.;

  sinTheta = std::sqrt( (1. - cosTheta)*(1. + cosTheta) );

  // Scattered photon angles. ( Z - axis along the parent photon)

  G4double phi = twopi * G4UniformRand() ;
  G4double dirX = sinTheta*std::cos(phi);
  G4double dirY = sinTheta*std::sin(phi);
  G4double dirZ = cosTheta;

  // Update G4VParticleChange for the scattered photon

  G4ThreeVector photonDirection1(dirX, dirY, dirZ);
  photonDirection1.rotateUz(photonDirection0);
  fParticleChange->ProposeMomentumDirection(photonDirection1);

  fParticleChange->SetProposedKineticEnergy(photonEnergy0); 
}


