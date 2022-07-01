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
// 14.10.12 V.Grichine, update of xsc and angular distribution
// 25.05.2011   first implementation

#include "G4XrayRayleighModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//////////////////////////////////////////////////////////////////////////////////

const G4double G4XrayRayleighModel::fCofA = 2.*pi2*Bohr_radius*Bohr_radius;

const G4double G4XrayRayleighModel::fCofR = 8.*pi*classic_electr_radius*classic_electr_radius/3.;

//////////////////////////////////////////////////////////////////////////////////

G4XrayRayleighModel::G4XrayRayleighModel(const G4ParticleDefinition*,
						   const G4String& nam)
  :G4VEmModel(nam),isInitialised(false)
{
  fParticleChange = nullptr;
  lowEnergyLimit  = 250*eV; 
  highEnergyLimit = 10.*MeV;
  fFormFactor     = 0.0;
  
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

G4XrayRayleighModel::~G4XrayRayleighModel() = default;

//////////////////////////////////////////////////////////////////////////////////

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

//////////////////////////////////////////////////////////////////////////////////

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
  G4double k   = gammaEnergy/hbarc;
           k  *= Bohr_radius;
  G4double p0  =  0.680654;  
  G4double p1  = -0.0224188;
  G4double lnZ = std::log(Z);    

  G4double lna = p0 + p1*lnZ; 

  G4double  alpha = std::exp(lna);

  G4double fo   = std::pow(k, alpha); 

  p0 = 3.68455;
  p1 = -0.464806;
  lna = p0 + p1*lnZ; 

  fo *= 0.01*std::exp(lna);

  fFormFactor = fo;

  G4double b    = 1. + 2.*fo;
  G4double b2   = b*b;
  G4double b3   = b*b2;

  G4double xsc  = fCofR*Z*Z/b3;
           xsc *= fo*fo + (1. + fo)*(1. + fo);  


  return   xsc;   

}

//////////////////////////////////////////////////////////////////////////////////

void G4XrayRayleighModel::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,  
                                            const G4MaterialCutsCouple* couple,
                                            const G4DynamicParticle* aDPGamma,
                                            G4double,
                                            G4double)
{
  if ( verboseLevel > 3)
  {
    G4cout << "Calling SampleSecondaries() of G4XrayRayleighModel" << G4endl;
  }
  G4double photonEnergy0 = aDPGamma->GetKineticEnergy();

  G4ParticleMomentum photonDirection0 = aDPGamma->GetMomentumDirection();


  // Sample the angle of the scattered photon
  // according to 1 + cosTheta*cosTheta distribution

  G4double cosDipole, cosTheta, sinTheta;
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
  cosDipole = cofA - 1./cofA;

  // select atom
  const G4Element* elm = SelectTargetAtom(couple, aDPGamma->GetParticleDefinition(),
                                          photonEnergy0,aDPGamma->GetLogKineticEnergy());
  G4double Z = elm->GetZ();

  G4double k   = photonEnergy0/hbarc;
           k  *= Bohr_radius;
  G4double p0  =  0.680654;  
  G4double p1  = -0.0224188;
  G4double lnZ = std::log(Z);    

  G4double lna = p0 + p1*lnZ; 

  G4double  alpha = std::exp(lna);

  G4double fo   = std::pow(k, alpha); 

  p0 = 3.68455;
  p1 = -0.464806;
  lna = p0 + p1*lnZ; 

  fo *= 0.01*pi*std::exp(lna);

  
  G4double beta = fo/(1 + fo);

  cosTheta = (cosDipole + beta)/(1. + cosDipole*beta);


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


