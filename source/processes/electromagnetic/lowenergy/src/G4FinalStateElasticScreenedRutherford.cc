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
// $Id: G4FinalStateElasticScreenedRutherford.cc,v 1.2 2007/10/12 23:10:33 pia Exp $
// GEANT4 tag $Name: geant4-09-01 $
// 
// Contact Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// Reference: TNS Geant4-DNA paper
// Reference for implementation model: NIM. 155, pp. 145-156, 1978

// History:
// -----------
// Date         Name              Modification
// 28 Apr 2007  M.G. Pia          Created in compliance with design described in TNS paper
//
// -------------------------------------------------------------------

// Class description:
// Reference: TNS Geant4-DNA paper
// S. Chauvie et al., Geant4 physics processes for microdosimetry simulation:
// design foundation and implementation of the first set of models,
// IEEE Trans. Nucl. Sci., vol. 54, no. 6, Dec. 2007.
// Further documentation available from http://www.ge.infn.it/geant4/dna

// -------------------------------------------------------------------


#include "G4FinalStateElasticScreenedRutherford.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4DynamicParticle.hh"
#include "Randomize.hh"

#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleMomentum.hh"

G4FinalStateElasticScreenedRutherford::G4FinalStateElasticScreenedRutherford()
{
  // These data members will be used in the next implementation iteration, 
  // when the enriched PhysicsModel policy is implemented
  name = "FinalStateElasticScreenedRutherford";
  lowEnergyLimit = 7.4 * eV;
  highEnergyLimit = 10 * MeV;
}


G4FinalStateElasticScreenedRutherford::~G4FinalStateElasticScreenedRutherford()
{ 
  // empty
  // G4DynamicParticle objects produced are owned by client
}
 

const G4FinalStateProduct& G4FinalStateElasticScreenedRutherford::GenerateFinalState(const G4Track& track, const G4Step& step)
{
  // Clear previous secondaries, energy deposit and particle kill status
  product.Clear();

  // Kinetic energy of primary particle
  G4double k = track.GetDynamicParticle()->GetKineticEnergy();

  // Assume material = water; H2O number of electrons
  // ---- MGP ---- To be generalized later
  const G4int z = 10; 

  G4double cosTheta = RandomizeCosTheta(k, z);
  
  G4double phi = 2. * pi * G4UniformRand();

  // G4cout << "cosTheta in GenerateFinalState = " << cosTheta << ", phi = " << phi << G4endl;

  G4ThreeVector zVers = track.GetDynamicParticle()->GetMomentumDirection();
  G4ThreeVector xVers = zVers.orthogonal();
  G4ThreeVector yVers = zVers.cross(xVers);

  G4double xDir = std::sqrt(1. - cosTheta*cosTheta);
  G4double yDir = xDir;
  xDir *= std::cos(phi);
  yDir *= std::sin(phi);

  // G4cout << "xDir, yDir = " << xDir <<", " << yDir << G4endl;

  // G4ThreeVector zPrimeVers((xDir*xVers + yDir*yVers + cosTheta*zVers).unit());
  G4ThreeVector zPrimeVers((xDir*xVers + yDir*yVers + cosTheta*zVers));

  // G4cout << "zPrimeVers = (" << zPrimeVers.x() << ", "<< zPrimeVers.y() << ", "<< zPrimeVers.z() << ") " << G4endl;

  //  product.ModifyPrimaryParticle(zPrimeVers.x(),zPrimeVers.y(),zPrimeVers.z(),k);
  product.ModifyPrimaryParticle(zPrimeVers,k);

  //  this->aParticleChange.ProposeEnergy(k);
  //  this->aParticleChange.ProposeMomentumDirection(zPrimeVers);
  //  this->aParticleChange.SetNumberOfSecondaries(0);

  return product;
}

G4double G4FinalStateElasticScreenedRutherford::RandomizeCosTheta(G4double k, G4int z) const
{

 //  d sigma_el                sigma_Ruth(K)
 // ------------ (K) ~ -----------------------------
 //   d Omega           (1 + 2 n(K) - cos(theta))^2
 //
 // We extract cos(theta) distributed as (1 + 2 n(K) - cos(theta))^-2
 //
 // Maximum is for theta=0: 1/(4 n(K)^2) (When n(K) is positive, that is always satisfied within the validity of the process)
 //
 // Phys. Med. Biol. 45 (2000) 3171-3194

 G4double n = ScreeningFactor(k, z);

 G4double oneOverMax = (4. * n*n);

 G4double cosTheta;
 G4double fCosTheta;

 do 
   { 
     cosTheta = 2. * G4UniformRand() - 1.;
     fCosTheta = (1 + 2.*n - cosTheta);
     fCosTheta = oneOverMax / (fCosTheta*fCosTheta);
   }
 while (fCosTheta < G4UniformRand());
 
 return cosTheta;
}

G4double G4FinalStateElasticScreenedRutherford::ScreeningFactor(G4double k, G4int z) const
{
  //
  //         alpha_1 + beta_1 ln(K/eV)   constK Z^(2/3)
  // n(T) = -------------------------- -----------------
  //              K/(m_e c^2)            2 + K/(m_e c^2)
  //
  // Where K is the electron non-relativistic kinetic energy
  //
  // n(T) > 0 for T < ~ 400 MeV
  //
  // Nucl. Instr. Meth. 155 (1978) 145-156
  
  const G4double alpha_1 = 1.64;
  const G4double beta_1 = -0.0825;
  const G4double constK = 1.7E-5;
  
  G4double numerator = (alpha_1 + beta_1 * std::log(k/eV)) * constK * std::pow(static_cast<double>(z), 2./3.);
  
  k /= electron_mass_c2;
  
  G4double denominator;
  denominator = k * (2 + k);
  
  G4double result = 0.;
  if (denominator != 0.) 
    {
      result = numerator / denominator;
    }
  else
    {
      // Throw an exception
      G4Exception("G4FinalStateElasticScreenedRutherford::ScreeningFactor - denominator = 0");
    }
  return result;

}
