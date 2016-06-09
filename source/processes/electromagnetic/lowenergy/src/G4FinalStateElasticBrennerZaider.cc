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
// $Id: G4FinalStateElasticBrennerZaider.cc,v 1.1 2007/10/12 23:11:41 pia Exp $
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


#include "G4FinalStateElasticBrennerZaider.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4DynamicParticle.hh"
#include "Randomize.hh"

#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4Electron.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleMomentum.hh"

G4FinalStateElasticBrennerZaider::G4FinalStateElasticBrennerZaider()
{
  // These data members will be used in the next implementation iteration, 
  // when the enriched PhysicsModel policy is implemented
  name = "FinalStateElasticBrennerZaider";
  lowEnergyLimit = 7.4 * eV;
  highEnergyLimit = 10 * MeV;

  betaCoeff.push_back(7.51525);
  betaCoeff.push_back(-0.41912);    
  betaCoeff.push_back(7.2017E-3);
  betaCoeff.push_back(-4.646E-5);    
  betaCoeff.push_back(1.02897E-7);

  deltaCoeff.push_back(2.9612); 
  deltaCoeff.push_back(-0.26376); 
  deltaCoeff.push_back(4.307E-3); 
  deltaCoeff.push_back(-2.6895E-5);
  deltaCoeff.push_back(5.83505E-8);

  gamma035_10Coeff.push_back(-1.7013); 
  gamma035_10Coeff.push_back(-1.48284); 
  gamma035_10Coeff.push_back(0.6331); 
  gamma035_10Coeff.push_back(-0.10911); 
  gamma035_10Coeff.push_back(8.358E-3); 
  gamma035_10Coeff.push_back(-2.388E-4);

  gamma10_100Coeff.push_back(-3.32517); 
  gamma10_100Coeff.push_back(0.10996); 
  gamma10_100Coeff.push_back(-4.5255E-3); 
  gamma10_100Coeff.push_back(5.8372E-5); 
  gamma10_100Coeff.push_back(-2.4659E-7);

  gamma100_200Coeff.push_back(2.4775E-2);
  gamma100_200Coeff.push_back(-2.96264E-5);
  gamma100_200Coeff.push_back(-1.20655E-7);
}


G4FinalStateElasticBrennerZaider::~G4FinalStateElasticBrennerZaider()
{ 
  // empty
  // G4DynamicParticle objects produced are owned by client
}
 

const G4FinalStateProduct& G4FinalStateElasticBrennerZaider::GenerateFinalState(const G4Track& track, const G4Step& /* step */)
{
  // Clear previous secondaries, energy deposit and particle kill status
  product.Clear();

  // Kinetic energy of primary particle
  G4double k = track.GetDynamicParticle()->GetKineticEnergy();

  // Assume material = water; H2O number of electrons
  // ---- MGP ---- To be generalized later
  // const G4int z = 10; 

  G4double cosTheta = RandomizeCosTheta(k);
  
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

G4double G4FinalStateElasticBrennerZaider::RandomizeCosTheta(G4double k)
{
  //  d sigma_el                         1                                 beta(K)
  // ------------ (K) ~ --------------------------------- + ---------------------------------
  //   d Omega           (1 + 2 gamma(K) - cos(theta))^2     (1 + 2 delta(K) + cos(theta))^2
  //
  // Maximum is < 1/(4 gamma(K)^2) + beta(K)/(4 delta(K)^2)
  //
  // Phys. Med. Biol. 29 N.4 (1983) 443-447
  
  // gamma(K), beta(K) and delta(K) are polynomials with coefficients for energy measured in eV
  k /= eV;
  
  G4double beta = std::exp(CalculatePolynomial(k,betaCoeff)); 
  G4double delta = std::exp(CalculatePolynomial(k,deltaCoeff)); 
  
  G4double gamma;
  if (k > 100.)
    {
      gamma = CalculatePolynomial(k, gamma100_200Coeff); // Only in this case it is not the exponent of the polynomial
    }  
  else 
    {

      if (k>10)
	{
	  gamma = std::exp(CalculatePolynomial(k, gamma10_100Coeff));
	}
      else
	{
	  gamma = std::exp(CalculatePolynomial(k, gamma035_10Coeff));
	}
    }
  
  // G4cout << "beta = " << beta << ", gamma = " << gamma << ", delta = " << delta << G4endl;

  G4double oneOverMax = 1. / (1./(4.*gamma*gamma) + beta/(4.*delta*delta));
  
  G4double cosTheta = 0.;
  G4double leftDenominator = 0.;
  G4double rightDenominator = 0.;
  G4double fCosTheta = 0.;
  
  do
    {
      cosTheta = 2. * G4UniformRand() - 1.;
      leftDenominator = (1 + 2.*gamma - cosTheta);
      rightDenominator = (1 + 2.*delta + cosTheta);
      fCosTheta = oneOverMax * (1./(leftDenominator*leftDenominator) + beta/(rightDenominator*rightDenominator));
    }
  while (fCosTheta < G4UniformRand());

  //  G4cout << "cosTheta = " << cosTheta << G4endl;
  
  return cosTheta; 
}

G4double G4FinalStateElasticBrennerZaider::CalculatePolynomial(G4double k, std::vector<G4double>& vec)
{
  // Sum_{i=0}^{size-1} vector_i k^i
  //
  // Phys. Med. Biol. 29 N.4 (1983) 443-447

  G4double result = 0.;
  size_t size = vec.size();

  while (size>0)
    {
      size--; 
      
      result *= k;
      result += vec[size];
    }
  
  return result;
}

