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
// $Id: G4FinalStateElasticBrennerZaider.cc,v 1.11 2009-06-11 15:47:08 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4FinalStateElasticBrennerZaider.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4FinalStateElasticBrennerZaider::G4FinalStateElasticBrennerZaider()
{
  lowEnergyLimit = 8.23 * eV; // SI : i/o of 7.4 * eV;
  highEnergyLimit = 200 * eV;

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

   G4cout << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "   The class G4FinalStateElasticBrennerZaider is NOT SUPPORTED ANYMORE. " << G4endl;
   G4cout << "   It will be REMOVED with the next major release of Geant4. " << G4endl;
   G4cout << "   Please consult: https://twiki.cern.ch/twiki/bin/view/Geant4/LoweProcesses" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4FinalStateElasticBrennerZaider::~G4FinalStateElasticBrennerZaider()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4FinalStateProduct& G4FinalStateElasticBrennerZaider::GenerateFinalState(const G4Track& track, const G4Step& /* step */)
{
  product.Clear();

  G4double k = track.GetDynamicParticle()->GetKineticEnergy();
  
  if ( k>= lowEnergyLimit && k<=highEnergyLimit )
  {
    G4double cosTheta = RandomizeCosTheta(k);
  
    G4double phi = 2. * pi * G4UniformRand();

    G4ThreeVector zVers = track.GetDynamicParticle()->GetMomentumDirection();
    G4ThreeVector xVers = zVers.orthogonal();
    G4ThreeVector yVers = zVers.cross(xVers);

    G4double xDir = std::sqrt(1. - cosTheta*cosTheta);
    G4double yDir = xDir;
    xDir *= std::cos(phi);
    yDir *= std::sin(phi);

    G4ThreeVector zPrimeVers((xDir*xVers + yDir*yVers + cosTheta*zVers));

    product.ModifyPrimaryParticle(zPrimeVers,k);
  }

  if (k<lowEnergyLimit)
  {
    product.KillPrimaryParticle();
  }
  
  return product;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4FinalStateElasticBrennerZaider::RandomizeCosTheta(G4double k)
{
  //  d sigma_el                         1                                 beta(K)
  // ------------ (K) ~ --------------------------------- + ---------------------------------
  //   d Omega           (1 + 2 gamma(K) - cos(theta))^2     (1 + 2 delta(K) + cos(theta))^2
  //
  // Maximum is < 1/(4 gamma(K)^2) + beta(K)/((2+2delta(K))^2)
  //
  // Phys. Med. Biol. 29 N.4 (1983) 443-447
  
  // gamma(K), beta(K) and delta(K) are polynomials with coefficients for energy measured in eV

  k /= eV;
  
  G4double beta = std::exp(CalculatePolynomial(k,betaCoeff)); 
  G4double delta = std::exp(CalculatePolynomial(k,deltaCoeff)); 
  G4double gamma;
  
  if (k > 100.)
  {
      gamma = CalculatePolynomial(k, gamma100_200Coeff); 
      // Only in this case it is not the exponent of the polynomial
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

  G4double oneOverMax = 1. / (1./(4.*gamma*gamma) + beta/( (2.+2.*delta)*(2.+2.*delta) ));
  
  G4double cosTheta = 0.;
  G4double leftDenominator = 0.;
  G4double rightDenominator = 0.;
  G4double fCosTheta = 0.;
  
  do
  {
      cosTheta = 2. * G4UniformRand() - 1.;
      leftDenominator = (1. + 2.*gamma - cosTheta);
      rightDenominator = (1. + 2.*delta + cosTheta);
      if ( (leftDenominator * rightDenominator) != 0. )
      {
	  fCosTheta = oneOverMax * (1./(leftDenominator*leftDenominator) + beta/(rightDenominator*rightDenominator));
      }
  }
  while (fCosTheta < G4UniformRand());

  return cosTheta; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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

