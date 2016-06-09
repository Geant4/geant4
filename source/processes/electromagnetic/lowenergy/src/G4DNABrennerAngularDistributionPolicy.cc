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
// $Id: G4DNABrennerAngularDistributionPolicy.cc,v 1.3 2006/06/29 19:39:00 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//

#include "G4DNAElectronElasticBrenner.hh"
#include "G4DNABrennerAngularDistributionPolicy.hh"

#include "Randomize.hh"

G4double                                 G4DNABrennerAngularDistributionPolicy :: RandomizeCosTheta(G4double k, G4int /* z */) const
{
 //  d sigma_el                         1                                 beta(K)
 // ------------ (K) ~ --------------------------------- + ---------------------------------
 //   d Omega           (1 + 2 gamma(K) - std::cos(theta))^2     (1 + 2 delta(K) + std::cos(theta))^2
 //
 // Maximum is < 1/(4 gamma(K)^2) + beta(K)/(4 delta(K)^2)
 //
 // Phys. Med. Biol. 29 N.4 (1983) 443-447
 
 const G4double betaCoeff[5]         = { 7.51525,   -0.419122,    7.2017E-3, -4.646E-5,    1.02897E-7};  
 const G4double deltaCoeff[5]        = { 2.9612,    -0.26376,     4.307E-3,  -2.6895E-5,   5.83505E-8};
 const G4double gamma035_10Coeff[6]  = {-1.7013,    -1.48284,     0.6331,    -0.10911,     8.358E-3,  -2.388E-4};
 const G4double gamma10_100Coeff[5]  = {-3.32517,    0.10996,    -4.5255E-3,  5.8372E-5,  -2.4659E-7};
 const G4double gamma100_200Coeff[3] = { 2.4775E-2, -2.96264E-5, -1.20655E-7};

 // gamma(K), beta(K) and delta(K) are polynomials with coefficients for energy measured in eV
 k/=eV;
 
 G4double beta;
 beta=std::exp(CalculatePolynomial(k, betaCoeff, 5)); 

 G4double delta;
 delta=std::exp(CalculatePolynomial(k, deltaCoeff, 5)); 

 G4double gamma;
 
 if (k>100)
  gamma=CalculatePolynomial(k, gamma100_200Coeff, 3); // Only in this case it is not the exponent of the polynomial
 else if (k>10)
  gamma=std::exp(CalculatePolynomial(k, gamma10_100Coeff, 5));
 else
  gamma=std::exp(CalculatePolynomial(k, gamma035_10Coeff, 6));
  
 G4double oneOverMax;
 oneOverMax=1./(1./(4.*gamma*gamma) + beta/(4.*delta*delta));
 
 G4double cosTheta;
 G4double leftDenominator;
 G4double rightDenominator;
 G4double fCosTheta;
 
 do
 {
  cosTheta = 2.*G4UniformRand()-1.;
  leftDenominator = (1 + 2.*gamma - cosTheta);
  rightDenominator = (1 + 2*delta + cosTheta);
  fCosTheta = oneOverMax*(1./(leftDenominator*leftDenominator) + beta/(rightDenominator*rightDenominator));
 }
 while (fCosTheta < G4UniformRand());
 
 return cosTheta; 
}

G4double                                 G4DNABrennerAngularDistributionPolicy :: CalculatePolynomial(G4double k, const G4double *vector, G4int size) const
{
 // Sum_{i=0}^{size-1} vector_i k^i
 //
 // Phys. Med. Biol. 29 N.4 (1983) 443-447

 G4double result(0);

 while (size>0)
 {
  size--;
  
  result*=k;
  result+=vector[size];
 }
 
 return result;
}

const G4double                           G4DNABrennerAngularDistributionPolicy :: stopAndKillLowEnergyLimit(G4DNABrennerEnergyLimitsPolicy::lowEnergyLimit);
