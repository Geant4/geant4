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

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4BetaDecayType.hh"
#include "G4BetaDecayCorrections.hh"

G4BetaDecayCorrections::G4BetaDecayCorrections(const G4int theZ, const G4int theA)
 : Z(theZ), A(theA)
{
  // alphaZ = fine_structure_const*std::abs(Z);
  alphaZ = fine_structure_const*Z;

  // Nuclear radius in units of hbar/m_e/c
  Rnuc = 0.5*fine_structure_const*std::pow(A, 0.33333);

  // Electron screening potential in units of electrom mass
  V0 = 1.13*fine_structure_const*fine_structure_const
           *std::pow(std::abs(Z), 1.33333);

  gamma0 = std::sqrt(1. - alphaZ*alphaZ);

  // Coefficients for gamma function with real argument
  gc[0] = -0.1010678;
  gc[1] =  0.4245549;
  gc[2] = -0.6998588;
  gc[3] =  0.9512363;
  gc[4] = -0.5748646;
  gc[5] = 1.0;
}


G4double G4BetaDecayCorrections::FermiFunction(const G4double& W)
{
  // Calculate the relativistic Fermi function.  Argument W is the
  // total electron energy in units of electron mass.

  G4double Wprime;
  if (Z < 0) {
    Wprime = W + V0;
  } else {
    Wprime = W - V0;
//    if (Wprime < 1.) Wprime = W;
    if (Wprime <= 1.00001) Wprime = 1.00001;
  }

  G4double p_e = std::sqrt(Wprime*Wprime - 1.);
  G4double eta = alphaZ*Wprime/p_e;
  G4double epieta = std::exp(pi*eta);
  G4double realGamma = Gamma(2.*gamma0+1);
  G4double mod2Gamma = ModSquared(gamma0, eta);

  // Fermi function
  G4double factor1 = 2*(1+gamma0)*mod2Gamma/realGamma/realGamma;
  G4double factor2 = epieta*std::pow(2*p_e*Rnuc, 2*(gamma0-1) );

  // Electron screening factor
  G4double factor3 = (Wprime/W)*std::sqrt( (Wprime*Wprime - 1.)/(W*W - 1.) );

  return factor1*factor2*factor3;
}


G4double
G4BetaDecayCorrections::ModSquared(const G4double& re, const G4double& im)
{
  // Calculate the squared modulus of the Gamma function 
  // with complex argument (re, im) using approximation B 
  // of Wilkinson, Nucl. Instr. & Meth. 82, 122 (1970).
  // Here, choose N = 1 in Wilkinson's notation for approximation B

  G4double factor1 = std::pow( (1+re)*(1+re) + im*im, re+0.5);
  G4double factor2 = std::exp(2*im * std::atan(im/(1+re)));
  G4double factor3 = std::exp(2*(1+re));
  G4double factor4 = 2.*pi;
  G4double factor5 = std::exp( (1+re)/( (1+re)*(1+re) + im*im)/6 );
  G4double factor6 = re*re + im*im;
  return factor1*factor4*factor5/factor2/factor3/factor6;
}


G4double G4BetaDecayCorrections::Gamma(const G4double& arg)
{
  // Use recursion relation to get argument < 1
  G4double fac = 1.0;
  G4double x = arg - 1.;

  G4int loop = 0;
  G4ExceptionDescription ed;
  ed << " While count exceeded " << G4endl; 
  while (x > 1.0) { /* Loop checking, 01.09.2015, D.Wright */
    fac *= x;
    x -= 1.0;
    loop++;
    if (loop > 1000) {
      G4Exception("G4BetaDecayCorrections::Gamma()", "HAD_RDM_100", JustWarning, ed);
      break;
    }
  }

  // Calculation of Gamma function with real argument
  // 0 < arg < 1 using polynomial from Abramowitz and Stegun
  G4double sum = gc[0];
  for (G4int i = 1; i < 6; i++) sum = sum*x + gc[i];

  return sum*fac;
}


G4double
G4BetaDecayCorrections::ShapeFactor(const G4BetaDecayType& bdt,
                                    const G4double& p_e, const G4double& e_nu)
{
  G4double twoPR = 2.*p_e*Rnuc;
  G4double factor(1.);

  switch (bdt)
  {
    case (allowed) :
      break;

    case (firstForbidden) :
      {
        // Parameters for 1st forbidden shape determined from 210Bi data
        // Not valid for other 1st forbidden nuclei
        G4double c1 = 0.578;
        G4double c2 = 28.466;
        G4double c3 = -0.658;

        G4double w = std::sqrt(1. + p_e*p_e);
        factor = 1. + c1*w + c2/w + c3*w*w;
      }
      break;

    case (uniqueFirstForbidden) :
      {
        G4double eta = alphaZ*std::sqrt(1. + p_e*p_e)/p_e;
        G4double gamma1 = std::sqrt(4. - alphaZ*alphaZ);
        G4double gamterm1 = Gamma(2.*gamma0+1.)/Gamma(2.*gamma1+1.);
        G4double term1 = e_nu*e_nu*(1. + gamma0)/6.;
        G4double term2 = 12.*(2. + gamma1)*p_e*p_e
                            *std::pow(twoPR, 2.*(gamma1-gamma0-1) )
                            *gamterm1*gamterm1
                            *ModSquared(gamma1, eta)/ModSquared(gamma0, eta);
        factor = term1 + term2;
      }
      break;

    case (secondForbidden) :
      break;

    case (uniqueSecondForbidden) :
      {
        G4double eta = alphaZ*std::sqrt(1. + p_e*p_e)/p_e;
        G4double gamma1 = std::sqrt(4. - alphaZ*alphaZ);
        G4double gamma2 = std::sqrt(9. - alphaZ*alphaZ);
        G4double gamterm0 = Gamma(2.*gamma0+1.);
        G4double gamterm1 = gamterm0/Gamma(2.*gamma1+1.);
        G4double gamterm2 = gamterm0/Gamma(2.*gamma2+1.);
        G4double term1 = e_nu*e_nu*e_nu*e_nu*(1. + gamma0)/60.;

        G4double term2 = 4.*(2. + gamma1)*e_nu*e_nu*p_e*p_e
                           *std::pow(twoPR, 2.*(gamma1-gamma0-1.) )
                           *gamterm1*gamterm1
                           *ModSquared(gamma1, eta)/ModSquared(gamma0, eta);

        G4double term3 = 180.*(3.+gamma2)*p_e*p_e*p_e*p_e
                             *std::pow(twoPR, 2.*(gamma2-gamma0-2) )
                             *gamterm2*gamterm2
                             *ModSquared(gamma2, eta)/ModSquared(gamma0, eta);

        factor = term1 + term2 + term3;
      }
      break;

    case (thirdForbidden) :
      break;

    case (uniqueThirdForbidden) :
      {
        G4double eta = alphaZ*std::sqrt(1. + p_e*p_e)/p_e;
        G4double gamma1 = std::sqrt(4. - alphaZ*alphaZ);
        G4double gamma2 = std::sqrt(9. - alphaZ*alphaZ);
        G4double gamma3 = std::sqrt(16. - alphaZ*alphaZ);
        G4double gamterm0 = Gamma(2.*gamma0+1.);
        G4double gamterm1 = gamterm0/Gamma(2.*gamma1+1.);
        G4double gamterm2 = gamterm0/Gamma(2.*gamma2+1.);
        G4double gamterm3 = gamterm0/Gamma(2.*gamma3+1.);

        G4double term1 = e_nu*e_nu*e_nu*e_nu*e_nu*e_nu*(1. + gamma0)/1260.;

        G4double term2 = 2.*(2. + gamma1)*e_nu*e_nu*e_nu*e_nu*p_e*p_e
                           *std::pow(twoPR, 2.*(gamma1-gamma0-1.) )
                           *gamterm1*gamterm1
                           *ModSquared(gamma1, eta)/ModSquared(gamma0, eta)/5.;

        G4double term3 = 60.*(3.+gamma2)*p_e*p_e*p_e*p_e*e_nu*e_nu
                             *std::pow(twoPR, 2.*(gamma2-gamma0-2.) )
                             *gamterm2*gamterm2
                             *ModSquared(gamma2, eta)/ModSquared(gamma0, eta);

        G4double term4 = 2240.*p_e*p_e*p_e*p_e*p_e*p_e*(4. + gamma3)
                             *std::pow(twoPR, 2.*(gamma3-gamma0-3.) )
                             *gamterm3*gamterm3
                             *ModSquared(gamma3, eta)/ModSquared(gamma0, eta);

        factor = term1 + term2 + term3 + term4;
      }
      break;

    default:
      G4Exception("G4BetaDecayCorrections::ShapeFactor()","HAD_RDM_010",
                  JustWarning,
                  "Transition not yet implemented - using allowed shape");
      break;
  }
  return factor;
}

 
