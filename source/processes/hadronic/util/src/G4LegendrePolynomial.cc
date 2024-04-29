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

#include "G4ios.hh"
#include "G4LegendrePolynomial.hh"
#include "G4Pow.hh"
#include "G4Exp.hh"
#include "G4Log.hh"

using namespace std;

G4double G4LegendrePolynomial::GetCoefficient(size_t i, size_t order)
{
  if(order >= fCoefficients.size()) BuildUpToOrder(order);
  if(order >= fCoefficients.size() ||
     i/2 >= fCoefficients[order].size() ||
     (i%2) != order %2) return 0;
  return fCoefficients[order][i/2];
}

G4double G4LegendrePolynomial::EvalLegendrePoly(G4int order, G4double x)
{
  // Call EvalAssocLegendrePoly with m=0
  return (EvalAssocLegendrePoly(order,0,x));
}

G4double G4LegendrePolynomial::EvalAssocLegendrePoly(G4int l, G4int m, G4double x)
{
  // Invalid calls --> 0. (Keeping for backward compatibility; should we throw instead?)
  if (l<0 || m<-l || m>l) return 0.0;

  // New check: g4pow doesn't handle integer arguments over 512.
  // For us that means:
  if ((l+m) > 512 || (l-m) > 512 || (2*m) > 512) return 0.0; 

  G4Pow* g4pow = G4Pow::GetInstance();
  G4double x2 = x*x;

  // hard-code the first few orders for speed
  switch (l) {
    case 0 : 
      return 1;
    case 1 : 
      switch (m) {
        case -1 : return 0.5 * sqrt(1.-x2); 
        case 0 : return x;
        case 1 : return -sqrt(1.-x2);
      }; 
      break;
    case 2 :
      switch (m) {
        case -2 : return 0.125 * (1.0 - x2);
        case -1 : return 0.5 * x * sqrt(1.0 - x2);
        case 0 : return 0.5*(3.*x2 - 1.);
        case 1 : return -3.*x*sqrt(1.-x2); 
        case 2 : return 3.*(1.-x2);
      };
      break;
    case 3 : 
      switch (m) {
        case -3 : return (1.0/48.0) * (1.0 - x2) * sqrt(1.0 - x2);
        case -2 : return 0.125 * x * (1.0 - x2);
        case -1 : return 0.125 * (5.0 * x2 - 1.0) * sqrt(1.0 - x2);
        case 0 : return 0.5*(5.*x*x2 - 3.*x);
        case 1 : return -1.5*(5.*x2-1.)*sqrt(1.-x2);
        case 2 : return 15.*x*(1.-x2);
        case 3 : return -15.*(1.-x2)*sqrt(1.-x2);
      };
      break;
    case 4 :
      switch (m) {
        case -4 : return (105.0/40320.0)*(1. - 2.*x2 + x2*x2);
        case -3 : return (105.0/5040.0)*x*(1.-x2)*sqrt(1.-x2);
        case -2 : return (15.0/720.0)*(7.*x2-1.)*(1.-x2);
        case -1 : return 0.125*(7.*x*x2-3.*x)*sqrt(1.-x2);
        case 0 : return 0.125*(35.*x2*x2 - 30.*x2 + 3.);
        case 1 : return -2.5*(7.*x*x2-3.*x)*sqrt(1.-x2);
        case 2 : return 7.5*(7.*x2-1.)*(1.-x2);
        case 3 : return -105.*x*(1.-x2)*sqrt(1.-x2);
        case 4 : return 105.*(1. - 2.*x2 + x2*x2);
      };
      break;
  };

  // if m<0, compute P[l,-m,x] and correct
  if (m < 0)
  {
    G4double complementary_value = EvalAssocLegendrePoly(l, -m, x);
    return complementary_value * (m%2==0 ? 1.0 : -1.0) * g4pow->factorial(l+m)/g4pow->factorial(l-m);
  } 

  // Iteratively walk up from P[m,m,x] to P[l,m,x]

  // prime the pump: P[l<m,m,x] = 0
  G4double previous = 0.0;

  // prime the pump: P[m,m,x]
  G4double current; 
  if (m == 0) current = 1.0;
  else if (m == 1) current = -sqrt((1.0 - (x2)));
  else {
    current = (m%2==0 ? 1.0 : -1.0) * 
	     G4Exp(g4pow->logfactorial(2*m) - g4pow->logfactorial(m)) * 
	     G4Exp(G4Log((1.0-(x2))*0.25)*0.5*G4double(m));
  }
  
  // Work up to P[l,m,x] 
  for(G4int i=m+1; i<=l; i++)
  {
    G4double next = (-(G4double(i+m-1))*previous + x*G4double(2*i-1)*current )/G4double(i-m);
    //G4double next = (-(G4double(i+m-1))*previous - x*G4double(2*i-1)*current )/G4double(i-m);  // INCORRECT
    previous = current;
    current = next;
  }
  
  return current;
}

void G4LegendrePolynomial::BuildUpToOrder(size_t orderMax)
{
  if(orderMax > 30) {
    G4cout << "G4LegendrePolynomial::GetCoefficient(): "
           << "I refuse to make a Legendre Polynomial of order " 
           << orderMax << G4endl;
    return;
  }
  while(fCoefficients.size() < orderMax+1) {  /* Loop checking, 30-Oct-2015, G.Folger */
    size_t order = fCoefficients.size();
    fCoefficients.resize(order+1);
    if(order <= 1) fCoefficients[order].push_back(1.);
    else {
      for(size_t iCoeff = 0; iCoeff < order+1; ++iCoeff) {
        if((order % 2) == (iCoeff % 2)) {
          G4double coeff = 0;
          if(iCoeff <= order-2) coeff -= fCoefficients[order-2][iCoeff/2]*G4double(order-1);
          if(iCoeff > 0) coeff += fCoefficients[order-1][(iCoeff-1)/2]*G4double(2*order-1);
          coeff /= G4double(order);
          fCoefficients[order].push_back(coeff);
        }
      }
    }
  }
}

