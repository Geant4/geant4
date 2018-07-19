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

G4double G4LegendrePolynomial::EvalAssocLegendrePoly(G4int l, G4int m, G4double x, 
                                                     map<G4int, map<G4int, G4double> >* cache)
{
  // Calculate P_l^m(x).
  // If cache ptr is non-null, use cache[l][m] if it exists, otherwise compute
  // P_l^m(x) and cache it in that position. The cache speeds up calculations
  // where many P_l^m computations are need at the same value of x.

  if(l<0 || m<-l || m>l) return 0;
  G4Pow* g4pow =  G4Pow::GetInstance();

  // Use non-log factorial for low l, m: it is more efficient until 
  // l and m get above 10 or so.
  // FIXME: G4Pow doesn't check whether the argument gets too large, 
  // which is unsafe! Max is 512; VI: It is assume that Geant4 does not
  // need higher order
  if(m<0) {
    G4double value = (m%2 ? -1. : 1.) * EvalAssocLegendrePoly(l, -m, x);
    if(l < 10) return value * g4pow->factorial(l+m)/g4pow->factorial(l-m);
    else { return value * G4Exp(g4pow->logfactorial(l+m) - g4pow->logfactorial(l-m));
    }
  }

  // hard-code the first few orders for speed
  if(l==0)   return 1;
  if(l==1) {
    if(m==0) return x;
    /*m==1*/ return -sqrt(1.-x*x);
  }
  if(l<5) {
    G4double x2 = x*x;
    if(l==2) {
      if(m==0) return 0.5*(3.*x2 - 1.);
      if(m==1) return -3.*x*sqrt(1.-x2);
      /*m==2*/ return 3.*(1.-x2);
    }
    if(l==3) {
      if(m==0) return 0.5*(5.*x*x2 - 3.*x);
      if(m==1) return -1.5*(5.*x2-1.)*sqrt(1.-x2);
      if(m==2) return 15.*x*(1.-x2);
      /*m==3*/ return -15.*(1.-x2)*sqrt(1.-x2);
    }
    if(l==4) {
      if(m==0) return 0.125*(35.*x2*x2 - 30.*x2 + 3.);
      if(m==1) return -2.5*(7.*x*x2-3.*x)*sqrt(1.-x2);
      if(m==2) return 7.5*(7.*x2-1.)*(1.-x2);
      if(m==3) return -105.*x*(1.-x2)*sqrt(1.-x2);
      /*m==4*/ return 105.*(1. - 2.*x2 + x2*x2);
    }
  }

  // Easy special cases
  // FIXME: G4Pow doesn't check whether the argument gets too large, which is unsafe! Max is 512. 
  if(m==l) return (l%2 ? -1. : 1.) * 
	     G4Exp(g4pow->logfactorial(2*l) - g4pow->logfactorial(l)) * 
	     G4Exp(G4Log((1.-x*x)*0.25)*0.5*G4double(l));
  if(m==l-1) return x*(2.*G4double(m)+1.)*EvalAssocLegendrePoly(m,m,x);

  // See if we have this value cached.
  if(cache != NULL && cache->count(l) > 0 && (*cache)[l].count(m) > 0) {
    return (*cache)[l][m];
  }

  // Otherwise calculate recursively
  G4double value = (x*G4double(2*l-1)*EvalAssocLegendrePoly(l-1,m,x) - 
	           (G4double(l+m-1))*EvalAssocLegendrePoly(l-2,m,x))/G4double(l-m);

  // If we are working with a cache, cache this value.
  if(cache != NULL) {
    (*cache)[l][m] = value;
  }
  return value;
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

