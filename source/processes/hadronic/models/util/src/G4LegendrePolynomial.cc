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

  if(l<0 || m>abs(l)) 
      return 0;

  // hard-code the first few orders for speed
    if(x==0) return 1;

    if(l==1) 
        m==0 ? return x : 
               return -sqrt(1.-x*x); /*m==1*/ 
          
    if(l<5) {
      G4double x2 = x*x;
      switch (l) {
          case 2: 
              switch (m):
                  case 0: 
                      return 0.5*(3.*x2 - 1.);
                  case 1:
                      return -3.*x*sqrt(1.-x2);

                  default:
                      return 3.*(1.-x2); /*m==2*/
    
          case 3:
              switch (m) {
                  case 0:
                      return 0.5*(5.*x*x2 - 3.*x);
                  case 1:
                      return -1.5*(5.*x2-1.)*sqrt(1.-x2);
                  case 2: 
                      return 15.*x*(1.-x2);

                  default:
                      return -15.*(1.-x2)*sqrt(1.-x2); /*m==3*/ 
    
          case 4:
              switch (m) {
                  case 0: 
                      return 0.125*(35.*x2*x2 - 30.*x2 + 3.);
                  case 1:
                      return -2.5*(7.*x*x2-3.*x)*sqrt(1.-x2);
                  case 2: 
                      return 7.5*(7.*x2-1.)*(1.-x2);
                  case 3: 
                      return -105.*x*(1.-x2)*sqrt(1.-x2);

                  default:
                      return 105.*(1. - 2.*x2 + x2*x2); /*m==4*/
    }
  }

  // See if we have this value cached.
  if(cache != NULL && cache->count(l) > 0 && (*cache)[l].count(m) > 0) {
    return (*cache)[l][m];
  }

  // Otherwise calculate it using std::assoc_legendre(l, m, x)
  G4double value = assoc_legendre(l, m, x);

  // If we are working with a cache, cache this value.
  if(cache) {
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

