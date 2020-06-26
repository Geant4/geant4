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
// ------------------------------------------------------------
// GEANT 4 class implementation
// ------------------------------------------------------------

#include <cmath>
#include <string.h>
#include "Gamma.hh"

MyGamma::MyGamma(){}

MyGamma::~MyGamma(){}

//____________________________________________________________________________
double MyGamma::Gamma(double z)
{
  if (z <= 0)
      return 0;

  return std::tgamma(z);
}

//____________________________________________________________________________
double MyGamma::Gamma(double a,double x)
{
  // Computation of the incomplete gamma function P(a,x)
  //
  // The algorithm is based on the formulas and code as denoted in
  // Numerical Recipes 2nd ed. on p. 210-212 (W.H.Press et al.).
  //
  //--- Nve 14-nov-1998 UU-SAP Utrecht
  
  if (a <= 0 || x <= 0) return 0;
  
  if (x < (a+1)) return GamSer(a,x);
  else           return GamCf(a,x);
}

//____________________________________________________________________________
double MyGamma::GamCf(double a,double x)
{
  // Computation of the incomplete gamma function P(a,x)
  // via its continued fraction representation.
  //
  // The algorithm is based on the formulas and code as denoted in
  // Numerical Recipes 2nd ed. on p. 210-212 (W.H.Press et al.).
  //
  //--- Nve 14-nov-1998 UU-SAP Utrecht
  
  int itmax    = 100;      // Maximum number of iterations
  double eps   = 3.e-7;    // Relative accuracy
  double fpmin = 1.e-30;   // Smallest double value allowed here
  
  if (a <= 0 || x <= 0) return 0;
  
  double gln = LnGamma(a);
  double b   = x+1-a;
  double c   = 1/fpmin;
  double d   = 1/b;
  double h   = d;
  double an,del;
  for (int i=1; i<=itmax; i++) {
    an = double(-i)*(double(i)-a);
    b += 2;
    d  = an*d+b;
    if (Abs(d) < fpmin) d = fpmin;
    c = b+an/c;
    if (Abs(c) < fpmin) c = fpmin;
    d   = 1/d;
    del = d*c;
    h   = h*del;
    if (Abs(del-1) < eps) break;
    //if (i==itmax) cout << "*GamCf(a,x)* a too large or itmax too small" << endl;
  }
  double v = Exp(-x+a*Log(x)-gln)*h;
  return (1-v);
}

//____________________________________________________________________________
double MyGamma::GamSer(double a,double x)
{
  // Computation of the incomplete gamma function P(a,x)
  // via its series representation.
  //
  // The algorithm is based on the formulas and code as denoted in
  // Numerical Recipes 2nd ed. on p. 210-212 (W.H.Press et al.).
  //
  //--- Nve 14-nov-1998 UU-SAP Utrecht
  
  int itmax  = 100;   // Maximum number of iterations
  double eps = 3.e-7; // Relative accuracy
  
  if (a <= 0 || x <= 0) return 0;
  
  double gln = LnGamma(a);
  double ap  = a;
  double sum = 1/a;
  double del = sum;
  for (int n=1; n<=itmax; n++) {
    ap  += 1;
    del  = del*x/ap;
    sum += del;
    if (MyGamma::Abs(del) < Abs(sum*eps)) break;
    //if (n==itmax) cout << "*GamSer(a,x)* a too large or itmax too small" << endl;
  }
  double v = sum*Exp(-x+a*Log(x)-gln);
  return v;
}


double MyGamma::LnGamma(double z)
{
  if (z <= 0)
      return 0;

  return std::lgamma(z);
}
