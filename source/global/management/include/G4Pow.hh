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
// $Id: G4Pow.hh 69386 2013-05-02 10:35:42Z vnivanch $
//
//
// -------------------------------------------------------------------
//
// Class G4Pow
//
// Class description:
//
// Utility singleton class for the fast computation of log and pow
// functions. Integer argument should in the interval 0-512, no
// check is performed inside these methods for performance reasons.
// For factorial integer argument should be in the interval 0-170
// Computations with double arguments are fast for the interval
// 0.5-255.5, standard library is used in the opposite case

// Author: Vladimir Ivanchenko
//
// Creation date: 23.05.2009
// -------------------------------------------------------------------

#ifndef G4Pow_h
#define G4Pow_h 1

#include "globals.hh"
#include "G4DataVector.hh"

class G4Pow
{

  public:

    static G4Pow* GetInstance();

    // Fast computation of Z^1/3
    //
    inline G4double Z13(G4int Z);
    inline G4double A13(G4double A);

    // Fast computation of Z^2/3
    //
    inline G4double Z23(G4int Z);
    inline G4double A23(G4double A);

    // Fast computation of log(Z)
    //
    inline G4double logZ(G4int Z);
    inline G4double logA(G4double A);

    // Fast computation of log10(Z)
    //
    inline G4double log10Z(G4int Z);
    inline G4double log10A(G4double A);

    // Fast computation of pow(Z,X)
    //
    inline G4double powZ(G4int Z, G4double y);
    inline G4double powA(G4double A, G4double y);
           G4double powN(G4double x, G4int n);

    // Fast factorial
    //
    inline G4double factorial(G4int Z);
    inline G4double logfactorial(G4int Z);

  private:

    G4Pow();
   ~G4Pow();

  private:

    static G4Pow* fpInstance;

    const G4double onethird;
    const G4double minA;
    const G4double maxA;

    G4DataVector pz13;
    G4DataVector lz;
    G4DataVector fact;
    G4DataVector logfact;
};

// -------------------------------------------------------------------

inline G4double G4Pow::Z13(G4int Z)
{
  return pz13[Z];
}

inline G4double G4Pow::A13(G4double A)
{
  G4double res;
  G4double a = A;
  if(1.0 > A) { a = 1.0/A; }
  if(a <= maxA)
  {
    G4int i = G4int(a + 0.5);
    G4double x = (a/G4double(i) - 1.0)*onethird;
    res = pz13[i]*(1.0 + x - x*x*(1.0 - 1.66666666*x));
    if(1.0 > A) { res = 1.0/res; }
  }
  else
  {
    res = std::pow(A, onethird); 
  }
  return res;
}

inline G4double G4Pow::Z23(G4int Z)
{
  G4double x = Z13(Z);
  return x*x;
}

inline G4double G4Pow::A23(G4double A)
{
  G4double x = A13(A);
  return x*x;
}

inline G4double G4Pow::logZ(G4int Z)
{
  return lz[Z];
}

inline G4double G4Pow::logA(G4double A)
{
  G4double res;
  G4double a = A;
  if(1.0 > A) { a = 1.0/A; }
  if(a <= maxA)
  {
    G4int i = G4int(a + 0.5);
    G4double x = a/G4double(i) - 1;
    res = lz[i] + x*(1.0 - (0.5 - onethird*x)*x);
    if(1.0 > A) { res = -res; }
  }
  else
  {
    res = std::log(A);
  }
  return res;
}

inline G4double G4Pow::log10Z(G4int Z)
{
  return lz[Z]/lz[10];
}

inline G4double G4Pow::log10A(G4double A)
{
  return logA(A)/lz[10];
}

inline G4double G4Pow::powZ(G4int Z, G4double y)
{
  return std::exp(y*lz[Z]);
}

inline G4double G4Pow::powA(G4double A, G4double y)
{
  return std::exp(y*logA(A));
}

inline G4double G4Pow::factorial(G4int Z)
{
  return fact[Z];
}

inline G4double G4Pow::logfactorial(G4int Z)
{
  return logfact[Z];
}

// -------------------------------------------------------------------

#endif

