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
// G4Pow
//
// Class description:
//
// Utility singleton class for the fast computation of log and pow
// functions. Integer argument should be in the interval 0-512, no
// check is performed inside these methods for performance reasons.
// For factorial integer argument should be in the interval 0-170
// Computations with double arguments are fast for the interval
// 0.002-511.5 for all functions except exponent, which is computed
// for the interval 0-84.4, standard library is used in the opposite case

// Author: Vladimir Ivanchenko, 23.05.2009
// --------------------------------------------------------------------
#ifndef G4Pow_hh
#define G4Pow_hh 1

#include "G4DataVector.hh"
#include "G4Exp.hh"
#include "G4Log.hh"
#include "globals.hh"

class G4Pow
{
 public:
  static G4Pow* GetInstance();
  ~G4Pow() = default;

  // Fast computation of Z^1/3
  //
  inline G4double Z13(G4int Z) const;
  G4double A13(G4double A) const;

  // Fast computation of Z^2/3
  //
  inline G4double Z23(G4int Z) const;
  inline G4double A23(G4double A) const;

  // Fast computation of log(Z)
  //
  inline G4double logZ(G4int Z) const;
  inline G4double logA(G4double A) const;
  inline G4double logX(G4double x) const;

  // Fast computation of log10(Z)
  //
  inline G4double log10Z(G4int Z) const;
  inline G4double log10A(G4double A) const;

  // Fast computation of exp(X)
  //
  inline G4double expA(G4double A) const;

  // Fast computation of pow(Z,X)
  //
  inline G4double powZ(G4int Z, G4double y) const;
  inline G4double powA(G4double A, G4double y) const;
  G4double powN(G4double x, G4int n) const;

  // Fast factorial
  //
  inline G4double factorial(G4int Z) const;
  inline G4double logfactorial(G4int Z) const;

 private:
  G4Pow();

  G4double A13Low(const G4double, const G4bool) const;
  G4double A13High(const G4double, const G4bool) const;

  inline G4double logBase(G4double x) const;

  static G4Pow* fpInstance;

  const G4double onethird = 1.0 / 3.0;
  const G4int max2        = 5;

  G4double maxA;
  G4double maxLowA;
  G4double maxA2;
  G4double maxAexp;

  G4DataVector ener;
  G4DataVector logen;
  G4DataVector pz13;
  G4DataVector lowa13;
  G4DataVector lz;
  G4DataVector lz2;
  G4DataVector fexp;
  G4DataVector fact;
  G4DataVector logfact;
};

// -----------------------------
// Inline methods implementation
// -----------------------------

inline G4double G4Pow::Z13(G4int Z) const { return pz13[Z]; }

inline G4double G4Pow::Z23(G4int Z) const
{
  G4double x = Z13(Z);
  return x * x;
}

inline G4double G4Pow::A23(G4double A) const
{
  G4double x = A13(A);
  return x * x;
}

inline G4double G4Pow::logZ(G4int Z) const { return lz[Z]; }

inline G4double G4Pow::logBase(G4double a) const
{
  G4double res;
  if(a <= maxA2)
  {
    G4int i = G4int(max2 * (a - 1) + 0.5);
    if(i > max2)
    {
      i = max2;
    }
    G4double x = a / (G4double(i) / max2 + 1) - 1;
    res        = lz2[i] + x * (1.0 - (0.5 - onethird * x) * x);
  }
  else if(a <= maxA)
  {
    G4int i    = G4int(a + 0.5);
    G4double x = a / G4double(i) - 1;
    res        = lz[i] + x * (1.0 - (0.5 - onethird * x) * x);
  }
  else
  {
    res = G4Log(a);
  }
  return res;
}

inline G4double G4Pow::logA(G4double A) const
{
  return (1.0 <= A ? logBase(A) : -logBase(1. / A));
}

inline G4double G4Pow::logX(G4double x) const
{
  G4double res = 0.0;
  G4double a   = (1.0 <= x) ? x : 1.0 / x;

  if(a <= maxA)
  {
    res = logBase(a);
  }
  else if(a <= ener[2])
  {
    res = logen[1] + logBase(a / ener[1]);
  }
  else if(a <= ener[3])
  {
    res = logen[2] + logBase(a / ener[2]);
  }
  else
  {
    res = G4Log(a);
  }

  if(1.0 > x)
  {
    res = -res;
  }
  return res;
}

inline G4double G4Pow::log10Z(G4int Z) const { return lz[Z] / lz[10]; }

inline G4double G4Pow::log10A(G4double A) const { return logX(A) / lz[10]; }

inline G4double G4Pow::expA(G4double A) const
{
  G4double res;
  G4double a = (0.0 <= A) ? A : -A;

  if(a <= maxAexp)
  {
    G4int i    = G4int(2 * a + 0.5);
    G4double x = a - i * 0.5;
    res        = fexp[i] * (1.0 + x * (1.0 + 0.5 * (1.0 + onethird * x) * x));
  }
  else
  {
    res = G4Exp(a);
  }
  if(0.0 > A)
  {
    res = 1.0 / res;
  }
  return res;
}

inline G4double G4Pow::powZ(G4int Z, G4double y) const
{
  return expA(y * lz[Z]);
}

inline G4double G4Pow::powA(G4double A, G4double y) const
{
  return (0.0 == A ? 0.0 : expA(y * logX(A)));
}

inline G4double G4Pow::factorial(G4int Z) const { return fact[Z]; }

inline G4double G4Pow::logfactorial(G4int Z) const { return logfact[Z]; }

#endif
