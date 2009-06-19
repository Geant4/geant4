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
// $Id: G4pow.hh,v 1.1 2009-06-19 14:18:18 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:     G4pow
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 23.05.2009
//
// Modifications:
//
// Class Description:
//
// A utility singleton class for the fast computation of log and pow
// functions 
//

#ifndef G4pow_h
#define G4pow_h 1

#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4int ZMAXNUM = 256; 

class G4pow
{

public:

  static G4pow* pow();
  ~G4pow();

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

private:

  G4pow();
  static G4pow* instance;

  G4double onethird;
  G4double* POWERZ13;
  G4double* LOGZ;
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4pow::Z13(G4int Z)
{
  return POWERZ13[Z];
}

inline G4double G4pow::A13(G4double A)
{
  G4double res;
  G4int i = G4int(A + 0.5);
  if(i >= 1 && i < ZMAXNUM) {
    G4double x = (1.0 - A/G4double(i))*onethird;
    res = POWERZ13[i]*(1.0 + x - x*x*(1.0 - 1.66666666*x));
  } else res = std::pow(A, onethird); 
  return res;
}

inline G4double G4pow::Z23(G4int Z)
{
  G4double x = Z13(Z);
  return x*x;
}

inline G4double G4pow::A23(G4double A)
{
  G4double x = A13(A);
  return x*x;
}

inline G4double G4pow::logZ(G4int Z)
{
  return LOGZ[Z];
}

inline G4double G4pow::logA(G4double A)
{
  G4double res;
  G4int i = G4int(A + 0.5);
  if(i >= 1 && i < ZMAXNUM) {
    G4double x = 1.0 - A/G4double(i);
    res = LOGZ[i] + x*(1.0 - 0.5*x + onethird*x*x);
  } else res = std::log(A);   
  return res;
}

inline G4double G4pow::log10Z(G4int Z)
{
  return LOGZ[Z]/LOGZ[10];
}

inline G4double G4pow::log10A(G4double A)
{
  return logA(A)/LOGZ[10];
}

inline G4double G4pow::powZ(G4int Z, G4double y)
{
  return std::exp(y*LOGZ[Z]);
}

inline G4double G4pow::powA(G4double A, G4double y)
{
  return std::exp(y*logA(A));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

