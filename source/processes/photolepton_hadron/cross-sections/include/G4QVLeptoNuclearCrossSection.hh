//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4QVLeptoNuclearCrossSection.hh,v 1.1 2004-03-05 13:23:03 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 physics class: G4QVLeptoNuclearCrossSection -- header file
// M.V. Kossov, CERN-ITEP(Moscow), 4-FEB-2004
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 4-Feb-04
//
// Short description: this G4Q interface-class defines what must be calculated for
// lepto-nuclear interaction with particular Element.

#ifndef G4QVLeptoNuclearCrossSection_h
#define G4QVLeptoNuclearCrossSection_h 1

#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include "G4NucleiPropertiesTable.hh"
#include <vector>
#include "Randomize.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"

class G4QVLeptoNuclearCrossSection
{
public:

  G4QVLeptoNuclearCrossSection(){}
  virtual ~G4QVLeptoNuclearCrossSection() {}

  virtual G4double GetCrossSection(G4double Energy, G4int Z, G4int N) = 0; // Cross-section
  virtual G4double GetEquivalentPhotonEnergy() = 0; // Simulate the nu value
  G4double GetVirtualFactor(G4double nu, G4double Q2); // Form-factor modification
  virtual G4double GetEquivalentPhotonQ2(G4double nu) =0; // Simulate Q2 for known nu
  G4int    GetFunctions(G4double a, G4double* x, G4double* y, G4double* z);
  G4double HighEnergyJ1(G4double lE);
  G4double HighEnergyJ2(G4double lE);
  G4double HighEnergyJ3(G4double lE);

private:
  virtual G4double SolveTheEquation(G4double f) = 0;
  virtual G4double Fun(G4double x) = 0;
  virtual G4double DFun(G4double x) = 0;
};

inline G4double G4QVLeptoNuclearCrossSection::HighEnergyJ1(G4double lEn)
{
  static const G4double le=log(50000.); // log(E0)
  static const G4double le2=le*le;      // log(E0)^2
  static const G4double a=.0375;        // a
  static const G4double ha=a*.5;        // a/2
  static const G4double ab=a*16.5;      // a*b
  static const G4double d=0.11;         // d
  static const G4double cd=1.0734/d;    // c/d
  static const G4double ele=exp(-d*le); // E0^(-d)
  return ha*(lEn*lEn-le2)-ab*(lEn-le)-cd*(exp(-d*lEn)-ele);
}

inline G4double G4QVLeptoNuclearCrossSection::HighEnergyJ2(G4double lEn)
{
  static const G4double e=50000.;       // E0
  static const G4double le=log(e);      // log(E0)
  static const G4double le1=(le-1.)*e;  // (log(E0)-1)*E0
  static const G4double a=.0375;        // a
  static const G4double ab=a*16.5;      // a*b
  static const G4double d=1.-0.11;      // 1-d
  static const G4double cd=1.0734/d;    // c/(1-d)
  static const G4double ele=exp(d*le);  // E0^(1-d)
  G4double En=exp(lEn);
  return a*((lEn-1.)*En-le1)-ab*(En-e)+cd*(exp(d*lEn)-ele);
}

inline G4double G4QVLeptoNuclearCrossSection::HighEnergyJ3(G4double lEn)
{
  static const G4double e=50000.;       // E0
  static const G4double le=log(e);      // log(E0)
  static const G4double e2=e*e;         // E0^2
  static const G4double leh=(le-.5)*e2; // (log(E0)-.5)*E0^2
  static const G4double ha=.0375*.5;    // a/2
  static const G4double hab=ha*16.5;    // a*b/2
  static const G4double d=2.-.11;       // 2-d
  static const G4double cd=1.0734/d;    // c/(2-d)
  static const G4double ele=exp(d*le);  // E0^(2-d)
  G4double lastE2=exp(lEn+lEn);
  return ha*((lEn-.5)*lastE2-leh)-hab*(lastE2-e2)+cd*(exp(d*lEn)-ele);
}

#endif
