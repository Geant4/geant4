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
// $Id: G4QTauNuclearCrossSection.hh,v 1.1 2004-03-05 13:23:03 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 physics class: G4QTauNuclearCrossSection -- header file
// M.V. Kossov, CERN-ITEP(Moscow), 4-FEB-2004
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 4-Feb-04
//
// Short description: this G4 class calculates tauNuclear cross section for
// particular Element (GetCrossSection member function)

#ifndef G4QTauNuclearCrossSection_h
#define G4QTauNuclearCrossSection_h 1

#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4QVLeptoNuclearCrossSection.hh"
#include "G4QPhotoNuclearCrossSection.hh"
#include <vector>
#include "Randomize.hh"
#include "G4TauPlus.hh"
#include "G4TauMinus.hh"

class G4QTauNuclearCrossSection : public G4QVLeptoNuclearCrossSection
{
public:

  G4QTauNuclearCrossSection() : G4QVLeptoNuclearCrossSection()  {}

  ~G4QTauNuclearCrossSection() {}

  G4double GetCrossSection(G4double Energy, G4int Z, G4int N);
  G4double GetEquivalentPhotonEnergy();
  G4double GetEquivalentPhotonQ2(G4double nu);

private:
  G4double SolveTheEquation(G4double f);
  G4double Fun(G4double x);
  G4double DFun(G4double x);

// Body
private:
  static G4int     lastN;    // The last N of calculated nucleus
  static G4int     lastZ;    // The last Z of calculated nucleus
  static G4int     lastF;    // Last used in the cross section TheFirstBin
  static G4double* lastJ1;   // Pointer to the last array of the J1 function
  static G4double* lastJ2;   // Pointer to the last array of the J2 function
  static G4double* lastJ3;   // Pointer to the last array of the J3 function
  static G4int     lastL;    // Last used in the cross section TheLastBin
  static G4double  lastE;    // Last used in the cross section Energy
  static G4double  lastTH;   // Last value of the Energy Threshold
  static G4double  lastSig;  // Last value of the Cross Section
  static G4double  lastG;    // Last value of gamma=lnE-ln(me)
  static G4double  lastH;    // Last value of the High energy A-dependence
  static G4double  ml;       // Mass of Tau-lepton
  static G4double  ml2;      // Squared Tau-lepton mass
  static G4double  lml;      // Logarithm of the Tau-lepton mass
};

inline G4double G4QTauNuclearCrossSection::DFun(G4double x) // Parametrization of PhNuCS
{
  static const G4double shd=1.0734;                    // HE PomShadowing(D)
  static const G4double poc=0.0375;                    // HE Pomeron coefficient
  static const G4double pos=16.5;                      // HE Pomeron shift
  static const G4double reg=.11;                       // HE Reggeon slope
  G4double y=exp(x-lastG-lml);                        // y for the x
  G4double flux=lastG*(2.-y*(2.-y))-1.;                // flux factor
  return (poc*(x-pos)+shd*exp(-reg*x))*flux;
}

inline G4double G4QTauNuclearCrossSection::Fun(G4double x) // Integrated PhotroNuc CS
{
  G4double dlg1=lastG+lastG-1.;
  G4double lgoe=lastG/lastE;
  G4double HE2=HighEnergyJ2(x);
  return dlg1*HighEnergyJ1(x)-lgoe*(HE2+HE2-HighEnergyJ3(x)/lastE);
}

#endif
