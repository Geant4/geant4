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
// $Id$
//
//
// GEANT4 physics class: G4QMuonNuclearCrossSection -- header file
// M.V. Kossov, CERN-ITEP(Moscow), 4-FEB-2004
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 4-Feb-04
//
// Short description: this G4 singletone class calculates muonNuclear cross section for
// the particular isotope (GetCrossSection member function)
// ****************************************************************************************
// ********* This HEADER is temporary moved from the photolepton_hadron directory *********
// ******* DO NOT MAKE ANY CHANGE! With time it'll move back to photolepton...(M.K.) ******
// ****************************************************************************************
// Short description: reaction cross-sections for muon-nuclear reactions, which
// are integrals over virtual equivalent photons photons. The muon-nuclear GHAD
// model (not CHIPS) gives 2-3 times smaller scattering angle and deposited energy.
// --------------------------------------------------------------------------------

#ifndef G4QMuonNuclearCrossSection_h
#define G4QMuonNuclearCrossSection_h 1

#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include <vector>
#include "Randomize.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4VQCrossSection.hh"

class G4QMuonNuclearCrossSection : public G4VQCrossSection
{
protected:

  G4QMuonNuclearCrossSection()  {};

public:

  ~G4QMuonNuclearCrossSection();

  static G4VQCrossSection* GetPointer(); // Gives a pointer to this singletone

  G4double ThresholdEnergy(G4int Z, G4int N, G4int PDG=13);

  // At present momentum (pMom) must be in GeV (@@ Units)
  virtual G4double GetCrossSection(G4bool fCS, G4double pMom, G4int tgZ, G4int tgN,
                                                                             G4int pPDG=0);

  G4double CalculateCrossSection(G4bool CS, G4int F, G4int I, G4int PDG, G4int Z, G4int N,
                                                                        G4double Momentum);

  G4int    GetExchangePDGCode();

  G4double GetExchangeEnergy();

  G4double GetVirtualFactor(G4double nu, G4double Q2);

  G4double GetExchangeQ2(G4double nu);

private:
  G4int    GetFunctions(G4double a, G4double* x, G4double* y, G4double* z);
  G4double HighEnergyJ1(G4double lE);
  G4double HighEnergyJ2(G4double lE);
  G4double HighEnergyJ3(G4double lE);
  G4double SolveTheEquation(G4double f);
  G4double Fun(G4double x);
  G4double DFun(G4double x);

// Body
private:
  static G4bool    onlyCS;   // flag to calculate only CS
  static G4double  lastSig;  // Last calculated cross section
  static G4int     lastL;    // Last used in the cross section TheLastBin
  static G4double  lastE;    // Last used in the cross section Energy
  static G4int     lastF;    // Last used in the cross section TheFirstBin
  static G4double  lastG;    // Last value of gamma=lnE-ln(m)
  static G4double  lastH;    // Last value of the High energy A-dependence
  static G4double* lastJ1;   // Pointer to the last array of the J1 function
  static G4double* lastJ2;   // Pointer to the last array of the J2 function
  static G4double* lastJ3;   // Pointer to the last array of the J3 function
  static G4int     lastPDG;  // The last projectile PDG
  static G4int     lastN;    // The last N of calculated nucleus
  static G4int     lastZ;    // The last Z of calculated nucleus
  static G4double  lastP;    // Last used in the cross section Momentum
  static G4double  lastTH;   // Last value of the Momentum Threshold
  static G4double  lastCS;   // Last value of the Cross Section
  static G4int     lastI;    // The last position in the DAMDB
  static std::vector <G4double*>* J1; // Vector of pointers to the J1 tabulated functions
  static std::vector <G4double*>* J2; // Vector of pointers to the J2 tabulated functions
  static std::vector <G4double*>* J3; // Vector of pointers to the J3 tabulated functions
};

inline G4double G4QMuonNuclearCrossSection::DFun(G4double x)// Parametrization of PhotNucCS
{
  static const G4double shd=1.0734;                    // HE PomShadowing(D)
  static const G4double poc=0.0375;                    // HE Pomeron coefficient
  static const G4double pos=16.5;                      // HE Pomeron shift
  static const G4double reg=.11;                       // HE Reggeon slope
  static const G4double mmu=105.65839;                 // Mass of a muon in MeV
  static const G4double lmmu=std::log(mmu);                 // Log of a muon mass
  G4double y=std::exp(x-lastG-lmmu);                        // y for the x
  G4double flux=lastG*(2.-y*(2.-y))-1.;                // flux factor
  return (poc*(x-pos)+shd*std::exp(-reg*x))*flux;
}

inline G4double G4QMuonNuclearCrossSection::Fun(G4double x) // Integrated PhoNucCS
{
  G4double dlg1=lastG+lastG-1.;
  G4double lgoe=lastG/lastE;
  G4double HE2=HighEnergyJ2(x);
  return dlg1*HighEnergyJ1(x)-lgoe*(HE2+HE2-HighEnergyJ3(x)/lastE);
}

inline G4double G4QMuonNuclearCrossSection::HighEnergyJ1(G4double lEn)
{
  static const G4double le=std::log(50000.); // log(E0)
  static const G4double le2=le*le;           // log(E0)^2
  static const G4double a=.0375;             // a
  static const G4double ha=a*.5;             // a/2
  static const G4double ab=a*16.5;           // a*b
  static const G4double d=0.11;              // d
  static const G4double cd=1.0734/d;         // c/d
  static const G4double ele=std::exp(-d*le); // E0^(-d)
  return ha*(lEn*lEn-le2)-ab*(lEn-le)-cd*(std::exp(-d*lEn)-ele);
}

inline G4double G4QMuonNuclearCrossSection::HighEnergyJ2(G4double lEn)
{
  static const G4double e=50000.;            // E0
  static const G4double le=std::log(e);      // log(E0)
  static const G4double le1=(le-1.)*e;       // (log(E0)-1)*E0
  static const G4double a=.0375;             // a
  static const G4double ab=a*16.5;           // a*b
  static const G4double d=1.-0.11;           // 1-d
  static const G4double cd=1.0734/d;         // c/(1-d)
  static const G4double ele=std::exp(d*le);  // E0^(1-d)
  G4double En=std::exp(lEn);
  return a*((lEn-1.)*En-le1)-ab*(En-e)+cd*(std::exp(d*lEn)-ele);
}

inline G4double G4QMuonNuclearCrossSection::HighEnergyJ3(G4double lEn)
{
  static const G4double e=50000.;            // E0
  static const G4double le=std::log(e);      // log(E0)
  static const G4double e2=e*e;              // E0^2
  static const G4double leh=(le-.5)*e2;      // (log(E0)-.5)*E0^2
  static const G4double ha=.0375*.5;         // a/2
  static const G4double hab=ha*16.5;         // a*b/2
  static const G4double d=2.-.11;            // 2-d
  static const G4double cd=1.0734/d;         // c/(2-d)
  static const G4double ele=std::exp(d*le);  // E0^(2-d)
  G4double lastE2=std::exp(lEn+lEn);
  return ha*((lEn-.5)*lastE2-leh)-hab*(lastE2-e2)+cd*(std::exp(d*lEn)-ele);
}

#endif
