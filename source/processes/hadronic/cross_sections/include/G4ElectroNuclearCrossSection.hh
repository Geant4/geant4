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
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 physics class: G4ElectroNuclearCrossSection -- header file
// M.V. Kossov, ITEP(Moscow), 24-OCT-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 25-Sept-03
//

#ifndef G4ElectroNuclearCrossSection_h
#define G4ElectroNuclearCrossSection_h 1

#include "G4VCrossSectionDataSet.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include <vector>
#include "Randomize.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

class G4ElectroNuclearCrossSection : public G4VCrossSectionDataSet
{
public:

  G4ElectroNuclearCrossSection();
  virtual ~G4ElectroNuclearCrossSection();

  G4bool IsApplicable(const G4DynamicParticle* aParticle, const G4Element* )
  {
    return IsIsoApplicable(aParticle, 0, 0);
  }

  G4bool IsIsoApplicable(const G4DynamicParticle* aParticle,
                         G4int /*ZZ*/, G4int /*AA*/)
  {
    G4bool result = false;
    if (aParticle->GetDefinition() == G4Electron::ElectronDefinition())
       result = true;
    if (aParticle->GetDefinition() == G4Positron::PositronDefinition())
       result = true;
    return result;
  }


  G4double GetCrossSection(const G4DynamicParticle* aParticle, 
                           const G4Element* anElement, G4double T=0.);

  G4double GetZandACrossSection(const G4DynamicParticle* aParticle, 
                                G4int ZZ, G4int AA, G4double T=0.);

  G4double GetEquivalentPhotonEnergy();

  G4double GetVirtualFactor(G4double nu, G4double Q2);

  G4double GetEquivalentPhotonQ2(G4double nu);

  void BuildPhysicsTable(const G4ParticleDefinition&) {}

  void DumpPhysicsTable(const G4ParticleDefinition&) {}

private:
  G4int    GetFunctions(G4double a, G4double* x, G4double* y, G4double* z);

  G4double ThresholdEnergy(G4int Z, G4int N);
  G4double HighEnergyJ1(G4double lE);
  G4double HighEnergyJ2(G4double lE);
  G4double HighEnergyJ3(G4double lE);
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

  // Vector of pointers to the J1 tabulated functions
  static std::vector <G4double*> J1;

  // Vector of pointers to the J2 tabulated functions
  static std::vector <G4double*> J2;

  // Vector of pointers to the J3 tabulated functions
  static std::vector <G4double*> J3;
};


inline G4double
G4ElectroNuclearCrossSection::DFun(G4double x)
{
  // Parametrization of the PhotoNucCS
  static const G4double shd=1.0734;              // HE PomShadowing(D)
  static const G4double poc=0.0375;              // HE Pomeron coefficient
  static const G4double pos=16.5;                // HE Pomeron shift
  static const G4double reg=.11;                 // HE Reggeon slope
  static const G4double mel=0.5109989;           // Mass of an electron in MeV
  static const G4double lmel=std::log(mel);      // Log of an electron mass
  G4double y=std::exp(x-lastG-lmel);             // y for the x
  G4double flux=lastG*(2.-y*(2.-y))-1.;          // flux factor
  return (poc*(x-pos)+shd*std::exp(-reg*x))*flux;
}


inline G4double
G4ElectroNuclearCrossSection::Fun(G4double x)
{
  // Integrated PhoNuc cross section
  G4double dlg1=lastG+lastG-1.;
  G4double lgoe=lastG/lastE;
  G4double HE2=HighEnergyJ2(x);
  return dlg1*HighEnergyJ1(x)-lgoe*(HE2+HE2-HighEnergyJ3(x)/lastE);
}


inline G4double
G4ElectroNuclearCrossSection::HighEnergyJ1(G4double lEn)
{
  static const G4double le=std::log(50000.); // std::log(E0)
  static const G4double le2=le*le;      // std::log(E0)^2
  static const G4double a=.0375;        // a
  static const G4double ha=a*.5;        // a/2
  static const G4double ab=a*16.5;      // a*b
  static const G4double d=0.11;         // d
  static const G4double cd=1.0734/d;    // c/d
  static const G4double ele=std::exp(-d*le); // E0^(-d)
  return ha*(lEn*lEn-le2)-ab*(lEn-le)-cd*(std::exp(-d*lEn)-ele);
}


inline G4double
G4ElectroNuclearCrossSection::HighEnergyJ2(G4double lEn)
{
  static const G4double e=50000.;       // E0
  static const G4double le=std::log(e);      // std::log(E0)
  static const G4double le1=(le-1.)*e;  // (std::log(E0)-1)*E0
  static const G4double a=.0375;        // a
  static const G4double ab=a*16.5;      // a*b
  static const G4double d=1.-0.11;      // 1-d
  static const G4double cd=1.0734/d;    // c/(1-d)
  static const G4double ele=std::exp(d*le);  // E0^(1-d)
  G4double En=std::exp(lEn);
  return a*((lEn-1.)*En-le1)-ab*(En-e)+cd*(std::exp(d*lEn)-ele);
}


inline G4double
G4ElectroNuclearCrossSection::HighEnergyJ3(G4double lEn)
{
  static const G4double e=50000.;       // E0
  static const G4double le=std::log(e);      // std::log(E0)
  static const G4double e2=e*e;         // E0^2
  static const G4double leh=(le-.5)*e2; // (std::log(E0)-.5)*E0^2
  static const G4double ha=.0375*.5;    // a/2
  static const G4double hab=ha*16.5;    // a*b/2
  static const G4double d=2.-.11;       // 2-d
  static const G4double cd=1.0734/d;    // c/(2-d)
  static const G4double ele=std::exp(d*le);  // E0^(2-d)
  G4double lastE2=std::exp(lEn+lEn);
  return ha*((lEn-.5)*lastE2-leh)-hab*(lastE2-e2)+cd*(std::exp(d*lEn)-ele);
}

#endif
