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
// GEANT4 tag $Name: not supported by cvs2svn $
//
// GEANT4 physics class: G4ElectroNuclearCrossSection -- header file
// M.V. Kossov, ITEP(Moscow), 24-OCT-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 25-Sept-03

//A. Dotti March 11 2014
//Back porting optimizations from G4 Version 10
//Using tags: hadr-cross-V09-06-28, -26, -21, -18
//     Not applied removing of check if particle is e-/e+, implemented by Witek
//     Not 100% sure this is actually valid in 9.6. This should be a minor correction

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

//A Cache element
struct cacheEl_t {
    G4int F;
    G4double* J1;
    G4double* J2;
    G4double* J3;
    G4double H;
    G4double TH;
};

class G4ElectroNuclearCrossSection : public G4VCrossSectionDataSet
{
public:

  G4ElectroNuclearCrossSection(const G4String& name = "ElectroNuclearXS");
  virtual ~G4ElectroNuclearCrossSection();

  virtual void CrossSectionDescription(std::ostream&) const;

  virtual G4bool
  IsIsoApplicable(const G4DynamicParticle* aParticle, G4int /*Z*/,
                  G4int /*A*/, const G4Element*, const G4Material*);

  virtual G4double
  GetIsoCrossSection(const G4DynamicParticle* aParticle, 
		     G4int /*Z*/, G4int /*A*/, 
		     const G4Isotope*, const G4Element*, const G4Material*);

  G4double GetEquivalentPhotonEnergy();

  G4double GetVirtualFactor(G4double nu, G4double Q2);

  G4double GetEquivalentPhotonQ2(G4double nu);

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
    G4int currentN;
    G4int currentZ;
     std::map<G4int,cacheEl_t> cache;
     G4int lastUsedKey;
     cacheEl_t* lastUsedCacheEl;
     G4double lastE ; //Last used energy value
     G4double lastSig ; //Last used XS value
     G4double lastG   ; //Last value of gamma=lnE-ln(me)
     G4int    lastL   ; //Last used in the cross section TheLastBin

};



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
