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
// GEANT4 physics class: G4ElectroNuclearCrossSection -- header file
// M.V. Kossov, ITEP(Moscow), 24-OCT-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 25-Sept-03

// A. Dotti 2-May-2013: Re-write caching logic and optimizations

#ifndef G4ElectroNuclearCrossSection_h
#define G4ElectroNuclearCrossSection_h 1

#include "G4VCrossSectionDataSet.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include "G4NistManager.hh"
#include <vector>
#include "Randomize.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include <map>

//A cache element
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

  G4ElectroNuclearCrossSection();
  virtual ~G4ElectroNuclearCrossSection();
    
  static const char* Default_Name() {return "ElectroNuclearXS";}

  virtual void CrossSectionDescription(std::ostream&) const;

  virtual G4bool IsElementApplicable(const G4DynamicParticle*, G4int Z,
                                     const G4Material*);
  virtual G4double GetElementCrossSection(const G4DynamicParticle*, G4int Z,
                                    const G4Material* mat);

  G4double GetEquivalentPhotonEnergy();

  G4double GetVirtualFactor(G4double nu, G4double Q2);

  G4double GetEquivalentPhotonQ2(G4double nu);

private:
  G4int    GetFunctions(G4double a, G4double* x, G4double* y, G4double* z);

  G4double ThresholdEnergy(G4int Z, G4int N);
  G4double SolveTheEquation(G4double f);
  G4double Fun(G4double x);
  G4double DFun(G4double x);
    
  G4double HighEnergyJ1(G4double lE);
  G4double HighEnergyJ2(G4double lE, G4double E);
  G4double HighEnergyJ3(G4double lE, G4double E2);

// Body
private:
    G4int currentN;
    G4int currentZ;
    
    //Cache structure
    G4int lastZ;
    std::vector<cacheEl_t*> cache;
    cacheEl_t* lastUsedCacheEl;
    G4NistManager* nistmngr;
        
    //Cache values for XS
    G4double lastE ; //Last used energy value
    G4double lastSig; //Last used XS value
    G4double lastG; //Last value of gamma=lnE-ln(me)
    G4int    lastL; //Last used in the cross section TheLastBin
    
    const G4double mNeut;
    const G4double mProt;
};

#endif
