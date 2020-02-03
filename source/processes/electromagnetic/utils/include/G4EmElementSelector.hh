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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4EmElementSelector
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 29.05.2008
//
// Modifications:
//
//
// Class Description:
//
// Generic helper class for the random selection of an element

// -------------------------------------------------------------------
//

#ifndef G4EmElementSelector_h
#define G4EmElementSelector_h 1

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4PhysicsLogVector.hh"
#include "Randomize.hh"
#include <vector>

class G4VEmModel;

class G4EmElementSelector
{

public:

  G4EmElementSelector(G4VEmModel*, const G4Material*, G4int bins, 
                      G4double emin, G4double emax, 
                      G4bool spline = true);

  ~G4EmElementSelector();

  void Initialise(const G4ParticleDefinition*, G4double cut = 0.0);

  void Dump(const G4ParticleDefinition* p = nullptr);

  inline const G4Element* SelectRandomAtom(G4double kineticEnergy) const;
  inline const G4Element* SelectRandomAtom(const G4double kineticEnergy,
                                           const G4double logEKin) const;

  inline const G4Material* GetMaterial() const;

private:

  //  hide assignment operator
  G4EmElementSelector & operator=(const  G4EmElementSelector &right) = delete;
  G4EmElementSelector(const  G4EmElementSelector&) = delete;

  G4VEmModel*       model;
  const G4Material* material;
  const G4ElementVector* theElementVector;

  G4int    nElmMinusOne;
  G4int    nbins;

  G4double cutEnergy;
  G4double lowEnergy;
  G4double highEnergy;

  std::vector<G4PhysicsLogVector*> xSections;
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

inline const G4Element* G4EmElementSelector::SelectRandomAtom(G4double e) const
{
  const G4Element* element = (*theElementVector)[nElmMinusOne];
  if (nElmMinusOne > 0) {
    G4double x = G4UniformRand();
    size_t idx(0);
    for(G4int i=0; i<nElmMinusOne; ++i) {
      if (x <= (xSections[i])->Value(e, idx)) {
        element = (*theElementVector)[i];
        break;
      }
    }
  }
  return element;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
inline const G4Element*
G4EmElementSelector::SelectRandomAtom(const G4double e,const G4double loge)const
{
  const G4Element* element = (*theElementVector)[nElmMinusOne];
  if (nElmMinusOne > 0) {
    // 1. Determine energy index (only once)
    const size_t nBins  = (xSections[0])->GetVectorLength();
    // handle cases below/above the enrgy grid (by ekin, idx that gives a=0/1)
    // ekin = x[0]   if e<=x[0]   and idx will be   0 ^ a=0 => so y=y0
    // ekin = x[N-1] if e>=x[N-1] and idx will be N-2 ^ a=1 => so y=y_{N-1}
    const G4double ekin = std::max((xSections[0])->Energy(0),
                                   std::min((xSections[0])->Energy(nBins-1),e));
    // compute the lower index of the bin (idx \in [0,N-2] will be guaranted)
    const size_t    idx = (xSections[0])->ComputeLogVectorBin(loge);
    // 2. Do the linear interp.(robust for corner cases through ekin, idx and a)
    const G4double   x1 = (xSections[0])->Energy(idx);
    const G4double   x2 = (xSections[0])->Energy(idx+1);
    // note: all corner cases of the previous methods are covered and eventually
    //       gives a=0/1 that results in y=y0\y_{N-1} if e<=x[0]/e>=x[N-1] or
    //       y=y_i/y_{i+1} if e<x[i]/e>=x[i+1] due to small numerical errors
    const G4double    a = std::max(0., std::min(1., (ekin - x1)/(x2 - x1)));
    const G4double urnd = G4UniformRand();
    for (G4int i = 0; i < nElmMinusOne; ++i) {
      const G4double  y1 = (*xSections[i])[idx];
      const G4double  y2 = (*xSections[i])[idx+1];
      if (urnd <= y1 + a*(y2-y1)) {
        element = (*theElementVector)[i];
        break;
      }
    }
  }
  return element;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

inline const G4Material* G4EmElementSelector::GetMaterial() const
{
  return material;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#endif

