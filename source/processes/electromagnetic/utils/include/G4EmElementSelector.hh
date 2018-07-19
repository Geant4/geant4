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
// $Id: G4EmElementSelector.hh 95657 2016-02-17 13:03:36Z gcosmo $
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
    for(G4int i=0; i<nElmMinusOne; ++i) {
      if (x <= (xSections[i])->Value(e)) {
        element = (*theElementVector)[i];
        break;
      }
    }
  }
  return element;
}

inline const G4Material* G4EmElementSelector::GetMaterial() const
{
  return material;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#endif

