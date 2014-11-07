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
// $Id: $
//
// ---------------------------------------------------------------
//
// G4ParticleChangeForOccurenceBiasing
//
// Class Description:
//     A G4VParticleChange dedicated to occurence biasing : it
// applies weights for non-interaction over a step, and for
// interaction at the end of the step (if interaction occurs) on
// top of a given particle change, that is the one produced by a
// physics process (biased in its final state or not).
//
// ---------------------------------------------------------------
//   Initial version                         Nov. 2013 M. Verderi
//
#ifndef G4ParticleChangeForOccurenceBiasing_hh
#define G4ParticleChangeForOccurenceBiasing_hh 1

#include "G4VParticleChange.hh"

class G4ParticleChangeForOccurenceBiasing : public G4VParticleChange {
public:
  G4ParticleChangeForOccurenceBiasing(G4String name);
  ~G4ParticleChangeForOccurenceBiasing();

public:
  void     SetOccurenceWeightForNonInteraction(G4double w)       {fOccurenceWeightForNonInteraction = w;}
  G4double GetOccurenceWeightForNonInteraction()           const {return fOccurenceWeightForNonInteraction;}
  void     SetOccurenceWeightForInteraction(G4double w)          {fOccurenceWeightForInteraction = w;}
  G4double GetOccurenceWeightForInteraction()              const {return fOccurenceWeightForInteraction;}

public:
  // -- set a wrapped particle change AND USE IT TO UPDATE this occurence particle change state:
  void               SetWrappedParticleChange(G4VParticleChange* wpc);
  G4VParticleChange* GetWrappedParticleChange()                      const {return fWrappedParticleChange;}
public:
  // -- collect the secondaries from the wrapped particle change, apply weight correction, and clear wrapped particle change:
  void StealSecondaries();
  
public:
  // -- from base class G4VParticleChange:
  virtual G4Step* UpdateStepForAtRest   (G4Step* step);
  virtual G4Step* UpdateStepForAlongStep(G4Step* step);
  virtual G4Step* UpdateStepForPostStep (G4Step* step);

public:
  const G4String& GetName() const {return fName;}

private:
  G4String                                       fName;
  G4VParticleChange*            fWrappedParticleChange;
  G4double           fOccurenceWeightForNonInteraction;
  G4double              fOccurenceWeightForInteraction;

  
};

#endif
