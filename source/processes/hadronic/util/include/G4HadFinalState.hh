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
#ifndef G4HadFinalState_hh
#define G4HadFinalState_hh

// Modifications:
// 20110810  M. Kelsey -- Store secondaries by value, not by pointer.
//		Improve constness of argument passing.  Drop stale flag.
// 20110907  M. Kelsey -- Improve constness of functions, discard add pointer,
//		add function to concatenate vectors

#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4HadSecondary.hh"
#include "G4LorentzRotation.hh"
#include "G4ThreeVector.hh"
#include <vector>

enum G4HadFinalStateStatus{isAlive, stopAndKill, suspend};


class G4HadFinalState
{
public:
  G4HadFinalState();

  inline
  G4int GetNumberOfSecondaries() const           { return theSecs.size(); }
  void SetEnergyChange(G4double anEnergy);
  inline
  G4double GetEnergyChange() const               { return theEnergy; }
  inline
  void SetMomentumChange(const G4ThreeVector& aV){ theDirection=aV; }
  void SetMomentumChange(G4double x, G4double y, G4double z);
  inline
  const G4ThreeVector& GetMomentumChange() const { return theDirection; }
  inline
  void AddSecondary(G4DynamicParticle* aP, G4int mod=-1) {
    theSecs.push_back(G4HadSecondary(aP, theW, mod));
  };
  inline
  void AddSecondary(const G4HadSecondary& aHS)   { theSecs.push_back(aHS); }
  inline
  void SetStatusChange(G4HadFinalStateStatus aS) { theStat=aS; }
  inline
  G4HadFinalStateStatus GetStatusChange() const  { return theStat; }
  void Clear();
  inline
  const G4LorentzRotation& GetTrafoToLab() const { return theT; }
  inline
  void SetTrafoToLab(const G4LorentzRotation & aT) { theT = aT; }
  inline
  void SetWeightChange(G4double aW)              { theW=aW; }
  inline
  G4double GetWeightChange() const               { return theW; }
  G4HadSecondary* GetSecondary(size_t i);
  const G4HadSecondary* GetSecondary(size_t i) const;
  inline
  void SetLocalEnergyDeposit(G4double aE)        { theEDep=aE; }
  inline
  G4double GetLocalEnergyDeposit() const         { return theEDep;}
  //  void SecondariesAreStale();    // Deprecated; not needed for values
  inline
  void ClearSecondaries()                        { theSecs.clear(); }

  // Concatenate lists efficiently
  void AddSecondaries(const std::vector<G4HadSecondary>& addSecs);
  inline
  void AddSecondaries(const G4HadFinalState& addHFS) {
    AddSecondaries(addHFS.theSecs);
  }
  inline
  void AddSecondaries(const G4HadFinalState* addHFS) {
    if (addHFS) AddSecondaries(addHFS->theSecs);
  }

private:
  G4ThreeVector theDirection;
  G4double theEnergy;
  std::vector<G4HadSecondary> theSecs;
  G4HadFinalStateStatus theStat;
  G4LorentzRotation theT;
  G4double theW;
  G4double theEDep;
};

#endif
