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
// Modifications:
// 20110810  M. Kelsey -- Store secondaries by value, not by pointer.
//		Improve constness of argument passing.  Fix up some
//		functions to avoid creating temporaries.

#include "G4HadFinalState.hh"
#include "G4HadronicException.hh"


G4HadFinalState::G4HadFinalState()
  : theDirection(0,0,1), theEnergy(-1), theStat(isAlive), 
    theW(1.), theEDep(0.) {}


G4int G4HadFinalState::GetNumberOfSecondaries() const {return theSecs.size();}

void G4HadFinalState::SetEnergyChange(G4double anEnergy) 
{
  theEnergy=anEnergy;
  if(theEnergy<0) 
    {
      std::cout << "Final state energy was: E = "<<theEnergy<<G4endl;
      throw G4HadronicException(__FILE__, __LINE__, "G4HadFinalState: fatal - negative energy");
    }
}

G4double G4HadFinalState::GetEnergyChange() const {return theEnergy;}

void G4HadFinalState::SetMomentumChange(const G4ThreeVector& aV) {
  theDirection=aV;
}

void G4HadFinalState::SetMomentumChange(G4double x, G4double y, G4double z) 
{
  theDirection.set(x,y,z);
  if(std::fabs(theDirection.mag()-1)>0.001) 
    {
      G4cout <<"We have negative theDirection.mag() = "<<theDirection.mag()<<G4endl;
      throw G4HadronicException(__FILE__, __LINE__, "G4HadFinalState: fatal - negative direction.mag().");
    }
}

const G4ThreeVector& G4HadFinalState::GetMomentumChange() const {return theDirection;}

void G4HadFinalState::AddSecondary(G4DynamicParticle *aP) {
  // NOTE:  In-situ constructor will be optimized away (no copying)
  if (aP) theSecs.push_back(G4HadSecondary(aP));
}

// Concatenate lists efficiently
void G4HadFinalState::AddSecondaries(const std::vector<G4HadSecondary>& addSecs)
{
  theSecs.insert(theSecs.end(),addSecs.begin(),addSecs.end());
}

void G4HadFinalState::SetStatusChange(G4HadFinalStateStatus aS){theStat=aS;}

G4HadFinalStateStatus G4HadFinalState::GetStatusChange() const {return theStat;}

void G4HadFinalState::ClearSecondaries() { 
  theSecs.clear();
}

void G4HadFinalState::Clear()
{
  theDirection.set(0,0,1);
  theEnergy = -1;
  theStat = isAlive;
  theW = 1.;
  theEDep = 0.;
  ClearSecondaries();
}

void G4HadFinalState::SecondariesAreStale() { /*DEPRECATED*/ }

const G4LorentzRotation& G4HadFinalState::GetTrafoToLab() const {return theT;}

void G4HadFinalState::SetTrafoToLab(const G4LorentzRotation & aT) {theT = aT;}

void G4HadFinalState::SetWeightChange(G4double aW){ theW=aW;}

G4double G4HadFinalState::GetWeightChange() const {return theW;}

G4HadSecondary * G4HadFinalState::GetSecondary(size_t i) 
{
  if(i>theSecs.size())
    {
      throw G4HadronicException(__FILE__, __LINE__, 
				"Trying direct access to secondary beyond end of list");
    }
  return &theSecs[i];
}

const G4HadSecondary* G4HadFinalState::GetSecondary(size_t i) const
{
  if(i>theSecs.size())
    {
      throw G4HadronicException(__FILE__, __LINE__, 
				"Trying direct access to secondary beyond end of list");
    }
  return &theSecs[i];
}

void G4HadFinalState::SetLocalEnergyDeposit(G4double aE) {theEDep=aE;}

G4double G4HadFinalState::GetLocalEnergyDeposit() const {return theEDep;}
