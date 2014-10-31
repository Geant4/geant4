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

void G4HadFinalState::SetEnergyChange(G4double anEnergy) 
{
  theEnergy=anEnergy;
  if(theEnergy<0) 
    {
      std::cout << "Final state energy was: E = "<<theEnergy<<G4endl;
      throw G4HadronicException(__FILE__, __LINE__, 
      "G4HadFinalState: fatal - negative energy");
    }
}

void G4HadFinalState::SetMomentumChange(G4double x, G4double y, G4double z) 
{
  theDirection.set(x,y,z);
  if(std::fabs(x*x + y*y + z*z - 1.0)>0.001) {
    G4cout <<"We have negative theDirection.mag() = "<<theDirection.mag()
	   <<G4endl;
    throw G4HadronicException(__FILE__, __LINE__, 
	  "G4HadFinalState: fatal - negative direction.mag().");
  }
}

// Concatenate lists efficiently
void G4HadFinalState::AddSecondaries(const std::vector<G4HadSecondary>& addSecs)
{
  theSecs.insert(theSecs.end(),addSecs.begin(),addSecs.end());
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

//void G4HadFinalState::SecondariesAreStale() { /*DEPRECATED*/ }

G4HadSecondary * G4HadFinalState::GetSecondary(size_t i) 
{
  if(i>theSecs.size()) {
    throw G4HadronicException(__FILE__, __LINE__, 
	  "Trying direct access to secondary beyond end of list");
  }
  return &theSecs[i];
}

const G4HadSecondary* G4HadFinalState::GetSecondary(size_t i) const
{
  if(i>theSecs.size()) {
    throw G4HadronicException(__FILE__, __LINE__, 
	  "Trying direct access to secondary beyond end of list");
  }
  return &theSecs[i];
}
