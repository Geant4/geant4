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
//J.M. Quesada (August2008). Based on:
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// 17-11-2010 V.Ivanchenko in constructor replace G4VEmissionProbability by 
//            G4EvaporationProbability and do not new and delete probability
//            object at each call; use G4Pow
// 16-07-2019 V.Ivanchenko use C++11 

#ifndef G4EvaporationChannel_h
#define G4EvaporationChannel_h 1

#include "G4VEvaporationChannel.hh"

class G4EvaporationProbability;
class G4CoulombBarrier;
class G4NuclearLevelData;

class G4EvaporationChannel : public G4VEvaporationChannel
{
public:

  explicit G4EvaporationChannel(G4int A, G4int Z, 
		                G4EvaporationProbability*);

  ~G4EvaporationChannel() override;

  void Initialise() override;

  G4double GetEmissionProbability(G4Fragment* fragment) override; 
  
  G4Fragment* EmittedFragment(G4Fragment* theNucleus) override;

private: 
  
  G4EvaporationChannel(const G4EvaporationChannel & right);
  const G4EvaporationChannel & operator=(const G4EvaporationChannel & right);
  G4bool operator==(const G4EvaporationChannel & right) const;
  G4bool operator!=(const G4EvaporationChannel & right) const;

  // This data member define the channel. 
  // They are intializated at object creation (constructor) time.
  G4int theA;
  G4int theZ;
  G4int resA;
  G4int resZ;

  G4double mass;
  G4double resMass;

  G4double evapMass;
  G4double evapMass2;

  // For evaporation probability calcualation
  G4EvaporationProbability* theProbability;

  // For Coulomb Barrier calculation
  G4CoulombBarrier* theCoulombBarrier;

  G4NuclearLevelData* theLevelData;
};

#endif
