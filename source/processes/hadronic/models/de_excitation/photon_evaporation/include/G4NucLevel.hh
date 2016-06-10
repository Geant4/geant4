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
// $Id: G4NucLevel.hh 67983 2013-03-13 10:42:03Z gcosmo $
//
// -------------------------------------------------------------------
//
//      GEANT4 header file 
//
//      File name:     G4NucLevel
//
//      Author:        V.Ivanchenko
// 
//      Creation date: 4 January 2012
//
//      Modifications:
//      
// -------------------------------------------------------------------
//
//  Container class keeping information about gamma transition
//  and for a given nuclear level
//

#ifndef G4NUCLEVEL_HH
#define G4NUCLEVEL_HH 1

#include "globals.hh"
#include "Randomize.hh"
#include <vector>

class G4NucLevel 
{
public:

  G4NucLevel(G4double energy, G4double halfLife,
	     const std::vector<G4double>& eGamma,
	     const std::vector<G4double>& wGamma);

  ~G4NucLevel();

  inline G4double LevelEnergy() const;

  inline G4double LevelHalfLife() const;

  inline G4double SampleEnergy() const;

private:  

  G4NucLevel(const G4NucLevel &right);
  G4bool operator==(const G4NucLevel &right) const;
  G4bool operator!=(const G4NucLevel &right) const;
  G4bool operator<(const G4NucLevel &right) const;
  const G4NucLevel& operator=(const G4NucLevel &right);
  
  std::vector<G4double> fTransitionEnergy;
  std::vector<G4double> fCumProbability;

  G4double fEnergy;
  G4double fHalfLifeTime;
  size_t   nTransitions;
};

inline G4double G4NucLevel::LevelEnergy() const
{
  return fEnergy;
}

inline G4double G4NucLevel::LevelHalfLife() const
{
  return fHalfLifeTime;
}

inline G4double G4NucLevel::SampleEnergy() const
{
  G4double e = 0.0;
  G4double x = G4UniformRand();
  for(size_t i=0; i<nTransitions; ++i) {
    if(x < fCumProbability[i]) {
      e = fTransitionEnergy[i];
      break;
    }
  }
  return e;
}


#endif






