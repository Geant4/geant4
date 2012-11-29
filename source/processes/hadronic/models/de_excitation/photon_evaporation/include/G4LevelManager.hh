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
// $Id$
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
// Nuclear level manager for photon de-excitation process 
// Providing management of levels for high energy generators
// The internal conversion is excluded because there are no atomic shell
// 

#ifndef G4LEVELMANAGER_HH
#define G4LEVELMANAGER_HH 1

#include "globals.hh"
#include "G4NucLevel.hh"
#include <vector>

class G4LevelReader;

class G4LevelManager 
{

public:

  G4LevelManager(G4int Z, G4int A, G4LevelReader& reader,
		 const G4String& filename);   

  ~G4LevelManager();
  
  inline G4int NumberOfLevels() const;

  inline const G4NucLevel* GetLevel(G4int i) const;

  inline const G4NucLevel* NearestLevel(G4double energy) const;
  
private:

  G4LevelManager(const G4LevelManager & right);  
  const G4LevelManager& operator=(const G4LevelManager &right);
  G4bool operator==(const G4LevelManager &right) const;
  G4bool operator!=(const G4LevelManager &right) const;

  std::vector<G4NucLevel*> fLevel;
  G4int    nLevels;
  G4int    theZ;
  G4int    theA;
  G4double fEdiffMax;
};

inline G4int G4LevelManager::NumberOfLevels() const
{
  return nLevels;
}

inline const G4NucLevel* G4LevelManager::GetLevel(G4int i) const
{
  return fLevel[i];
}

inline const G4NucLevel* 
G4LevelManager::NearestLevel(G4double energy) const
{   
  const G4NucLevel* p = 0;
  if(energy < fLevel[nLevels-1]->LevelEnergy() + fEdiffMax) {
    for(G4int i=nLevels-1; i>=0; --i) {
      p = fLevel[i];
      G4double lEnergy = p->LevelEnergy();
      if(0 == i || energy > lEnergy) { break; }
      if(lEnergy - energy <= 0.5*(lEnergy - fLevel[i-1]->LevelEnergy()))  
	{ break; }
    }
  }
  return p;
}

#endif
