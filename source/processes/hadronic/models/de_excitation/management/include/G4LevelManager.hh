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
// $Id: G4LevelManager.hh 88516 2015-02-25 11:00:16Z vnivanch $
//
// -------------------------------------------------------------------
//
//      GEANT4 header file 
//
//      File name:     G4LevelManager
//
//      Author:        V.Ivanchenko
// 
//      Creation date: 4 January 2012
//
//      Modifications:
//  13.02.2015 Design change for gamma de-excitation 
//      
// -------------------------------------------------------------------
//
// Nuclear level manager for photon de-excitation process 
// 

#ifndef G4LEVELMANAGER_HH
#define G4LEVELMANAGER_HH 1

#include "globals.hh"
#include "G4NucLevel.hh"
#include <assert.h>
#include <vector>

static const G4float tolerance = 0.0001f;

class G4LevelManager 
{

public:
  // levels - vector of nuclear level objects, ground state 
  //          level has NULL pointer
  // energies - list of excitation energies of nuclear levels starting
  //            from the ground state with energy zero 
  G4LevelManager(const std::vector<G4float>& energies,
		 const std::vector<G4float>& lifetime,
		 const std::vector<G4float>& lifetimegamma,
		 const std::vector<G4int>& spin,
		 const std::vector<const G4NucLevel*>& levels); 

  ~G4LevelManager();

  //===================================================================
  // run time inlined const functions
  //===================================================================

  inline size_t NumberOfTransitions() const;

  inline const G4NucLevel* GetLevel(size_t i) const;

  inline G4float LevelEnergy(size_t i) const;

  inline G4float MaxLevelEnergy() const;

  inline size_t NearestLevelIndex(G4double energy, size_t index=0) const;

  inline const G4NucLevel* NearestLevel(G4double energy, size_t index=0) const;

  inline G4float NearestLevelEnergy(G4double energy, size_t index=0) const;

  inline G4float LifeTime(size_t i) const;

  inline G4float LifeTimeGamma(size_t i) const;
   
  inline G4int Spin(size_t i) const;

private:

  G4LevelManager(const G4LevelManager & right);  
  const G4LevelManager& operator=(const G4LevelManager &right);
  G4bool operator==(const G4LevelManager &right) const;
  G4bool operator!=(const G4LevelManager &right) const;

  const std::vector<G4float>  fLevelEnergy;
  const std::vector<G4float>  fLifeTime;
  const std::vector<G4float>  fLifeTimeGamma;
  const std::vector<G4int>    fSpin;
  const std::vector<const G4NucLevel*> fLevels;
  size_t nTransitions;

};

inline size_t G4LevelManager::NumberOfTransitions() const
{
  return nTransitions;
}

inline const G4NucLevel* G4LevelManager::GetLevel(size_t i) const
{
  assert(i <= nTransitions);
  return fLevels[i]; 
}

inline G4float G4LevelManager::LevelEnergy(size_t i) const
{
  assert(i <= nTransitions);
  return fLevelEnergy[i]; 
}

inline G4float G4LevelManager::MaxLevelEnergy() const
{
  return fLevelEnergy[nTransitions]; 
}

inline size_t  
G4LevelManager::NearestLevelIndex(G4double ener, size_t index) const
{
  //G4cout << "index= " << index << " max= " << nTransitions << " exc= " << ener 
  //	 << " Emax= " << fLevelEnergy[nTransitions] << G4endl;
  size_t idx = std::min(index, nTransitions);
  G4float energy = (G4float)ener;
  if(0 < nTransitions && std::fabs(energy - fLevelEnergy[idx]) > tolerance) {
    // ground state
    if(energy <= fLevelEnergy[1]*0.5f)  
      { idx = 0; }
    else if((fLevelEnergy[nTransitions] + fLevelEnergy[nTransitions-1])*0.5f <= energy) 
      { idx = nTransitions; }

    // if shortcuts are not working, make binary search
    else {
      idx = std::lower_bound(fLevelEnergy.begin(), fLevelEnergy.end(), energy)
          - fLevelEnergy.begin() - 1;
      if(energy - fLevelEnergy[idx] > fLevelEnergy[idx+1] - energy) { ++idx; }
      //G4cout << "E= " << energy << " " << fLevelEnergy[idx-1] 
      //<< " " << fLevelEnergy[idx] << G4endl;
    }
  }
  return idx;
}

inline const G4NucLevel* 
G4LevelManager::NearestLevel(G4double energy, size_t index) const
{   
  return GetLevel(NearestLevelIndex(energy, index));
}

inline G4float
G4LevelManager::NearestLevelEnergy(G4double energy, size_t index) const
{   
  return LevelEnergy(NearestLevelIndex(energy, index));
}

inline G4float G4LevelManager::LifeTime(size_t i) const
{
  assert(i <= nTransitions);
  return fLifeTime[i]; 
}

inline G4float G4LevelManager::LifeTimeGamma(size_t i) const
{
  assert(i <= nTransitions);
  return fLifeTimeGamma[i]; 
}
   
inline G4int G4LevelManager::Spin(size_t i) const
{
  assert(i <= nTransitions);
  return fSpin[i]; 
}

#endif
