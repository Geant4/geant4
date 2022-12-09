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
#include <vector>
#include <iostream>

class G4LevelManager 
{

public:
  // levels - vector of nuclear level objects, ground state 
  //          level has NULL pointer
  // energies - list of excitation energies of nuclear levels starting
  //            from the ground state with energy zero 
  // spin - 2J, where J is the full angular momentum of the state 
  explicit G4LevelManager(G4int Z, G4int A, std::size_t nlev,
			  const std::vector<G4double>& energies,
			  const std::vector<G4int>& spin,
			  const std::vector<const G4NucLevel*>& levels); 

  ~G4LevelManager();

  //===================================================================
  // run time inlined const functions
  //===================================================================

  // only in this method there is a check on the vector boundary
  std::size_t NearestLevelIndex(const G4double energy, const std::size_t index=0) const;

  inline std::size_t NumberOfTransitions() const;

  inline const G4NucLevel* GetLevel(const std::size_t i) const;

  inline G4double LevelEnergy(const std::size_t i) const;

  inline G4double MaxLevelEnergy() const;

  inline std::size_t NearestLowEdgeLevelIndex(const G4double energy) const;

  inline const G4NucLevel* NearestLevel(const G4double energy, 
                                        const std::size_t index=0) const;

  inline G4double NearestLevelEnergy(const G4double energy, 
                                     const std::size_t index=0) const;

  inline G4double NearestLowEdgeLevelEnergy(const G4double energy) const;

  // for stable isotopes life time is -1
  inline G4double LifeTime(const std::size_t i) const;

  inline G4int SpinTwo(const std::size_t i) const;

  inline G4int Parity(const std::size_t i) const;

  inline G4int FloatingLevel(const std::size_t i) const;

  inline G4double ShellCorrection() const;

  inline G4double LevelDensity(const G4double U) const;

  const G4String& FloatingType(const std::size_t i) const;

  void StreamInfo(std::ostream& os) const;

  G4LevelManager(const G4LevelManager & right) = delete;  
  const G4LevelManager& operator=(const G4LevelManager &right) = delete;
  G4bool operator==(const G4LevelManager &right) const = delete;
  G4bool operator!=(const G4LevelManager &right) const = delete;

private:

  std::vector<G4double> fLevelEnergy;
  std::vector<G4int> fSpin;
  std::vector<const G4NucLevel*> fLevels;
  
  G4double fShellCorrection;
  G4double fLevelDensity;

  std::size_t nTransitions;

  static const G4int nfloting = 13;
  static G4String fFloatingLevels[nfloting];

};

inline std::size_t G4LevelManager::NumberOfTransitions() const
{
  return nTransitions;
}

inline const G4NucLevel* G4LevelManager::GetLevel(const std::size_t i) const
{
  return fLevels[i]; 
}

inline G4double G4LevelManager::LevelEnergy(const std::size_t i) const
{
  return fLevelEnergy[i]; 
}

inline G4double G4LevelManager::MaxLevelEnergy() const
{
  return fLevelEnergy[nTransitions]; 
}

inline std::size_t 
G4LevelManager::NearestLowEdgeLevelIndex(const G4double energy) const
{
  std::size_t idx = nTransitions;
  if(energy < fLevelEnergy[nTransitions]) {
    idx = std::lower_bound(fLevelEnergy.begin(), fLevelEnergy.end(), energy)
      - fLevelEnergy.begin() - 1;
  }
  return idx;
}

inline const G4NucLevel* 
G4LevelManager::NearestLevel(const G4double energy, const std::size_t index) const
{   
  return GetLevel(NearestLevelIndex(energy, index));
}

inline G4double
G4LevelManager::NearestLevelEnergy(const G4double energy, 
                                   const std::size_t index) const
{   
  return LevelEnergy(NearestLevelIndex(energy, index));
}

inline G4double 
G4LevelManager::NearestLowEdgeLevelEnergy(const G4double energy) const
{   
  return LevelEnergy(NearestLowEdgeLevelIndex(energy));
}

inline G4double G4LevelManager::LifeTime(const std::size_t i) const
{
  return (fLevels[i]) ? fLevels[i]->GetTimeGamma() : 0.0;
}
   
inline G4int G4LevelManager::SpinTwo(const std::size_t i) const
{
  return std::abs(fSpin[i]%100000 - 100); 
}

inline G4int G4LevelManager::Parity(const std::size_t i) const
{
  return (fSpin[i]%100000 - 100 > 0) ? 1 : -1; 
}

inline G4int G4LevelManager::FloatingLevel(const std::size_t i) const
{
  return fSpin[i]/100000; 
}

inline G4double G4LevelManager::ShellCorrection() const
{
  return fShellCorrection;
}

inline G4double G4LevelManager::LevelDensity(const G4double) const
{
  return fLevelDensity;
}

#endif
