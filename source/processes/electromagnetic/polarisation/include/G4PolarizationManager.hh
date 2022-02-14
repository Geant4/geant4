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
// Geant44 Class header file
//
// File name:     G4PolarizationManager
//
// Author:        Andreas Schaelicke
//
// Class Description:
//   Provides polarization information for logical volumes

#ifndef G4PolarizationManager_h
#define G4PolarizationManager_h 1

#include "globals.hh"
#include "G4StokesVector.hh"
#include "G4ThreeVector.hh"

#include <vector>
#include <map>

class G4LogicalVolume;
class G4PolarizationMessenger;

typedef std::map<G4LogicalVolume*, G4ThreeVector> PolarizationMap;

class G4PolarizationManager
{
 public:
  ~G4PolarizationManager();
  static G4PolarizationManager* GetInstance();
  static void Dispose();

  void ListVolumes();
  inline void Clean();

  void SetVolumePolarization(G4LogicalVolume* lVol, const G4ThreeVector& pol);
  void SetVolumePolarization(const G4String& lVolName,
                             const G4ThreeVector& pol);

  inline const G4StokesVector GetVolumePolarization(
    G4LogicalVolume* lVol) const;
  inline bool IsPolarized(G4LogicalVolume* lVol) const;

  inline void SetVerbose(G4int val);
  inline G4int GetVerbose() const;

  inline void SetActivated(G4bool val);
  inline bool IsActivated() const;

  G4PolarizationManager& operator=(const G4PolarizationManager& right) = delete;
  G4PolarizationManager(const G4PolarizationManager&)                  = delete;

 private:
  G4PolarizationManager();
  static G4ThreadLocal G4PolarizationManager* fInstance;

  G4PolarizationMessenger* fMessenger;

  PolarizationMap fVolumePolarizations;

  G4int fVerboseLevel;

  G4bool fActivated;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline const G4StokesVector G4PolarizationManager::GetVolumePolarization(
  G4LogicalVolume* lVol) const
{
  if(!fActivated)
    return G4StokesVector::ZERO;
  PolarizationMap::const_iterator cit = fVolumePolarizations.find(lVol);
  if(cit != fVolumePolarizations.end())
  {
    const G4StokesVector vec = G4StokesVector(cit->second);
    return vec;
  }
  return G4StokesVector::ZERO;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline bool G4PolarizationManager::IsPolarized(G4LogicalVolume* lVol) const
{
  if(!fActivated)
    return false;
  PolarizationMap::const_iterator cit = fVolumePolarizations.find(lVol);
  return (cit != fVolumePolarizations.end());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline void G4PolarizationManager::SetVerbose(G4int val)
{
  fVerboseLevel = val;
}
inline G4int G4PolarizationManager::GetVerbose() const { return fVerboseLevel; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline void G4PolarizationManager::SetActivated(G4bool val)
{
  fActivated = val;
}
inline bool G4PolarizationManager::IsActivated() const { return fActivated; }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4PolarizationManager::Clean() { fVolumePolarizations.clear(); }

#endif
