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
// $Id: G4PolarizationManager.hh 97384 2016-06-02 09:59:17Z gcosmo $
//
// GEANT4 Class header file
//
//
// File name:     G4PolarizationManager
//
// Author:        Andreas Schaelicke
//
// Creation date: 01.05.2005
//
// Modifications:
// 12-08-06   Helper routines are moved to G4PolarizationHelper
//
// Class Description:
//
// Provides polarization information for logical volumes, and some basic 
// transformation routines.

#ifndef G4PolarizationManager_h
#define G4PolarizationManager_h 1

#include "globals.hh"
#include <vector>
#include "G4ThreeVector.hh"

#include <map>

class G4LogicalVolume;
class G4PolarizationMessenger;

typedef std::map<G4LogicalVolume*,G4ThreeVector> PolarizationMap;

class G4PolarizationManager {
public:

  virtual ~G4PolarizationManager();
  static G4PolarizationManager* GetInstance();
  static void Dispose();

  void ListVolumes();
  inline void Clean();

  void SetVolumePolarization(G4LogicalVolume* lVol, const G4ThreeVector & pol);
  void SetVolumePolarization(const G4String & lVolName, const G4ThreeVector & pol);

  inline const G4ThreeVector & GetVolumePolarization(G4LogicalVolume* lVol) const;
  inline bool IsPolarized(G4LogicalVolume* lVol) const;

  inline void SetVerbose(G4int val);
  inline G4int GetVerbose() const;

  inline void SetActivated(G4bool val);
  inline bool IsActivated() const;

private:
  G4PolarizationManager();

  G4ThreeVector zeroPolarization;
  PolarizationMap volumePolarizations;
  G4PolarizationMessenger * messenger;

  G4int verboseLevel;
  G4bool activated;

  static G4ThreadLocal G4PolarizationManager* instance;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline const G4ThreeVector & G4PolarizationManager::GetVolumePolarization(G4LogicalVolume* lVol) const
{
  if (!activated) return zeroPolarization;
  PolarizationMap::const_iterator cit=volumePolarizations.find(lVol);
  if (cit!=volumePolarizations.end()) return cit->second;
  return zeroPolarization;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline bool G4PolarizationManager::IsPolarized(G4LogicalVolume* lVol) const
{
  if (!activated) return false;
  PolarizationMap::const_iterator cit=volumePolarizations.find(lVol);
  return (cit!=volumePolarizations.end());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4PolarizationManager::SetVerbose(G4int val) { verboseLevel = val; }
inline G4int G4PolarizationManager::GetVerbose() const { return verboseLevel; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4PolarizationManager::SetActivated(G4bool val) { activated = val; }
inline bool G4PolarizationManager::IsActivated() const { return activated; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4PolarizationManager::Clean() { volumePolarizations.clear(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

