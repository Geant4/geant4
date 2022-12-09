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
// Geant4 Class file
//
// File name:     G4PolarizationManager
//
// Author:        Andreas Schaelicke
//
// Class Description:
//   Provides polarization information for logical volumes, and some basic
//   transformation routines.

#include "G4PolarizationManager.hh"

#include "G4LogicalVolume.hh"
#include "G4PolarizationMessenger.hh"
#include "G4StokesVector.hh"

G4ThreadLocal G4PolarizationManager* G4PolarizationManager::fInstance = nullptr;

G4PolarizationManager* G4PolarizationManager::GetInstance()
{
  if(fInstance == nullptr)
    fInstance = new G4PolarizationManager();
  return fInstance;
}

void G4PolarizationManager::Dispose()
{
  if(fInstance != nullptr)
  {
    delete fInstance;
    fInstance = nullptr;
  }
}

G4PolarizationManager::G4PolarizationManager()
  : fMessenger(nullptr)
  , fVerboseLevel(0)
  , fActivated(true)
{
  fMessenger = new G4PolarizationMessenger(this);
}

G4PolarizationManager::~G4PolarizationManager() {}

void G4PolarizationManager::ListVolumes()
{
  if(fVolumePolarizations.empty())
    return;
  G4cout << " Polarization for " << fVolumePolarizations.size()
         << " registered volume(s) : " << G4endl;
  if(!fActivated)
    G4cout << " but polarization deactivated " << G4endl;
  for(auto& vp : fVolumePolarizations)
  {
    G4cout << vp.first->GetName() << " : " << vp.second << G4endl;
  }
}

void G4PolarizationManager::SetVolumePolarization(G4LogicalVolume* lVol,
                                                  const G4ThreeVector& pol)
{
  fVolumePolarizations[lVol] = pol;
  if(fVerboseLevel >= 1)
    G4cout << " SetVolumePolarization " << lVol->GetName() << " " << pol
           << G4endl;
}

void G4PolarizationManager::SetVolumePolarization(const G4String& lVolName,
                                                  const G4ThreeVector& pol)
{
  for(auto& vp : fVolumePolarizations)
  {
    if(vp.first->GetName() == lVolName)
    {
      vp.second = pol;
      if(fVerboseLevel >= 1)
        G4cout << " SetVolumePolarization " << lVolName << " " << pol << G4endl;
      return;
    }
  }
  G4ExceptionDescription ed;
  ed << " Logical volume '" << lVolName << "'not registered yet.\n"
     << " Please register before using '/polarization/volume/set'\n";
  G4Exception("G4PolarizationManager::SetVolumePolarization", "pol040",
              FatalException, ed);
}
