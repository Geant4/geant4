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
// $Id: G4PolarizationManager.cc 97384 2016-06-02 09:59:17Z gcosmo $
//
// GEANT4 Class file
//
//
// File name:     G4PolarizationManager
//
// Author:        Andreas Schaelicke
//
// Creation date: 01.05.2005
//
// Modifications:
//
// Class Description:
//
// Provides polarization information for logical volumes, and some basic 
// transformation routines.
//
#include "G4PolarizationManager.hh"
#include "G4PolarizationMessenger.hh"
#include "G4StokesVector.hh"

#include "G4LogicalVolume.hh"

G4ThreadLocal G4PolarizationManager * G4PolarizationManager::instance = nullptr;

G4PolarizationManager* G4PolarizationManager::GetInstance()
{
  if (instance == nullptr) instance = new G4PolarizationManager();
  return instance;
}

void G4PolarizationManager::Dispose()
{
  if (instance != nullptr)
  {
    delete instance;
    instance = nullptr;
  }
}

G4PolarizationManager::G4PolarizationManager()
  : messenger(nullptr), verboseLevel(0), activated(true)
{
  messenger = new G4PolarizationMessenger(this);  
}

G4PolarizationManager::~G4PolarizationManager()
{  
}

void G4PolarizationManager::ListVolumes()
{
  if (volumePolarizations.size()==0) return;
  G4cout<<" Polarization for "<<volumePolarizations.size()
	<<" registered volume(s) : "<<G4endl;
  if (!activated) 
    G4cout<<" but polarization deactivated "<<G4endl;
  for (auto vp : volumePolarizations) {
    G4cout << vp.first->GetName() << " : " << vp.second << G4endl;
  }
}

void G4PolarizationManager::SetVolumePolarization(G4LogicalVolume* lVol, const G4ThreeVector & pol)
{
  volumePolarizations[lVol]=pol;
  if (verboseLevel>=1) G4cout<<" SetVolumePolarization "
			     <<lVol->GetName()<<" "
			     <<pol<<G4endl;
}

void G4PolarizationManager::SetVolumePolarization(const G4String & lVolName, const G4ThreeVector & pol)
{
  for (auto& vp : volumePolarizations) {
    if (vp.first->GetName()==lVolName) {
      vp.second=pol;
      if (verboseLevel>=1) G4cout<<" SetVolumePolarization "
				 <<lVolName<<" "
				 <<pol<<G4endl;
      return;
    }
  }
  G4cout<<" logical volume '"<<lVolName<<"'not registerd yet \n"
	<<" please register before using '/polarization/volume/set' "<<G4endl;
}

