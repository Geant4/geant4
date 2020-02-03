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
//
//
// John Allison  22nd August 2018
// Touchable utilities

#include "G4TouchableUtils.hh"

#include "G4TouchablePropertiesScene.hh"
#include "G4TransportationManager.hh"

G4PhysicalVolumeModel::TouchableProperties G4TouchableUtils::FindTouchableProperties
(G4ModelingParameters::PVNameCopyNoPath path)
{
  // If the search is not successful properties.fpTouchablePV will be null.
  G4PhysicalVolumeModel::TouchableProperties properties;
  G4TransportationManager* transportationManager =
  G4TransportationManager::GetTransportationManager ();
  size_t nWorlds = transportationManager->GetNoWorlds();
  std::vector<G4VPhysicalVolume*>::iterator iterWorld =
  transportationManager->GetWorldsIterator();
  for (size_t i = 0; i < nWorlds; ++i, ++iterWorld) {
    G4PhysicalVolumeModel pvModel (*iterWorld);  // Unlimited depth.
    G4ModelingParameters mp;  // Default - no culling.
    pvModel.SetModelingParameters (&mp);
    G4TouchablePropertiesScene scene (&pvModel,path);
    pvModel.DescribeYourselfTo (scene);  // Initiate geometry tree traversal.
    if (scene.GetFoundTouchableProperties().fpTouchablePV) {
      properties = scene.GetFoundTouchableProperties();
      break; // Found, so no need to scan more worlds.
    }
  }
  return properties;
}
