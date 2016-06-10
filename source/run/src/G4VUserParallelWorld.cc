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
// $Id: G4VUserParallelWorld.cc 69917 2013-05-17 13:28:02Z gcosmo $
//

#include "G4VUserParallelWorld.hh"

#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VSensitiveDetector.hh"
#include "G4SDManager.hh"

G4VUserParallelWorld::G4VUserParallelWorld(G4String worldName)
{ fWorldName = worldName; }

G4VUserParallelWorld::~G4VUserParallelWorld()
{ ; }

void G4VUserParallelWorld::ConstructSD()
{ ; }

G4VPhysicalVolume* G4VUserParallelWorld::GetWorld()
{
  G4VPhysicalVolume* pWorld
       = G4TransportationManager::GetTransportationManager()
         ->GetParallelWorld(fWorldName);
  pWorld->SetName(fWorldName);
  return pWorld;
}

void G4VUserParallelWorld::SetSensitiveDetector
(const G4String& logVolName, G4VSensitiveDetector* aSD, G4bool multi)
{
  G4bool found = false;
  G4LogicalVolumeStore* store = G4LogicalVolumeStore::GetInstance();
  for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++)
  {
    if((*pos)->GetName()==logVolName)
    {
      if(found && !multi)
      {
        G4String eM = "More than one logical volumes of the name <";
        eM += (*pos)->GetName();
        eM += "> are found and thus the sensitive detector <";
        eM += aSD->GetName();
        eM += "> cannot be uniquely assigned.";
        G4Exception("G4VUserParallelWorld::SetSensitiveDetector",
                    "Run5052",FatalErrorInArgument,eM);
      }
      found = true;
      SetSensitiveDetector(*pos,aSD);
    }
  }
  if(!found)
  {
    G4String eM2 = "No logical volume of the name <";
    eM2 += logVolName;
    eM2 += "> is found. The specified sensitive detector <";
    eM2 += aSD->GetName();
    eM2 += "> couldn't be assigned to any volume.";
    G4Exception("G4VUserParallelWorld::SetSensitiveDetector",
                    "Run5053",FatalErrorInArgument,eM2);
  }
}

void G4VUserParallelWorld::SetSensitiveDetector
(G4LogicalVolume* logVol, G4VSensitiveDetector* aSD)
{
  G4SDManager::GetSDMpointer()->AddNewDetector(aSD);
  logVol->SetSensitiveDetector(aSD);
}

