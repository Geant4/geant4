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
#include "G4ScoringRealWorld.hh"

#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4Region.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"

#include "G4ScoringManager.hh"
#include "G4StatDouble.hh"

#include "G4SystemOfUnits.hh"

G4ScoringRealWorld::G4ScoringRealWorld(G4String lvName)
  : G4VScoringMesh(lvName)
{
  fShape          = MeshShape::realWorldLogVol;
  logVolName      = lvName;
  G4double size[] = { 0., 0., 0. };
  SetSize(size);
  G4int nBin[] = { 1, 1, 1 };
  SetNumberOfSegments(nBin);
}

void G4ScoringRealWorld::List() const
{
  G4cout << "G4ScoringRealWorld : " << logVolName << G4endl;
  G4VScoringMesh::List();
}

#include "G4AutoLock.hh"
namespace
{
  G4Mutex logvolmutex = G4MUTEX_INITIALIZER;
}
void G4ScoringRealWorld::SetupGeometry(G4VPhysicalVolume*)
{
  G4AutoLock l(&logvolmutex);
  auto store = G4LogicalVolumeStore::GetInstance();
  auto itr   = store->begin();
  for(; itr != store->end(); itr++)
  {
    if((*itr)->GetName() == logVolName)
    {
      fMeshElementLogical = (*itr);  // Logical volume to score
      G4int nb            = 0;
      auto pvStore        = G4PhysicalVolumeStore::GetInstance();
      auto pvItr          = pvStore->begin();
      for(; pvItr != pvStore->end(); pvItr++)
      {
        if((*pvItr)->GetLogicalVolume() == (*itr))
        {
          nb += (*pvItr)->GetMultiplicity();
        }
      }
      G4int nBin[] = { nb, 1, 1 };
      SetNumberOfSegments(nBin);
      // check if this logical volume belongs to the real world
      auto region = (*itr)->GetRegion();
      if((region != nullptr) && !(region->IsInMassGeometry()))
      {
        G4ExceptionDescription ed;
        ed << "Logical Volume with name <" << logVolName
           << "> is not used in the mass world.";
        G4Exception("G4ScoringRealWorld", "SWV0001", FatalException, ed);
      }
      // set the sensitive detector
      fMeshElementLogical->SetSensitiveDetector(fMFD);
      return;
    }
  }
  G4ExceptionDescription ed;
  ed << "Logical Volume with name <" << logVolName << "> is not found";
  G4Exception("G4ScoringRealWorld", "SWV0000", FatalException, ed);
}
