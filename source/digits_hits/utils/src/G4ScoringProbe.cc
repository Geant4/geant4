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
#include "G4ScoringProbe.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4Threading.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"

#include "G4ScoringManager.hh"
#include "G4StatDouble.hh"

#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"

G4ScoringProbe::G4ScoringProbe(G4String lvName, G4double half_size,
                               G4bool checkOverlap)
  : G4VScoringMesh(lvName)
  , chkOverlap(checkOverlap)
  , layeredMaterialName("none")
  , layeredMaterial(nullptr)
{
  fShape        = MeshShape::probe;
  logVolName    = lvName;
  probeSize     = half_size;
  G4double hs[] = { half_size, half_size, half_size };
  SetSize(hs);
  G4int nBin[] = { 1, 1, 1 };
  SetNumberOfSegments(nBin);
  regName = lvName + "_region";
  if(G4Threading::IsMasterThread())
  {
    new G4Region(regName);
  }
}

void G4ScoringProbe::List() const
{
  G4cout << "G4ScoringProbe : " << logVolName << G4endl;
  std::size_t np = posVec.size();
  for(std::size_t i = 0; i < np; ++i)
  {
    G4cout << " >> probe #" << i << " at " << posVec[i] << G4endl;
  }
  G4VScoringMesh::List();
}

#include "G4AutoLock.hh"
namespace
{
  G4Mutex logvolmutex = G4MUTEX_INITIALIZER;
}

void G4ScoringProbe::SetupGeometry(G4VPhysicalVolume* worldPhys)
{
  if(G4Threading::IsMasterThread())
  {
    auto worldLog = worldPhys->GetLogicalVolume();
    auto region   = G4RegionStore::GetInstance()->GetRegion(regName);
    assert(region != nullptr);
    region->AddRootLogicalVolume(worldLog);
    region->SetWorld(worldPhys);

    auto boxSolid =
      new G4Box(logVolName + "_solid", probeSize, probeSize, probeSize);
    fMeshElementLogical =
      new G4LogicalVolume(boxSolid, layeredMaterial, logVolName + "_log");

    std::size_t np = posVec.size();
    for(std::size_t i = 0; i < np; ++i)
    {
      new G4PVPlacement(nullptr, posVec[i], fMeshElementLogical, logVolName + "_phy",
                        worldLog, false, (G4int)i, chkOverlap);
    }

    auto  wisatt = new G4VisAttributes(G4Colour(.5, .5, .5));
    wisatt->SetVisibility(false);
    worldLog->SetVisAttributes(wisatt);
    auto  visatt = new G4VisAttributes(G4Colour(.5, .5, .5));
    visatt->SetVisibility(true);
    fMeshElementLogical->SetVisAttributes(visatt);
  }
  else
  {
    G4AutoLock l(&logvolmutex);
    fMeshElementLogical =
      G4LogicalVolumeStore::GetInstance()->GetVolume(logVolName, false);
    assert(fMeshElementLogical != nullptr);
    l.unlock();
  }

  fMeshElementLogical->SetSensitiveDetector(fMFD);
}

G4bool G4ScoringProbe::SetMaterial(G4String val)
{
  if(G4Threading::IsMasterThread())
  {
    if(val == "none")
    {
      layeredMaterialName = val;
      layeredMassFlg      = false;
      layeredMaterial     = nullptr;
    }
    else
    {
      auto mat = G4NistManager::Instance()->FindOrBuildMaterial(val);
      if(mat == nullptr)
      {
        return false;
      }
      layeredMaterialName = val;
      layeredMassFlg      = true;
      layeredMaterial     = mat;
    }
    auto region = G4RegionStore::GetInstance()->GetRegion(regName);
    assert(region != nullptr);
    region->UpdateMaterialList();
  }
  return true;
}
