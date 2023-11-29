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
// John Allison  5th April 2001
// A base class for a scene handler to dump geometry hierarchy.
// Based on a provisional G4VTreeGraphicsScene (was in modeling).

#include "G4VTreeSceneHandler.hh"

#include "G4VSolid.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

#define G4warn G4cout

G4int G4VTreeSceneHandler::fSceneIdCount = 0;
// Counter for Tree scene handlers.

G4VTreeSceneHandler::G4VTreeSceneHandler(G4VGraphicsSystem& system,
					 const G4String& name):
  G4VSceneHandler(system, fSceneIdCount++, name),
  fpCurrentObjectTransformation (0)
{}

G4VTreeSceneHandler::~G4VTreeSceneHandler () {}

void G4VTreeSceneHandler::BeginModeling() {
  G4VSceneHandler::BeginModeling();  // Required: see G4VSceneHandler.hh.
}

void G4VTreeSceneHandler::EndModeling() {
  fDrawnLVStore.clear();
  G4VSceneHandler::EndModeling();  // Required: see G4VSceneHandler.hh.
}

void G4VTreeSceneHandler::PreAddSolid
(const G4Transform3D& objectTransformation, const G4VisAttributes& visAttribs)
{
  G4VSceneHandler::PreAddSolid (objectTransformation, visAttribs);

  G4PhysicalVolumeModel* pPVModel =
    dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
  if (!pPVModel) return;  // Not from a G4PhysicalVolumeModel.

  // This call comes from a G4PhysicalVolumeModel, drawnPVPath is
  // the path of the current drawn (non-culled) volume in terms of
  // drawn (non-culled) ancestors.  Each node is identified by a
  // PVNodeID object, which is a physical volume and copy number.  It
  // is a vector of PVNodeIDs corresponding to the geometry hierarchy
  // actually selected, i.e., not culled.
  typedef G4PhysicalVolumeModel::G4PhysicalVolumeNodeID PVNodeID;
  typedef std::vector<PVNodeID> PVPath;
  const PVPath& drawnPVPath = pPVModel->GetDrawnPVPath();
  //G4int currentDepth = pPVModel->GetCurrentDepth();
  //G4VPhysicalVolume* pCurrentPV = pPVModel->GetCurrentPV();
  //G4LogicalVolume* pCurrentLV = pPVModel->GetCurrentLV();
  //G4Material* pCurrentMaterial = pPVModel->GetCurrentMaterial();

  // Actually, it is enough to store the logical volume of current
  // physical volume...
  fDrawnLVStore.insert
    (drawnPVPath.back().GetPhysicalVolume()->GetLogicalVolume());

  // Find mother.  ri points to drawn mother, if any.
  PVPath::const_reverse_iterator ri = ++drawnPVPath.rbegin();
  if (ri != drawnPVPath.rend()) {
    // This volume has a mother.
    G4LogicalVolume* drawnMotherLV =
      ri->GetPhysicalVolume()->GetLogicalVolume();
    if (fDrawnLVStore.find(drawnMotherLV) != fDrawnLVStore.end()) {
      // Mother previously encountered.  Add this volume to
      // appropriate node in scene graph tree.
      // ...
    } else {
      // Mother not previously encountered.  Shouldn't happen, since
      // G4PhysicalVolumeModel sends volumes as it encounters them,
      // i.e., mothers before daughters, in its descent of the
      // geometry tree.  Error!
      G4warn << "ERROR: G4VTreeSceneHandler::PreAddSolid: Mother "
	     << ri->GetPhysicalVolume()->GetName()
	     << ':' << ri->GetCopyNo()
	     << " not previously encountered."
	"\nShouldn't happen!  Please report to visualization coordinator."
	     << G4endl;
      // Continue anyway.  Add to root of scene graph tree.
      // ...
    }
  } else {
    // This volume has no mother.  Must be a top level un-culled
    // volume.  Add to root of scene graph tree.
    // ...
  }
}
