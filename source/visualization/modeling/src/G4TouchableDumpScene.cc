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
// $Id: G4TouchableDumpScene.cc 66773 2013-01-12 14:48:08Z allison $
//
// 
// John Allison  15th May 2014.
// An artificial scene to dump touchable attributes.

#include "G4TouchableDumpScene.hh"

#include "G4VSolid.hh"
#include "G4Vector3D.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeModel.hh"

#include "G4AttValue.hh"
#include "G4AttCheck.hh"

G4TouchableDumpScene::G4TouchableDumpScene
(std::ostream& os,
 G4PhysicalVolumeModel* pPVModel,
 const G4ModelingParameters::PVNameCopyNoPath& requiredTouchable)
:fos(os)
,fpPVModel(pPVModel)
,fRequiredTouchable(requiredTouchable)
,fFound(false)
{}

G4TouchableDumpScene::~G4TouchableDumpScene () {}

void G4TouchableDumpScene::ProcessVolume (const G4VSolid& solid) {

  const std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>&
  fullPVPath = fpPVModel->GetFullPVPath();

  if (fRequiredTouchable.size() == fullPVPath.size()) {
    // OK - there's a size match.  Check it out.
//    G4cout << "Size match" << G4endl;
    G4ModelingParameters::PVNameCopyNoPathConstIterator iNameCopyNo;
    std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>::const_iterator
    iPVNodeId;
    for (iNameCopyNo = fRequiredTouchable.begin(), iPVNodeId = fullPVPath.begin();
         iNameCopyNo != fRequiredTouchable.end();
         ++iNameCopyNo, ++iPVNodeId) {
//              G4cout
//              << iNameCopyNo->GetName()
//              << ',' << iNameCopyNo->GetCopyNo()
//              << "; " << iPVNodeId->GetPhysicalVolume()->GetName()
//              << ','  << iPVNodeId->GetPhysicalVolume()->GetCopyNo()
//              << G4endl;
      if (!(
            iNameCopyNo->GetName() ==
            iPVNodeId->GetPhysicalVolume()->GetName() &&
            iNameCopyNo->GetCopyNo() ==
            iPVNodeId->GetPhysicalVolume()->GetCopyNo()
            )) {
        break;
      }
    }
    if (iNameCopyNo == fRequiredTouchable.end()) {
//      G4cout << "Match found" << G4endl;
      fFound = true;

      const std::map<G4String,G4AttDef>* attDefs = fpPVModel->GetAttDefs();
      std::vector<G4AttValue>* attValues = fpPVModel->CreateCurrentAttValues();
      fos << G4AttCheck(attValues, attDefs);
      delete attValues;

      G4Polyhedron* polyhedron = solid.GetPolyhedron();
      fos << "\nLocal polyhedron coordinates:\n" << *polyhedron;
      G4Transform3D* transform = fpPVModel->GetCurrentTransform();
      polyhedron->Transform(*transform);
      fos << "\nGlobal polyhedron coordinates:\n" << *polyhedron;

      fpPVModel->Abort();  // No need to look further.
    }
  }
}
