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
// $Id$
//
// 
// John Allison  10th August 1998.
// An artificial scene to find physical volumes.

#include "G4PhysicalVolumeSearchScene.hh"

#include "G4VSolid.hh"
#include "G4Vector3D.hh"
#include "G4PhysicalVolumeModel.hh"

G4PhysicalVolumeSearchScene::G4PhysicalVolumeSearchScene
(G4PhysicalVolumeModel* pPVModel,
 const G4String& requiredPhysicalVolumeName,
 G4int requiredCopyNo):
  fpPVModel                     (pPVModel),
  fRequiredPhysicalVolumeName   (requiredPhysicalVolumeName),
  fRequiredCopyNo               (requiredCopyNo),
  fpCurrentObjectTransformation (0),
  fFoundDepth                   (0),
  fpFoundPV                     (0),
  fpFoundLV                     (0),
  fMultipleOccurrence           (false)
{}

G4PhysicalVolumeSearchScene::~G4PhysicalVolumeSearchScene () {}

void G4PhysicalVolumeSearchScene::FindVolume (const G4VSolid&) {

  typedef G4PhysicalVolumeModel::G4PhysicalVolumeNodeID PVNodeID;
  typedef std::vector<PVNodeID> PVPath;
  const PVPath& fullPVPath = fpPVModel->GetFullPVPath();
  G4int currentDepth = fpPVModel->GetCurrentDepth();
  G4VPhysicalVolume* pCurrentPV = fpPVModel->GetCurrentPV();
  G4LogicalVolume* pCurrentLV = fpPVModel->GetCurrentLV();
  //G4Material* pCurrentMaterial = fpPVModel->GetCurrentMaterial();
  // Note: pCurrentMaterial may be zero (parallel world).

  /**************************************************
  G4cout << "Required volume: \"" << fRequiredPhysicalVolumeName
	 << "\", copy no. " << fRequiredCopyNo << G4endl;
  G4cout << "PhysicalVolume:  \"" << pCurrentPV -> GetName ()
	 << "\", copy no. " << pCurrentPV -> GetCopyNo () << G4endl;
  *******************************************/

  if (fRequiredPhysicalVolumeName == pCurrentPV -> GetName () &&
      (fRequiredCopyNo             < 0 ||  // I.e., ignore negative request.
       fRequiredCopyNo             == pCurrentPV -> GetCopyNo ())) {
    // Current policy - take first one found!!
    if (!fpFoundPV) {  // i.e., if not already found.
      fFoundFullPVPath           = fullPVPath;
      fFoundDepth                = currentDepth;
      fpFoundPV                  = pCurrentPV;
      fpFoundLV                  = pCurrentLV;
      fFoundObjectTransformation = *fpCurrentObjectTransformation;
    }
    else {
      if (!fMultipleOccurrence) {
	fMultipleOccurrence = true;
	G4cout << "G4PhysicalVolumeSearchScene::FindVolume:"
	       << "\n  Required volume \""
	       << fRequiredPhysicalVolumeName
	       << "\"";
	if (fRequiredCopyNo >= 0) {
	  G4cout << ", copy no. " << fRequiredCopyNo << ",";
	}
	G4cout << " found more than once."
	  "\n  This function is not smart enough to distinguish identical"
	  "\n  physical volumes which have different parentage.  It is"
	  "\n  tricky to specify in general.  This function gives you access"
	  "\n  to the first occurrence only."
	       << G4endl;
      }
    }
  }      
}
