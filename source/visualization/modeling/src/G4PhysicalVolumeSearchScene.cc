//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4PhysicalVolumeSearchScene.cc,v 1.6 2001-07-11 10:09:23 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  10th August 1998.
// An artificial scene to find physical volumes.

#include "G4PhysicalVolumeSearchScene.hh"

#include "G4VSolid.hh"
#include "G4Vector3D.hh"
#include "G4PhysicalVolumeModel.hh"

G4PhysicalVolumeSearchScene::G4PhysicalVolumeSearchScene
(const G4String& requiredPhysicalVolumeName, G4int requiredCopyNo):
  fRequiredPhysicalVolumeName   (requiredPhysicalVolumeName),
  fRequiredCopyNo               (requiredCopyNo),
  fCurrentDepth                 (0),
  fpCurrentPV                   (0),
  fpCurrentLV                   (0),
  fpCurrentObjectTransformation (0),
  fFoundDepth                   (0),
  fpFoundPV                     (0),
  fpFoundLV                     (0),
  fMultipleOccurrence           (false)
{}

G4PhysicalVolumeSearchScene::~G4PhysicalVolumeSearchScene () {}

void G4PhysicalVolumeSearchScene::EstablishSpecials
(G4PhysicalVolumeModel& pvModel) {
  pvModel.DefinePointersToWorkingSpace (&fCurrentDepth,
					&fpCurrentPV,
					&fpCurrentLV);
}

void G4PhysicalVolumeSearchScene::FindVolume (const G4VSolid& solid) {

  /**************************************************
  G4cout << "Required volume: \"" << fRequiredPhysicalVolumeName
	 << "\", copy no. " << fRequiredCopyNo << G4endl;
  G4cout << "PhysicalVolume:  \"" << fpCurrentPV -> GetName ()
	 << "\", copy no. " << fpCurrentPV -> GetCopyNo () << G4endl;
  *******************************************/

  if (fRequiredPhysicalVolumeName == fpCurrentPV -> GetName () &&
      (fRequiredCopyNo             < 0 ||  // I.e., ignore negative request.
       fRequiredCopyNo             == fpCurrentPV -> GetCopyNo ())) {
    // Current policy - take first one found!!
    if (!fpFoundPV) {  // i.e., if not already found.
      fFoundDepth                = fCurrentDepth;
      fpFoundPV                  = fpCurrentPV;
      fpFoundLV                  = fpCurrentLV;
      fFoundObjectTransformation = *fpCurrentObjectTransformation;
    }
    else {
      if (!fMultipleOccurrence) {
	fMultipleOccurrence = true;
	G4cout << "G4PhysicalVolumeSearchScene::FindVolume:"
	       << "\n  Required volume: \""
	       << fRequiredPhysicalVolumeName
	       << "\"";
	if (fRequiredCopyNo >= 0) {
	  G4cout << ", copy no. " << fRequiredCopyNo << ",";
	}
	G4cout << " found more than once."
	  "\n  This function is not smart enough to distinguish identical"
	  "\n  physical volumes which have different parentage.  Besides,"
	  "\n  it's tricky to specify in general; we plan a GUI to do this."
	  "\n  This function gives you access to the first occurrence only."
	       << G4endl;
      }
    }
  }      
}
