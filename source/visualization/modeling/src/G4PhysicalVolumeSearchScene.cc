// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicalVolumeSearchScene.cc,v 1.3 1999-02-07 17:23:22 johna Exp $
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
	 << "\", copy no. " << fRequiredCopyNo << endl;
  G4cout << "PhysicalVolume:  \"" << fpCurrentPV -> GetName ()
	 << "\", copy no. " << fpCurrentPV -> GetCopyNo () << endl;
  *******************************************/

  if (fRequiredPhysicalVolumeName == "list") {
    for (G4int i = 0; i < fCurrentDepth; i++ ) G4cout << "  ";
    G4cout << "\"" << fpCurrentPV -> GetName ()
	   << "\", copy no. " << fpCurrentPV -> GetCopyNo ()
	   << endl;
    return;
  }

  if (fRequiredPhysicalVolumeName == fpCurrentPV -> GetName () &&
      fRequiredCopyNo             == fpCurrentPV -> GetCopyNo ()) {
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
	       << "\n  Required volume: \"" << fRequiredPhysicalVolumeName
	       << "\", copy no. " << fRequiredCopyNo
	       << " found more than once."
	  "\n  This function is not smart enough to distinguish identical"
	  "\n  physical volumes which have different parentage.  Besides,"
	  "\n  it's tricky to specify in general; we plan a GUI to do this."
	  "\n  This function gives you access to the first occurrence only."
	       << endl;
      }
    }
  }      
}
