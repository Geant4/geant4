// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhysicalVolumeSearchScene.hh,v 1.4 1999-05-10 14:04:25 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  10th August 1998.
// An artificial scene to find physical volumes.

#ifndef G4PHYSICALVOLUMESEARCHSCENE_HH
#define G4PHYSICALVOLUMESEARCHSCENE_HH

#include "G4VGraphicsScene.hh"
#include "G4VisExtent.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Sphere.hh"
#include "G4Para.hh"
#include "G4Torus.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

class G4PhysicalVolumeSearchScene: public G4VGraphicsScene {

public:
  G4PhysicalVolumeSearchScene
  (const G4String& requiredPhysicalVolumeName, G4int requiredCopyNo);
  virtual ~G4PhysicalVolumeSearchScene ();
  void AddThis (const G4Box& s) {FindVolume (s);}
  void AddThis (const G4Cons & s) {FindVolume (s);}
  void AddThis (const G4Tubs& s) {FindVolume (s);}
  void AddThis (const G4Trd& s) {FindVolume (s);}
  void AddThis (const G4Trap& s) {FindVolume (s);}
  void AddThis (const G4Sphere& s) {FindVolume (s);}
  void AddThis (const G4Para& s) {FindVolume (s);}
  void AddThis (const G4Torus& s) {FindVolume (s);}
  void AddThis (const G4Polycone& s) {FindVolume (s);}
  void AddThis (const G4Polyhedra& s) {FindVolume (s);}
  void AddThis (const G4VSolid& s) {FindVolume (s);}
  void PreAddThis (const G4Transform3D& objectTransformation,
		   const G4VisAttributes&);
  void PostAddThis ();
  void EstablishSpecials (G4PhysicalVolumeModel&);
  G4int                GetFoundDepth          () const;
  G4VPhysicalVolume*   GetFoundVolume         () const;
  const G4Transform3D& GetFoundTransformation () const;

private:
  void FindVolume (const G4VSolid&);
  G4String             fRequiredPhysicalVolumeName;
  G4int                fRequiredCopyNo;
  G4int                fCurrentDepth;  // Current depth of geom. hierarchy.
  G4VPhysicalVolume*   fpCurrentPV;    // Current physical volume.
  G4LogicalVolume*     fpCurrentLV;    // Current logical volume.
  const G4Transform3D* fpCurrentObjectTransformation;
  G4int                fFoundDepth;                  // Found depth.
  G4VPhysicalVolume*   fpFoundPV;                    // Found physical volume.
  G4LogicalVolume*     fpFoundLV;                    // Found logical volume.
  G4Transform3D        fFoundObjectTransformation;   // Found transformation.
  G4bool               fMultipleOccurrence;
};

#include "G4PhysicalVolumeSearchScene.icc"

#endif
