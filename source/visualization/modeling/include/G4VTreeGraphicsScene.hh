// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VTreeGraphicsScene.hh,v 1.1 2000-05-22 07:39:25 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  18th May 2000
// An artificial graphics scene to find dump geometry hierarchy.

#ifndef G4VTREEGRAPHICSSCENE_HH
#define G4VTREEGRAPHICSSCENE_HH

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

class G4VTreeGraphicsScene: public G4VGraphicsScene {

public:
  G4VTreeGraphicsScene();
  virtual ~G4VTreeGraphicsScene ();
  void AddThis (const G4Box& s) {Dump (s);}
  void AddThis (const G4Cons & s) {Dump (s);}
  void AddThis (const G4Tubs& s) {Dump (s);}
  void AddThis (const G4Trd& s) {Dump (s);}
  void AddThis (const G4Trap& s) {Dump (s);}
  void AddThis (const G4Sphere& s) {Dump (s);}
  void AddThis (const G4Para& s) {Dump (s);}
  void AddThis (const G4Torus& s) {Dump (s);}
  void AddThis (const G4Polycone& s) {Dump (s);}
  void AddThis (const G4Polyhedra& s) {Dump (s);}
  void AddThis (const G4VSolid& s) {Dump (s);}
  void PreAddThis (const G4Transform3D& objectTransformation,
		   const G4VisAttributes&);
  void PostAddThis ();
  void EstablishSpecials (G4PhysicalVolumeModel&);
  G4int                GetFoundDepth          () const;
  G4VPhysicalVolume*   GetFoundVolume         () const;
  const G4Transform3D& GetFoundTransformation () const;

  ////////////////////////////////////////////////////////////////
  // Functions not used but required by the abstract interface.

  virtual void BeginPrimitives (const G4Transform3D& objectTransformation) {}
  virtual void EndPrimitives () {}
  virtual void AddPrimitive (const G4Polyline&)   {}
  virtual void AddPrimitive (const G4Text&)       {}
  virtual void AddPrimitive (const G4Circle&)     {}
  virtual void AddPrimitive (const G4Square&)     {}
  virtual void AddPrimitive (const G4Polymarker&) {}
  virtual void AddPrimitive (const G4Polyhedron&) {}
  virtual void AddPrimitive (const G4NURBS&)      {}

protected:
  void Dump (const G4VSolid&);
  G4int                fCurrentDepth;  // Current depth of geom. hierarchy.
  G4VPhysicalVolume*   fpCurrentPV;    // Current physical volume.
  G4LogicalVolume*     fpCurrentLV;    // Current logical volume.
  const G4Transform3D* fpCurrentObjectTransformation;
};

#include "G4VTreeGraphicsScene.icc"

#endif
