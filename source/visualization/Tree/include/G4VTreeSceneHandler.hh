// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VTreeSceneHandler.hh,v 1.2 2001-06-05 09:59:15 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A base class for a scene handler to dump geometry hierarchy.

#ifndef G4VTREESCENEHANDLER_HH
#define G4VTREESCENEHANDLER_HH

#include "G4VSceneHandler.hh"
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

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4ModelingParameters;

class G4VTreeSceneHandler: public G4VSceneHandler {

public:
  G4VTreeSceneHandler(G4VGraphicsSystem& system,
		      const G4String& name);
  virtual ~G4VTreeSceneHandler ();
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

  virtual void BeginModeling();
  virtual void EndModeling();

  static G4int GetSceneCount();

protected:
  virtual void         Dump (const G4VSolid&) = 0;
  static G4int         fSceneIdCount;  // Counter for Tree scene handlers.
  static G4int         fSceneCount;    // No. of extanct scene handlers.
  G4int                fCurrentDepth;  // Current depth of geom. hierarchy.
  G4VPhysicalVolume*   fpCurrentPV;    // Current physical volume.
  G4LogicalVolume*     fpCurrentLV;    // Current logical volume.
  const G4Transform3D* fpCurrentObjectTransformation;
  const G4ModelingParameters* fpOriginalMP;  // Keeps pointer to original.
  G4ModelingParameters* fpNonCullingMP;      // For temporary non-culling.
};

#include "G4VTreeSceneHandler.icc"

#endif
