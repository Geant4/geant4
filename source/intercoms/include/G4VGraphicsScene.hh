// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VGraphicsScene.hh,v 1.5 2000-05-19 06:29:48 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// John Allison  19th July 1996
//
// Class Description:
// Abstract interface class for a graphics scene handler.
// It is a minimal scene handler for the GEANT4 kernel.
// See G4VSceneHandler for a fuller description.  G4VSceneHandler is
// the full abstract interface to graphics systems.

#ifndef G4VGRAPHICSSCENE_HH
#define G4VGRAPHICSSCENE_HH

class G4VisAttributes;
class G4VSolid;
class G4Box;
class G4Cons;
class G4Tubs;
class G4Trd;
class G4Trap;
class G4Sphere;
class G4Para;
class G4Torus;
class G4PhysicalVolumeModel;
class G4Polycone;
class G4Polyhedra;
class G4Polyline;
class G4Text;
class G4Circle;
class G4Square;
class G4Polymarker;
class G4Polyhedron;
class G4NURBS;

#include "G4Transform3D.hh"

class G4VGraphicsScene {

public: // With description

  ///////////////////////////////////////////////////////////////////
  // Functions for adding raw GEANT4 objects to the scene handler.
  // The concrete graphics scene handler has the option of
  // implementing its own model or asking the solid to provide a
  // G4Polyhedron or similar primitive - see, for example,
  // G4VSceneHandler in the Visualization Category.

  virtual void PreAddThis (const G4Transform3D& objectTransformation,
			   const G4VisAttributes& visAttribs) = 0;
  // objectTransformation is the transformation in the world
  // coordinate system of the object about to be added, and
  // visAttribs is its visualization attributes.

  virtual void PostAddThis () = 0;

  virtual void AddThis (const G4Box&)       = 0;
  virtual void AddThis (const G4Cons&)      = 0;
  virtual void AddThis (const G4Tubs&)      = 0;
  virtual void AddThis (const G4Trd&)       = 0;
  virtual void AddThis (const G4Trap&)      = 0;
  virtual void AddThis (const G4Sphere&)    = 0;
  virtual void AddThis (const G4Para&)      = 0;
  virtual void AddThis (const G4Torus&)     = 0;
  virtual void AddThis (const G4Polycone&)  = 0;
  virtual void AddThis (const G4Polyhedra&) = 0;
  virtual void AddThis (const G4VSolid&)    = 0;  // For solids not above.

  ///////////////////////////////////////////////////////////////////
  // Functions for adding graphics primitives to the scene handler.

  virtual void BeginPrimitives (const G4Transform3D& objectTransformation) = 0;
  // IMPORTANT: invoke this from your polymorphic versions, e.g.:
  // void MyXXXSceneHandler::BeginPrimitives
  // (const G4Transform3D& objectTransformation) {
  //   G4VSceneHandler::BeginPrimitives (objectTransformation);
  //   ...
  // }

  virtual void EndPrimitives () = 0;
  // IMPORTANT: invoke this from your polymorphic versions, e.g.:
  // void MyXXXSceneHandler::EndPrimitives () {
  //   ...
  //   G4VSceneHandler::EndPrimitives ();
  // }

  virtual void AddPrimitive (const G4Polyline&)   = 0;
  virtual void AddPrimitive (const G4Text&)       = 0;
  virtual void AddPrimitive (const G4Circle&)     = 0;
  virtual void AddPrimitive (const G4Square&)     = 0;
  virtual void AddPrimitive (const G4Polymarker&) = 0;
  virtual void AddPrimitive (const G4Polyhedron&) = 0;
  virtual void AddPrimitive (const G4NURBS&)      = 0;

  ///////////////////////////////////////////////////////////////////
  // Special functions for particular models.

  virtual void EstablishSpecials (G4PhysicalVolumeModel&) {}
  // Used to establish any special relationships between scene and
  // this particular type of model - non-pure, i.e., no requirement to
  // implement.  See G4PhysicalVolumeModel.hh for details.

  virtual void DecommissionSpecials (G4PhysicalVolumeModel&) {}
  // Used to reverse the effect of EstablishSpecials, if required.

};

#endif
