// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VGraphicsScene.hh,v 1.2 1999-02-07 17:02:03 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  19th July 1996
// Abstract interface class for the concept of a graphics scene.
// It is a minimal scene for the GEANT4 kernel.
// See G4VScene for a fuller description.  G4VScene is the full abstract
// interface to graphics systems.

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

#include "G4Transform3D.hh"

// Scene Interface - begin snippet.
class G4VGraphicsScene {

public:

  ///////////////////////////////////////////////////////////////////
  // Functions for adding raw GEANT4 objects.  The concrete graphics
  // scene has the option of implementing its own model or asking the
  // solid to provide a G4Polyhedron or similar primitive - see, for
  // example, G4VScene in the Visualization Category.

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
  // Other functions.

  virtual void EstablishSpecials (G4PhysicalVolumeModel&) {}
  // Used to establish any special relationships between scene and
  // this particular type of model - non-pure, i.e., no requirement to
  // implement.  See G4PhysicalVolumeModel.hh for details.

  virtual void DecommissionSpecials (G4PhysicalVolumeModel&) {}
  // Used to reverse the effect of EstablishSpecials, if required.

  virtual void PreAddThis (const G4Transform3D& objectTransformation,
			   const G4VisAttributes& visAttribs) = 0;
  // objectTransformation is the transformation in the world
  // coordinate system of the object about to be added, and
  // visAttribs is its visualization attributes.

  virtual void PostAddThis () = 0;

};
// Scene Interface - end snippet.

#endif
