// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VVisManager.hh,v 1.1 1999-04-28 14:19:32 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Abstract interface for GEANT4 Visualization Manager.
// John Allison 19/Oct/1996.

// This is a "Singleton", i.e., only one instance of it may exist.
// This is ensured by making the constructor private.

// It has only one public access function which is used to obtain a pointer
// to the concrete G4VisManager, should it exist.
// G4VVisManager* pVVMan =  G4VVisManager::GetConcreteInstance ();
// points to the real (concrete) G4VisManager, if a view is available for
// drawing, otherwise is zero.  Thus all code must be protected,
// for example, by:
//   if (pVVMan) pVVMan -> Draw (polyline);

#ifndef G4VVISMANAGER_HH
#define G4VVISMANAGER_HH

#include "G4Transform3D.hh"
#include "G4ThreeVector.hh"       // Just a typedef Hep3Vector.
#include "G4RotationMatrix.hh"    // Just a typedef HepRotation.

class G4Polyline;
class G4Text;
class G4Circle;
class G4Square;
class G4Polymarker;
class G4Polyhedron;
class G4NURBS;
class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VisAttributes;

class G4VVisManager {

public:

  static G4VVisManager* GetConcreteInstance ();
  // Returns pointer to actual visualization manager if a view is
  // available for drawing, else returns null.  Always check value.

  virtual ~G4VVisManager () {}

  ///////////////////////////////////////////////////////////////////
  // Functions to Draw "transient" objects, useful for hits, digis, etc.

  //1 Note that the {\tt G4Transform3D} objects refer to the
  //1 transformation of the {\em object} being drawn.  However, for
  //1 some functions, there is a version which takes a {\tt
  //1 G4Translation} and a {\tt G4RotationMatrix}, the latter being a
  //1 {\em system} rotation, as in {\tt G4PVPlacement}.
  //1
  //1 Note also that where a {\tt G4VisAttributes} argument is
  //1 specified, it overrides any attributes belonging to the object
  //1 itself.  Otherwise, the visualization attributes are assumed to
  //1 be those belonging to the object being drawn (you can set its
  //1 attributes --- see Section \ref{ap:setting_attribs}).

  // VVisManager Interface - begin snippet.
  virtual void Draw (const G4Polyline&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity) = 0;

  virtual void Draw (const G4Text&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity) = 0;

  virtual void Draw (const G4Circle&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity) = 0;

  virtual void Draw (const G4Square&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity) = 0;

  virtual void Draw (const G4Polymarker&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity) = 0;

  virtual void Draw (const G4Polyhedron&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity) = 0;

  virtual void Draw (const G4NURBS&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity) = 0;

  virtual void Draw (const G4VSolid&, const G4VisAttributes&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity) = 0;

  virtual void Draw (const G4LogicalVolume&, const G4VisAttributes&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity) = 0;

  virtual void Draw (const G4VPhysicalVolume&, const G4VisAttributes&,
    const G4Transform3D& objectTransformation = G4Transform3D::Identity) = 0;
  // VVisManager Interface - end snippet.

  // Other management functions...

  virtual void GeometryHasChanged () = 0;
  // This is used by the run manager to notify a change of geometry.

protected:

  static G4VVisManager* fpConcreteInstance;  // Pointer to real G4VisManager.

};

inline G4VVisManager* G4VVisManager::GetConcreteInstance () {
  return fpConcreteInstance;
}

#endif
