// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLScene.hh,v 1.1 1999-01-07 16:14:49 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  27th March 1996
// OpenGL scene - base for immediate mode and stored mode classes to
//                inherit from.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLSCENE_HH
#define G4OPENGLSCENE_HH

#include <rw/tvhdict.h>

#include "G4VScene.hh"
#include "G4OpenGLView.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

// Base class for various OpenGLScene classes.
class G4OpenGLScene: public G4VScene {

public:
  void AddPrimitive (const G4Polyline&);
  void AddPrimitive (const G4Text&);
  void AddPrimitive (const G4Circle&);
  void AddPrimitive (const G4Square&);
  void AddPrimitive (const G4Polyhedron&);
  void AddPrimitive (const G4NURBS&);
  void AddPrimitive (const G4Polymarker&);

  void AddThis (const G4Box&);
  void AddThis (const G4Cons&);
  void AddThis (const G4Tubs&);
  void AddThis (const G4Trd&);
  void AddThis (const G4Trap&);
  void AddThis (const G4Sphere&);
  void AddThis (const G4Para&);
  void AddThis (const G4Torus&);
  void AddThis (const G4VSolid&);

protected:
  G4OpenGLScene (G4VGraphicsSystem& system,
		 G4int id,
		 const G4String& name = "");
  ~G4OpenGLScene ();
  G4bool initialize_hlr;

private:
  GLdouble clear_colour[4];
};

#include "G4OpenGLScene.icc"

#endif

#endif
