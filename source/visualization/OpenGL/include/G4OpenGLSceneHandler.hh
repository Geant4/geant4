// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLSceneHandler.hh,v 1.5 1999-12-16 17:25:05 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  27th March 1996
// OpenGL scene handler - base for immediate mode and stored mode classes to
//                        inherit from.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLSCENEHANDLER_HH
#define G4OPENGLSCENEHANDLER_HH

#include "G4VSceneHandler.hh"
#include "G4OpenGLViewer.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

// Base class for various OpenGLScene classes.
class G4OpenGLSceneHandler: public G4VSceneHandler {

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
  void AddThis (const G4Polycone&);
  void AddThis (const G4Polyhedra&);
  void AddThis (const G4VSolid&);

protected:
  G4OpenGLSceneHandler (G4VGraphicsSystem& system,
		 G4int id,
		 const G4String& name = "");
  virtual ~G4OpenGLSceneHandler ();
  G4bool initialize_hlr;

private:
  GLdouble clear_colour[4];
};

#include "G4OpenGLSceneHandler.icc"

#endif

#endif
