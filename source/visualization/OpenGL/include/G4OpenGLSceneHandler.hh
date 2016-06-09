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
// $Id: G4OpenGLSceneHandler.hh,v 1.18 2005/01/27 20:04:52 johna Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// Andrew Walkden  27th March 1996
// OpenGL scene handler - base for immediate mode and stored mode classes to
//                        inherit from.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLSCENEHANDLER_HH
#define G4OPENGLSCENEHANDLER_HH

#include "globals.hh"
#include "G4RotationMatrix.hh"

#include "G4OpenGL.hh"

#include "G4VSceneHandler.hh"
#include "G4OpenGLViewer.hh"

// Base class for various OpenGLScene classes.
class G4OpenGLSceneHandler: public G4VSceneHandler {

public:
  void AddPrimitive (const G4Polyline&);
  void AddPrimitive (const G4Text&);
  void AddPrimitive (const G4Circle&);
  void AddPrimitive (const G4Square&);
  void AddPrimitive (const G4Polyhedron&);
  void AddPrimitive (const G4NURBS&);
  // Explicitly invoke base class methods to avoid warnings about
  // hiding of base class methods...
  void AddPrimitive(const G4Polymarker& polymarker) {
    G4VSceneHandler::AddPrimitive (polymarker);
  }
  void AddPrimitive (const G4Scale& scale) {
    G4VSceneHandler::AddPrimitive (scale);
  }

  void AddSolid (const G4Box&);
  void AddSolid (const G4Cons&);
  void AddSolid (const G4Tubs&);
  void AddSolid (const G4Trd&);
  void AddSolid (const G4Trap&);
  void AddSolid (const G4Sphere&);
  void AddSolid (const G4Para&);
  void AddSolid (const G4Torus&);
  void AddSolid (const G4Polycone&);
  void AddSolid (const G4Polyhedra&);
  void AddSolid (const G4VSolid&);
  void AddCompound (const G4VTrajectory&);
  void AddCompound (const G4VHit&);

protected:
  G4OpenGLSceneHandler (G4VGraphicsSystem& system,
			G4int id,
			const G4String& name = "");
  virtual ~G4OpenGLSceneHandler ();

private:
  void AddCircleSquare (const G4VMarker&, G4int nSides);
  /**************************************************
  Not needed - but see note on future development in .cc.
  void DrawScreenPolygon (G4double size,
		          const G4Point3D& centre,
		          G4int nSides);
  // Draws in screen coordinates.
  *********************************/

  void DrawXYPolygon (G4double size,
		      const G4Point3D& centre,
		      G4int nSides);
  // Draws in world coordinates a polygon in the screen plane knowing
  // viewpoint direction and up vector.

  static const GLubyte fStippleMaskHashed [128];
};

#include "G4OpenGLSceneHandler.icc"

#endif

#endif
