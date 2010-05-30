//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4OpenGLSceneHandler.hh,v 1.30 2010-05-30 09:53:05 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  27th March 1996
// OpenGL scene handler - base for immediate mode and stored mode classes to
//                        inherit from.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLSCENEHANDLER_HH
#define G4OPENGLSCENEHANDLER_HH

#include "globals.hh"

#include "G4VSceneHandler.hh"
#include "G4OpenGLBitMapStore.hh"
#include <GL/gl.h>
#include <GL/glu.h>

#include <map>

class G4AttHolder;

// Base class for various OpenGLSceneHandler classes.
class G4OpenGLSceneHandler: public G4VSceneHandler {

  friend class G4OpenGLViewer;

public:
  virtual void BeginPrimitives (const G4Transform3D& objectTransformation);
  virtual void EndPrimitives ();
  virtual void BeginPrimitives2D (const G4Transform3D& objectTransformation);
  virtual void EndPrimitives2D ();
  void AddPrimitive (const G4Polyline&);
  void AddPrimitive (const G4Polymarker&);
  void AddPrimitive (const G4Text&);
  void AddPrimitive (const G4Circle&);
  void AddPrimitive (const G4Square&);
  void AddPrimitives (std::vector <G4VMarker>);
  void AddPrimitive (const G4Scale&);
  void AddPrimitive (const G4Polyhedron&);
  void AddPrimitive (const G4NURBS&);

  void PreAddSolid (const G4Transform3D& objectTransformation,
		    const G4VisAttributes&);
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
  void AddCompound (const G4VDigi&);
  void AddCompound (const G4THitsMap<G4double>&);

protected:

  G4OpenGLSceneHandler (G4VGraphicsSystem& system,
			G4int id,
			const G4String& name = "");
  virtual ~G4OpenGLSceneHandler ();

  G4VSolid* CreateSectionSolid ();
  G4VSolid* CreateCutawaySolid ();

  void ClearAndDestroyAtts();  // Destroys att holders and clears pick map.

  GLuint fPickName;
  std::map<GLuint, G4AttHolder*> fPickMap;  // For picking.

  G4bool fProcessing2D;

private:

  void AddCircleSquare (const G4VMarker&, G4OpenGLBitMapStore::Shape);
  void AddCircleSquareVector (std::vector <G4VMarker>, G4OpenGLBitMapStore::Shape);

  void DrawXYPolygon
  (G4OpenGLBitMapStore::Shape,
   G4double size,
   const G4Point3D& centre,
   const G4VisAttributes* pApplicableVisAtts);
  // Draws in world coordinates a polygon in the screen plane knowing
  // viewpoint direction and up vector.

  static const GLubyte fStippleMaskHashed [128];

  G4bool fProcessingPolymarker;
};

#include "G4OpenGLSceneHandler.icc"

#endif

#endif
