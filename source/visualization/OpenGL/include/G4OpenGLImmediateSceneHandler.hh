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
// $Id: G4OpenGLImmediateSceneHandler.hh 80466 2014-04-23 12:07:25Z gcosmo $
//
// 
// Andrew Walkden  10th February 1997
// G4OpenGLImmediateSceneHandler - no Display lists.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLIMMEDIATESCENEHANDLER_HH
#define G4OPENGLIMMEDIATESCENEHANDLER_HH

#include "G4OpenGLSceneHandler.hh"
#include "globals.hh"

class G4OpenGLImmediate;

class G4OpenGLImmediateSceneHandler: public G4OpenGLSceneHandler {

public:

  G4OpenGLImmediateSceneHandler (G4VGraphicsSystem& system, const G4String& name);
  virtual ~G4OpenGLImmediateSceneHandler ();
  void BeginPrimitives (const G4Transform3D& objectTransformation);
  void EndPrimitives ();
  void BeginPrimitives2D (const G4Transform3D& objectTransformation);
  void EndPrimitives2D ();
  void BeginModeling ();
  void EndModeling ();
  void AddPrimitive (const G4Polyline&);
  void AddPrimitive (const G4Polymarker&);
  void AddPrimitive (const G4Text& text);
  void AddPrimitive (const G4Circle&);
  void AddPrimitive (const G4Square&);
  void AddPrimitive (const G4Scale& scale);
  void AddPrimitive (const G4Polyhedron&);

protected:

  G4bool AddPrimitivePreamble(const G4VMarker& visible);
  G4bool AddPrimitivePreamble(const G4Polyline& visible);
  G4bool AddPrimitivePreamble(const G4Polyhedron& visible);

  void ClearTransientStore ();

  static G4int    fSceneIdCount;  // static counter for OpenGLImmediate scenes.
  private:
  bool AddPrimitivePreambleInternal(const G4Visible& visible, bool isMarker, bool isPolyline);
};

#endif

#endif
