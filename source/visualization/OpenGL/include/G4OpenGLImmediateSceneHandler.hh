// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLImmediateSceneHandler.hh,v 1.2 1999-05-10 14:03:39 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  10th February 1997
// G4OpenGLImmediateSceneHandler - no Display lists.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLIMMEDIATESCENEHANDLER_HH
#define G4OPENGLIMMEDIATESCENEHANDLER_HH

#include "G4VSceneHandler.hh"
#include "G4OpenGLViewer.hh"
#include "G4OpenGLImmediateViewer.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"
#include "G4OpenGLSceneHandler.hh"

class G4OpenGLImmediate;

class G4OpenGLImmediateSceneHandler: public G4OpenGLSceneHandler {

public:
  G4OpenGLImmediateSceneHandler (G4VGraphicsSystem& system, const G4String& name);
  virtual ~G4OpenGLImmediateSceneHandler ();
  void BeginPrimitives (const G4Transform3D& objectTransformation);
  void EndPrimitives ();
  void BeginModeling ();
  void EndModeling ();
  static G4int GetSceneCount ();

private:
  static G4int    fSceneIdCount;  // static counter for OpenGLImmediate scenes.
  static G4int    fSceneCount;    // No. of extanct scenes.
};

#endif

#endif
