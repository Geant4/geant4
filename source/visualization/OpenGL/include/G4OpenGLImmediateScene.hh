// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLImmediateScene.hh,v 1.1 1999-01-07 16:14:47 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  10th February 1997
// G4OpenGLImmediateScene - no Display lists.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLIMMEDIATESCENE_HH
#define G4OPENGLIMMEDIATESCENE_HH

#include "G4VScene.hh"
#include "G4OpenGLView.hh"
#include "G4OpenGLImmediateView.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"
#include "G4OpenGLScene.hh"

class G4OpenGLImmediate;

class G4OpenGLImmediateScene: public G4OpenGLScene {

public:
  G4OpenGLImmediateScene (G4VGraphicsSystem& system, const G4String& name);
  ~G4OpenGLImmediateScene ();
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
