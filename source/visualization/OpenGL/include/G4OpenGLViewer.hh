// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLViewer.hh,v 1.4 1999-11-11 15:38:04 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  27th March 1996
// OpenGL viewer - opens window, hard copy, etc.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLVIEWER_HH
#define G4OPENGLVIEWER_HH

#include "G4VViewer.hh"
#include "globals.hh"
#include "g4rw/tvordvec.h"

class G4OpenGLSceneHandler;

// Base class for various OpenGLView classes.
class G4OpenGLViewer: virtual public G4VViewer {

public:
  void ClearView  ();

protected:
  G4OpenGLViewer (G4OpenGLSceneHandler& scene);
  virtual ~G4OpenGLViewer ();
  void SetView    ();
  void HaloingFirstPass ();
  void HaloingSecondPass ();
  void HLRFirstPass ();
  void HLRSecondPass ();
  void HLRThirdPass ();
  void InitializeGLView ();
  G4bool white_background,  //the OpenGL clear colour
    doublebuffer,           //are we using a double buffered visual?
    transparency_enabled,   //is alpha blending enabled?
    antialiasing_enabled,   //is antialiasing enabled?
    haloing_enabled;        //is haloing enabled for wireframe?

  static int snglBuf_RGBA[10];
  static int dblBuf_RGBA[11];
  G4OpenGLSceneHandler&                    fSceneHandler;  // Graphics Scene for this view.

private:
  //  G4OpenGLSceneHandler&                    fSceneHandler;  // Graphics Scene for this view.
};

class G4OpenGLImmediateSceneHandler;
class G4OpenGLStoredSceneHandler;

#endif

#endif
