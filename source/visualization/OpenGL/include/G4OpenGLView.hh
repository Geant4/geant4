// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLView.hh,v 1.1 1999-01-07 16:14:51 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  27th March 1996
// OpenGL view - opens window, hard copy, etc.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLVIEW_HH
#define G4OPENGLVIEW_HH

#include "G4VView.hh"
#include "globals.hh"
#include <rw/tvordvec.h>

class G4OpenGLScene;

// Base class for various OpenGLView classes.
class G4OpenGLView: virtual public G4VView {

public:
  void ClearView  ();

protected:
  G4OpenGLView (G4OpenGLScene& scene);
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
  G4OpenGLScene&                    fScene;  // Graphics Scene for this view.

private:
  //  G4OpenGLScene&                    fScene;  // Graphics Scene for this view.
};

class G4OpenGLImmediateScene;
class G4OpenGLStoredScene;

#endif

#endif
