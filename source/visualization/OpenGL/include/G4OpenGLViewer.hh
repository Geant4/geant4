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
// $Id: G4OpenGLViewer.hh,v 1.8 2002-02-24 01:47:51 johna Exp $
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

  G4OpenGLSceneHandler&                    fSceneHandler;  // Graphics Scene for this view.

private:
  //  G4OpenGLSceneHandler&                    fSceneHandler;  // Graphics Scene for this view.
};

class G4OpenGLImmediateSceneHandler;
class G4OpenGLStoredSceneHandler;

#endif

#endif
