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
// $Id: G4OpenGLWin32Viewer.cc,v 1.5 2002-10-16 10:44:16 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// G4OpenGLWin32Viewer : Class to provide WindowsNT specific
//                     functionality for OpenGL in GEANT4

#ifdef G4VIS_BUILD_OPENGLWIN32_DRIVER

#include "G4OpenGLWin32Viewer.hh"

#include <GL/gl.h>
//#include <GL/glx.h>  Change this for Win32

#include "G4ios.hh"

#include "G4VisExtent.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Point3D.hh"
#include "G4Normal3D.hh"

void G4OpenGLWin32Viewer::FinishView () {

  //Wait for GL commands to exit before going on...
  //Swap buffers if double buffered...
  //Flush GL commands...

}

void G4OpenGLWin32Viewer::GetWin32Connection () {
// get a connection.
}

void G4OpenGLWin32Viewer::CreateGLWin32Context () {
// create a GL context
// set attributes of window now the GL context has been created
    G4cout <<
      "ERROR: G4OpenGLWin32Viewer::CreateGLWin32Context:"
      "\n  **** WIN32 VIEWER NOT YET IMPLEMENTED ****"
           << G4endl;
}

void G4OpenGLWin32Viewer::CreateMainWindow () {
  
// create a window
// connect the context to a window

}

G4OpenGLWin32Viewer::G4OpenGLWin32Viewer (G4OpenGLSceneHandler& scene):
G4VViewer (scene, -1),
G4OpenGLViewer (scene)
{
  GetWin32Connection ();
  if (fViewId < 0) return;
  
  // Try for a visual suitable for OpenGLImmediate..
  // first try for a single buffered RGB window
  // next try for a double buffered RGB, but Draw to top buffer


  // Now try for a visual suitable for OpenGLStored...
  // Try for a double buffered RGB window
  
}

G4OpenGLWin32Viewer::~G4OpenGLWin32Viewer () {
  if (fViewId >= 0) {
    //Close a window from here
    //destroy GL context
    //make NULL context current
    //destroy the Win32 window
  }
}

#endif
