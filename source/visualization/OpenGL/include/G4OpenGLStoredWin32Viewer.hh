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
// $Id: G4OpenGLStoredWin32Viewer.hh,v 1.3 2001-07-11 10:08:50 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Class G4OpenGLStoredWin32Viewer : a class derived from
//   G4OpenGLWin32Viewer and G4OpenGLStoredViewer.

#ifdef G4VIS_BUILD_OPENGLWIN32_DRIVER

#ifndef G4OPENGLSTOREDWIN32VIEWER_HH
#define G4OPENGLSTOREDWIN32VIEWER_HH

#include "G4VViewer.hh"
#include "G4OpenGLStoredViewer.hh"
#include "G4OpenGLWin32Viewer.hh"

class G4OpenGLStoredSceneHandler;

class G4OpenGLStoredWin32Viewer:
public G4OpenGLWin32Viewer, public G4OpenGLStoredViewer{
  
public:
  G4OpenGLStoredWin32Viewer (G4OpenGLStoredSceneHandler& scene);
  void DrawView ();
};

#endif

#endif

