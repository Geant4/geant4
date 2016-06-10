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
// $Id: G4OpenGLStoredWin32Viewer.hh 66373 2012-12-18 09:41:34Z gcosmo $
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
  G4OpenGLStoredWin32Viewer (G4OpenGLStoredSceneHandler& scene,
			     const G4String& name = "");
  void Initialise ();
  void DrawView ();
  void FinishView ();
};

#endif

#endif

