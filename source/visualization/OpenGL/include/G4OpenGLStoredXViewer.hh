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
// $Id: G4OpenGLStoredXViewer.hh,v 1.5 2001-07-14 21:47:46 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  7th February 1997
// Class G4OpenGLStoredXViewer : a class derived from G4OpenGLXViewer and
//                               G4OpenGLStoredViewer.

#ifdef G4VIS_BUILD_OPENGLX_DRIVER

#ifndef G4OPENGLSTOREDXVIEWER_HH
#define G4OPENGLSTOREDXVIEWER_HH

#include "G4VViewer.hh"
#include "G4OpenGLStoredViewer.hh"
#include "G4OpenGLXViewer.hh"

class G4OpenGLStoredSceneHandler;

class G4OpenGLStoredXViewer:
public G4OpenGLXViewer, public G4OpenGLStoredViewer{
  
public:
  G4OpenGLStoredXViewer (G4OpenGLStoredSceneHandler& scene, const G4String& name = "");
  virtual ~G4OpenGLStoredXViewer ();
  void Initialise ();
  void DrawView ();
};

#endif

#endif

