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
// $Id: G4OpenGLImmediateViewer.hh,v 1.7 2002-10-16 10:44:13 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  7th February 1997
// Class G4OpenGLImmediateViewer : Encapsulates the `immediateness' of
//                                 an OpenGL viewer, for inheritance by
//                                 derived (X, Xm...) classes.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLIMMEDIATEVIEWER_HH
#define G4OPENGLIMMEDIATEVIEWER_HH

#include "G4VViewer.hh"
#include "G4OpenGLViewer.hh"
#include "G4OpenGLSceneHandler.hh"
#include "G4OpenGLImmediateSceneHandler.hh"
#include "G4OpenGLTransform3D.hh"
#include "globals.hh"

class G4OpenGLSceneHandler;
class G4OpenGLImmediateSceneHandler;

class G4OpenGLImmediateViewer: virtual public G4OpenGLViewer {
  
public:
  G4OpenGLImmediateViewer (G4OpenGLImmediateSceneHandler& scene);
};

#endif

#endif
