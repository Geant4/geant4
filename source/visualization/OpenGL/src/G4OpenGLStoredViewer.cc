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
// $Id: G4OpenGLStoredViewer.cc,v 1.9 2002-11-27 12:44:10 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  7th February 1997
// Class G4OpenGLStoredViewer : Encapsulates the `storedness' of
//                            an OpenGL view, for inheritance by
//                            derived (X, Xm...) classes.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#include "G4OpenGLStoredViewer.hh"

#include <GL/gl.h>

#include "G4ios.hh"

#include "G4OpenGLStoredSceneHandler.hh"

G4OpenGLStoredViewer::G4OpenGLStoredViewer
(G4OpenGLStoredSceneHandler& sceneHandler):
G4VViewer (sceneHandler, -1),  
G4OpenGLViewer (sceneHandler), 
fG4OpenGLStoredSceneHandler (sceneHandler)
{
  fLastVP = fDefaultVP; // Not sure if this gets executed before or
  // after G4VViewer::G4VViewer!!  Doesn't matter much.
}

G4OpenGLStoredViewer::~G4OpenGLStoredViewer () {}

void G4OpenGLStoredViewer::KernelVisitDecision () {
  
  // If there's a significant difference with the last view parameters
  // of either the scene handler or this viewer, trigger a rebuild.

  if (!fG4OpenGLStoredSceneHandler.fTopPODL ||
      CompareForKernelVisit(fG4OpenGLStoredSceneHandler.fLastVP)  ||
      CompareForKernelVisit(fLastVP)) {
    NeedKernelVisit ();
  }      
  fLastVP = fVP;
  fG4OpenGLStoredSceneHandler.fLastVP = fVP;
}

G4bool G4OpenGLStoredViewer::CompareForKernelVisit(G4ViewParameters& lastVP) {

  if (
      (lastVP.GetDrawingStyle ()    != fVP.GetDrawingStyle ())    ||
      (lastVP.GetRepStyle ()        != fVP.GetRepStyle ())        ||
      (lastVP.IsCulling ()          != fVP.IsCulling ())          ||
      (lastVP.IsCullingInvisible () != fVP.IsCullingInvisible ()) ||
      (lastVP.IsDensityCulling ()   != fVP.IsDensityCulling ())   ||
      (lastVP.IsCullingCovered ()   != fVP.IsCullingCovered ())   ||
      (lastVP.IsSection ()          != fVP.IsSection ())          ||
      // No need to visit kernel if section plane changes.
      (lastVP.IsCutaway ()          != fVP.IsCutaway ())          ||
      (lastVP.GetCutawayPlanes ().size () !=
                                 fVP.GetCutawayPlanes ().size ()) ||
      // No need to visit kernel if cutaway planes change.
      (lastVP.IsExplode ()          != fVP.IsExplode ())          ||
      (lastVP.GetNoOfSides ()       != fVP.GetNoOfSides ())
      ) {
    return true;
  }

  if (lastVP.IsDensityCulling () &&
      (lastVP.GetVisibleDensity () != fVP.GetVisibleDensity ()))
    return true;

  if (lastVP.IsExplode () &&
      (lastVP.GetExplodeFactor () != fVP.GetExplodeFactor ()))
    return true;

  return false;
}

void G4OpenGLStoredViewer::DrawDisplayLists () {
  
  if (fG4OpenGLStoredSceneHandler.fTopPODL)
    glCallList (fG4OpenGLStoredSceneHandler.fTopPODL);
  
  G4int nTODLs = fG4OpenGLStoredSceneHandler.fTODLList.size ();
  for (int i = 0; i < nTODLs; i++) {
    glPushMatrix();
    G4OpenGLTransform3D oglt
      (fG4OpenGLStoredSceneHandler.fTODLTransformList [i]);
    glMultMatrixd (oglt.GetGLMatrix ());
    glCallList (fG4OpenGLStoredSceneHandler.fTODLList[i]);
    glPopMatrix();
  }
}

#endif
