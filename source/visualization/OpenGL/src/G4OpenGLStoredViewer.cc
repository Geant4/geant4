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
// $Id: G4OpenGLStoredViewer.cc,v 1.14 2006/06/29 21:19:18 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// 
// Andrew Walkden  7th February 1997
// Class G4OpenGLStoredViewer : Encapsulates the `storedness' of
//                            an OpenGL view, for inheritance by
//                            derived (X, Xm...) classes.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#include "G4OpenGLStoredViewer.hh"

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
      CompareForKernelVisit(fLastVP)) {
    NeedKernelVisit ();
  }      
  fLastVP = fVP;
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
      (lastVP.GetNoOfSides ()       != fVP.GetNoOfSides ())       ||
      (lastVP.GetBackgroundColour ()!= fVP.GetBackgroundColour ())
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
