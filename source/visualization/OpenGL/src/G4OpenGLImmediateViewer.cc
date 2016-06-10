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
// $Id: G4OpenGLImmediateViewer.cc 66373 2012-12-18 09:41:34Z gcosmo $
//
// 
// Andrew Walkden  7th February 1997
// Class G4OpenGLImmediateViewer : Encapsulates the `immediateness' of
//                               an OpenGL view, for inheritance by
//                               derived (X, Xm...) classes.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#include "G4OpenGLImmediateViewer.hh"
#include "G4OpenGLImmediateSceneHandler.hh"

G4OpenGLImmediateViewer::G4OpenGLImmediateViewer (G4OpenGLImmediateSceneHandler& scene):
G4VViewer (scene, -1),
G4OpenGLViewer (scene)
{}

void G4OpenGLImmediateViewer::ProcessView ()
{
  const G4Planes& cutaways = fVP.GetCutawayPlanes();
  G4bool cutawayUnion = fVP.IsCutaway() &&
    fVP.GetCutawayMode() == G4ViewParameters::cutawayUnion;
  size_t nPasses = cutawayUnion? cutaways.size(): 1;
  for (size_t i = 0; i < nPasses; ++i) {

    if (cutawayUnion) {
      double a[4];
      a[0] = cutaways[i].a();
      a[1] = cutaways[i].b();
      a[2] = cutaways[i].c();
      a[3] = cutaways[i].d();
      glClipPlane (GL_CLIP_PLANE2, a);
      glEnable (GL_CLIP_PLANE2);
    }

    NeedKernelVisit ();  // Always need to visit G4 kernel.
    G4VViewer::ProcessView ();

    if (cutawayUnion) glDisable (GL_CLIP_PLANE2);
  }
}

#endif
