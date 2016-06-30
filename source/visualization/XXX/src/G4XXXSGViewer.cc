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
// $Id: G4XXXSGViewer.cc 95225 2016-02-01 09:15:32Z gcosmo $
//
// 
// John Allison  10th March 2006
// A template for a sophisticated graphics driver with a scene graph.
//?? Lines beginning like this require specialisation for your driver.

#include "G4XXXSGViewer.hh"

#include "G4VSceneHandler.hh"
#include "G4XXXSGSceneHandler.hh"

#include <fstream>

G4XXXSGViewer::G4XXXSGViewer
(G4VSceneHandler& sceneHandler, const G4String& name):
  G4VViewer(sceneHandler, sceneHandler.IncrementViewCount(), name)
{}

G4XXXSGViewer::~G4XXXSGViewer() {}

void G4XXXSGViewer::SetView() {
  //#ifdef G4XXXSGDEBUG
  G4cout << "G4XXXSGViewer::SetView() called." << G4endl;
  //#endif
}

void G4XXXSGViewer::ClearView() {
  //#ifdef G4XXXSGDEBUG
  G4cout << "G4XXXSGViewer::ClearView() called." << G4endl;
  //#endif
}

void G4XXXSGViewer::DrawView() {
  //#ifdef G4XXXSGDEBUG
  G4cout << "G4XXXSGViewer::DrawView() called." << G4endl;
  //#endif

  // First, a view should decide when to re-visit the G4 kernel.
  // Sometimes it might not be necessary, e.g., if the scene is stored
  // in a graphical database (e.g., OpenGL's display lists) and only
  // the viewing angle has changed.  But graphics systems without a
  // graphical database will always need to visit the G4 kernel.

  // The fNeedKernelVisit flag might have been set by the user in
  // /vis/viewer/rebuild, but if not, make decision and set flag only
  // if necessary...
  if (!fNeedKernelVisit) KernelVisitDecision();
  G4bool kernelVisitWasNeeded = fNeedKernelVisit; // Keep (ProcessView resets).

  ProcessView ();  // Clears store and processes scene only if necessary.

  if (kernelVisitWasNeeded) {
    // Some systems, notably OpenGL, can draw while re-building, so
    // there might not be a need to draw from store again here.  But
    // in this case...
    DrawFromStore("G4XXXSGViewer::DrawView");
  } else {
    DrawFromStore("G4XXXSGViewer::DrawView");
  }

  // ...before finally...
  FinishView ();       // Flush streams and/or swap buffers.
}

void G4XXXSGViewer::ShowView() {
  //#ifdef G4XXXSGDEBUG
  G4cout << "G4XXXSGViewer::ShowView() called." << G4endl;
  //#endif
  // This is what you should see...
  DrawFromStore("G4XXXSGViewer::ShowView");
}

void G4XXXSGViewer::KernelVisitDecision () {
  
  // If there's a significant difference with the last view parameters
  // of either the scene handler or this viewer, trigger a rebuild.

  SceneGraph& sceneGraph =
    static_cast<G4XXXSGSceneHandler&>(fSceneHandler).fSceneGraph;
  if (sceneGraph.fDaughters.size() == 3  // I.e., only the root nodes.
      // (The above needs re-thinking.)
      || CompareForKernelVisit(fLastVP)) {
    NeedKernelVisit ();  // Sets fNeedKernelVisit.
  }      
  fLastVP = fVP;
}

G4bool G4XXXSGViewer::CompareForKernelVisit(G4ViewParameters& lastVP)
{
  // Typical comparison.  Taken from OpenGL.
  if (
      (lastVP.GetDrawingStyle ()    != fVP.GetDrawingStyle ())    ||
      (lastVP.IsAuxEdgeVisible ()   != fVP.IsAuxEdgeVisible ())   ||
      (lastVP.IsCulling ()          != fVP.IsCulling ())          ||
      (lastVP.IsCullingInvisible () != fVP.IsCullingInvisible ()) ||
      (lastVP.IsDensityCulling ()   != fVP.IsDensityCulling ())   ||
      (lastVP.IsCullingCovered ()   != fVP.IsCullingCovered ())   ||
      // No need to visit kernel if section plane changes.
      // No need to visit kernel if cutaway planes change.
      (lastVP.IsExplode ()          != fVP.IsExplode ())          ||
      (lastVP.GetNoOfSides ()       != fVP.GetNoOfSides ())       ||
      (lastVP.IsMarkerNotHidden ()  != fVP.IsMarkerNotHidden ())  ||
      (lastVP.GetDefaultVisAttributes()->GetColour() !=
       fVP.GetDefaultVisAttributes()->GetColour())                ||
      (lastVP.GetDefaultTextVisAttributes()->GetColour() !=
       fVP.GetDefaultTextVisAttributes()->GetColour())            ||
      (lastVP.GetBackgroundColour ()!= fVP.GetBackgroundColour ())||
      (lastVP.GetVisAttributesModifiers() !=
       fVP.GetVisAttributesModifiers())
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

void G4XXXSGViewer::DrawFromStore(const G4String& source) {
  SceneGraph& sceneGraph =
    static_cast<G4XXXSGSceneHandler&>(fSceneHandler).fSceneGraph;
  // Write to a file for testing...
  static G4int iCount = 0;
  std::ostringstream oss;
  oss << source << '.' << fName << '.' << iCount++ << ".out";
  G4cout << "Writing " << oss.str() << G4endl;
  std::ofstream ofs(oss.str().c_str());
  JA::PrintTree(ofs,&sceneGraph);
  ofs.close();
}
