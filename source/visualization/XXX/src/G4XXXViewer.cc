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
// $Id: G4XXXViewer.cc,v 1.6 2003/11/06 15:06:51 johna Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $

#include "G4XXXViewer.hh"

#include "G4ios.hh"
#include <strstream>

#include "G4VSceneHandler.hh"

G4XXXViewer::G4XXXViewer
(G4VSceneHandler& sceneHandler, const G4String& name):
  G4VViewer(sceneHandler, sceneHandler.IncrementViewCount(), name) {}

G4XXXViewer::~G4XXXViewer() {}

void G4XXXViewer::SetView() {
#ifdef G4XXXDEBUG
  G4cout << "G4XXXViewer::SetView() called." << G4endl;
#endif
}

void G4XXXViewer::ClearView() {
#ifdef G4XXXDEBUG
  G4cout << "G4XXXViewer::ClearView() called." << G4endl;
#endif
}

void G4XXXViewer::DrawView() {
#ifdef G4XXXDEBUG
  G4cout << "G4XXXViewer::DrawView() called." << G4endl;
#endif

  // First, a view should decide when to re-visit the G4 kernel.
  // Sometimes it might not be necessary, e.g., if the scene is stored
  // in a graphical database (e.g., OpenGL's display lists) and only
  // the viewing angle has changed.  But graphics systems without a
  // graphical database will always need to visit the G4 kernel.

  NeedKernelVisit ();  // Default is - always visit G4 kernel.
  // Note: this routine sets the fNeedKernelVisit flag of *all* the
  // views of the scene.

  ProcessView ();      // The basic logic is here.

  // Then a view may have more to do, e.g., display the graphical
  // database.  That code should come here before finally...

  FinishView ();       // Flush streams and/or swap buffers.
}
