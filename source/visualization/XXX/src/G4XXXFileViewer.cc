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
// $Id: G4XXXFileViewer.cc,v 1.1 2006-03-28 17:16:41 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  7th March 2006
// A template for a file-writing graphics driver.
//?? Lines beginning like this require specialisation for your driver.

#include "G4XXXFileViewer.hh"

#include "G4VSceneHandler.hh"
#include "G4XXXFileSceneHandler.hh"

G4XXXFileViewer::G4XXXFileViewer
(G4VSceneHandler& sceneHandler, const G4String& name):
  G4VViewer(sceneHandler, sceneHandler.IncrementViewCount(), name),
  fFileWriter(name + ".out")
{}

G4XXXFileViewer::~G4XXXFileViewer() {}

void G4XXXFileViewer::SetView() {
#ifdef G4XXXFileDEBUG
  G4cout << "G4XXXFileViewer::SetView() called." << G4endl;
#endif
}

void G4XXXFileViewer::ClearView() {
#ifdef G4XXXFileDEBUG
  G4cout << "G4XXXFileViewer::ClearView() called." << G4endl;
#endif
  fFileWriter.Rewind();
}

void G4XXXFileViewer::DrawView() {
#ifdef G4XXXFileDEBUG
  G4cout << "G4XXXFileViewer::DrawView() called." << G4endl;
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
  // database.  That code should come here...

  // ...before finally...
  FinishView ();       // Flush streams and/or swap buffers.
}

void G4XXXFileViewer::ShowView() {
#ifdef G4XXXFileDEBUG
  G4cout << "G4XXXFileViewer::ShowView() called." << G4endl;
#endif
  fFileWriter.Close();
}
