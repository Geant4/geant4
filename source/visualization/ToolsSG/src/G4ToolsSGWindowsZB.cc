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
// Guy Barrand 18th April 2023


#include "G4ToolsSGWindowsZB.hh"

#include "G4ToolsSGViewer.hh"

#include <toolx/Windows/zb_viewer>

G4ToolsSGWindowsZB::G4ToolsSGWindowsZB():
parent
("TOOLSSG_WINDOWS_ZB",
 "TSG_WINDOWS_ZB",
 "TOOLSSG_WINDOWS_ZB is a graphics driver based on the g4tools tools/sg scene graph logic where\n\
 the rendering is done with the g4tools zbuffer and the windowing is done with Microsoft Windows.",
 parent::threeDInteractive)
,fSGSession(nullptr)
{}

G4ToolsSGWindowsZB::~G4ToolsSGWindowsZB() {
  delete fSGSession;
}

void G4ToolsSGWindowsZB::Initialise() {
  if(fSGSession) return; //done.
  fSGSession = new toolx::Windows::session(G4cout);
  if(!fSGSession->is_valid()) {
    G4cerr << "G4ToolsSGWindowsZB::Initialise : session::is_valid() failed." << G4endl;
    delete fSGSession;
    fSGSession = nullptr;
    return;
  }
}

G4VSceneHandler* G4ToolsSGWindowsZB::CreateSceneHandler(const G4String& a_name) {
  G4VSceneHandler* pScene = new G4ToolsSGSceneHandler(*this, a_name);
  return pScene;
}

G4VViewer* G4ToolsSGWindowsZB::CreateViewer(G4VSceneHandler& a_scene,const G4String& a_name) {
  if(!fSGSession) Initialise();
  if(!fSGSession) return nullptr;
  G4VViewer* pView =
    new G4ToolsSGViewer<toolx::Windows::session,toolx::Windows::zb_viewer>(*fSGSession,(G4ToolsSGSceneHandler&)a_scene,a_name);
  if (pView) {
    if (pView->GetViewId() < 0) {
      G4cerr <<
      "G4ToolsSGWindowsZB::CreateViewer: ERROR flagged by negative"
      " view id in G4ToolsSGViewer creation."
      "\n Destroying view and returning null pointer."
      << G4endl;
      delete pView;
      pView = nullptr;
    }
  }
  if (!pView) {
    G4cerr <<
    "G4ToolsSGWindowsZB::CreateViewer: ERROR: null pointer on new G4ToolsSGViewer."
    << G4endl;
  }
  return pView;
}
