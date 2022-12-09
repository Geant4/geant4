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
// Guy Barrand 09th June 2022

#include "G4ToolsSGOffscreen.hh"

#include "G4ToolsSGOffscreenViewer.hh"

#include <tools/offscreen/sg_viewer>

G4ToolsSGOffscreen::G4ToolsSGOffscreen():
parent
("TOOLSSG_OFFSCREEN",
 "TSG_OFFSCREEN",
 "TOOLSSG_OFFSCREEN is a graphics driver based on the g4tools tools/sg scene graph logic where\n\
 the rendering is done by using various offscreen library as tools/sg/zb, gl2ps, png, jpeg.",
 parent::threeDInteractive)
,fSGSession(nullptr)
{}

G4ToolsSGOffscreen::~G4ToolsSGOffscreen() {
  delete fSGSession;
}

void G4ToolsSGOffscreen::Initialise() {
  if(fSGSession) return; //done.
  fSGSession = new tools::offscreen::session(G4cout);
  if(!fSGSession->is_valid()) {
    G4cerr << "G4ToolsSGOffscreen::Initialise : session::is_valid() failed." << G4endl;
    delete fSGSession;
    fSGSession = nullptr;
    return;
  }
}

G4VSceneHandler* G4ToolsSGOffscreen::CreateSceneHandler(const G4String& a_name) {
  G4VSceneHandler* pScene = new G4ToolsSGSceneHandler(*this,a_name);
  return pScene;
}

G4VViewer* G4ToolsSGOffscreen::CreateViewer(G4VSceneHandler& a_scene,const G4String& a_name) {
  if(!fSGSession) Initialise();
  if(!fSGSession) return nullptr;
  G4VViewer* pView = new G4ToolsSGOffscreenViewer(*fSGSession,(G4ToolsSGSceneHandler&)a_scene,a_name);
  if (pView) {
    if (pView->GetViewId() < 0) {
      G4cerr <<
      "G4ToolsSGOffscreen::CreateViewer: ERROR flagged by negative"
      " view id in G4ToolsSGViewer creation."
      "\n Destroying view and returning null pointer."
      << G4endl;
      delete pView;
      pView = nullptr;
    }
  }
  if (!pView) {
    G4cerr <<
    "G4ToolsSGOffscreen::CreateViewer: ERROR: null pointer on new G4ToolsSGViewer."
    << G4endl;
  }
  return pView;
}

G4bool G4ToolsSGOffscreen::IsUISessionCompatible () const
{
  //G4bool isCompatible = true;
  //return isCompatible;
  return true;
}
