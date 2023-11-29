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
// Guy Barrand 12th March 2021

#include "G4ToolsSGXtGLES.hh"

#include "G4ToolsSGViewer.hh"
#include <toolx/Xt/sg_viewer>

#include "G4Xt.hh"

class session : public toolx::Xt::session {
  using parent = toolx::Xt::session;
public:
  session(std::ostream& a_out)
  :parent(a_out)
  {
    m_app_widget = (Widget)G4Xt::getInstance()->GetMainInteractor();
    m_app_context = ::XtWidgetToApplicationContext(m_app_widget);
    m_app_owner = false;
  }
  virtual ~session() {}
protected:
  session(const session& a_from):parent(a_from){}
  session& operator=(const session& a_from){parent::operator=(a_from);return *this;}
};

G4ToolsSGXtGLES::G4ToolsSGXtGLES():
parent
("TOOLSSG_XT_GLES",
 "TSG_XT_GLES",
 "TOOLSSG_XT_GLES is a graphics driver based on the g4tools tools/sg scene graph logic where\n\
 the rendering is done with GLES and the windowing is done with the Xt toolkit.",
 parent::threeDInteractive)
,fSGSession(nullptr)
{}

G4ToolsSGXtGLES::~G4ToolsSGXtGLES() {
  delete fSGSession;
}

void G4ToolsSGXtGLES::Initialise() {
  if(fSGSession) return; //done.
  fSGSession = new session(G4cout);
  if(!fSGSession->is_valid()) {
    G4cerr << "G4ToolsSGXtGLES::Initialise : session::is_valid() failed." << G4endl;
    delete fSGSession;
    fSGSession = nullptr;
    return;
  }
}

G4VSceneHandler* G4ToolsSGXtGLES::CreateSceneHandler(const G4String& a_name) {
  G4VSceneHandler* pScene = new G4ToolsSGSceneHandler(*this, a_name);
  return pScene;
}

G4VViewer* G4ToolsSGXtGLES::CreateViewer(G4VSceneHandler& a_scene,const G4String& a_name) {
  if(!fSGSession) Initialise();
  if(!fSGSession) return nullptr;
  G4VViewer* pView =
    new G4ToolsSGViewer<toolx::Xt::session,toolx::Xt::sg_viewer>(*fSGSession,(G4ToolsSGSceneHandler&)a_scene,a_name);
  if (pView) {
    if (pView->GetViewId() < 0) {
      G4cerr <<
      "G4ToolsSGXtGLES::CreateViewer: ERROR flagged by negative"
      " view id in G4ToolsSGViewer creation."
      "\n Destroying view and returning null pointer."
      << G4endl;
      delete pView;
      pView = nullptr;
    }
  }
  if (!pView) {
    G4cerr <<
    "G4ToolsSGXtGLES::CreateViewer: ERROR: null pointer on new G4ToolsSGViewer."
    << G4endl;
  }
  return pView;
}
