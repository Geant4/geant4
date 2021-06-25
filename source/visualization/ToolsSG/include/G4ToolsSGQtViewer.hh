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
// Guy Barrand 8th April 2021

#if defined (G4VIS_BUILD_TOOLSSG_QT_GLES_DRIVER) || defined (G4VIS_USE_TOOLSSG_QT_GLES)

#ifndef G4TOOLSSGQTVIEWER_HH
#define G4TOOLSSGQTVIEWER_HH

#include "G4ToolsSGViewer.hh"

#include <tools/Qt/sg_viewer>

#include "G4Qt.hh"
#include "G4UIQt.hh"
#include "G4UImanager.hh"

class G4ToolsSGQtViewer : public G4ToolsSGViewer<tools::Qt::session,tools::Qt::sg_viewer> {
  typedef G4ToolsSGViewer<tools::Qt::session,tools::Qt::sg_viewer> parent;
public:
  G4ToolsSGQtViewer(tools::Qt::session& a_session,G4ToolsSGSceneHandler& a_scene_handler, const G4String& a_name)
  :parent(a_session,a_scene_handler,a_name)
  ,fInTab(false)
  ,fSGQWidget(nullptr)
  {}
  virtual ~G4ToolsSGQtViewer() = default;
protected:
  G4ToolsSGQtViewer(const G4ToolsSGQtViewer& a_from):parent(a_from),fInTab(false){}
  G4ToolsSGQtViewer& operator=(const G4ToolsSGQtViewer&) {return *this;}
public:  
  virtual void Initialise() {
    if(fSGQWidget) return; //done.
    parent::Initialise();
    fInTab = false;
    if(!fSGViewer) {
      G4cerr << "G4ToolsSGQtViewer::Initialise: ERROR: G4ToolsSGQtViewer has no tools::Qt::sg_viewer." << G4endl;
      return;
    }
    fSGQWidget = fSGViewer->shell();
    if (!fSGQWidget) {
      G4cerr << "G4ToolsSGQtViewer::Initialise: ERROR: tools::Qt::sg_viewer has no QWidget shell." << G4endl;
      return;
    }
  
   {G4UImanager* ui = G4UImanager::GetUIpointer();
    G4UIsession* session = ui->GetG4UIWindow();
    G4UIQt* uiQt = session? dynamic_cast<G4UIQt*>(session) :nullptr;
    if(uiQt) {      
      G4Qt* interactorManager = G4Qt::getInstance ();
      if (!interactorManager->IsExternalApp()) {
        uiQt->AddTabWidget(fSGQWidget,QString(GetName()));
        fInTab = true;
        fSGViewer->set_own_shell(false);
      }
    }}

    fSGViewer->enable_keyboard_focus();
  }
  virtual void ShowView() {
    if(fInTab) {
      if(!fSGViewer) return;
      fSGViewer->win_render();
      fSGSession.sync();
    } else {
      parent::ShowView();
    }
  }

  virtual void FinishView() {
    if(fInTab) {
      if(!fSGViewer) return;
      fSGViewer->win_render();
      fSGSession.sync();
    } else {
      parent::FinishView();
    }
  }
protected:
  bool fInTab;
  QWidget* fSGQWidget;
};

#endif

#endif
