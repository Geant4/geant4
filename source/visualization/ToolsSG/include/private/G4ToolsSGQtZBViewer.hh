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
// Guy Barrand 13th April 2023

#ifndef G4TOOLSSGQTZBVIEWER_HH
#define G4TOOLSSGQTZBVIEWER_HH

#include "G4ToolsSGViewer.hh"

#include "G4Qt.hh"
#include "G4UIQt.hh"
#include "G4UImanager.hh"

#include <toolx/Qt/zb_viewer>

#include <qmainwindow.h>

class G4ToolsSGQtZBDestroyCallback : public QObject {
  Q_OBJECT
public:
  G4ToolsSGQtZBDestroyCallback(G4VViewer* a_viewer):fViewer(a_viewer) {}
  virtual ~G4ToolsSGQtZBDestroyCallback() {}
private:
  G4ToolsSGQtZBDestroyCallback(const G4ToolsSGQtZBDestroyCallback&);
  G4ToolsSGQtZBDestroyCallback& operator=(const G4ToolsSGQtZBDestroyCallback&);
public slots:
  void execute(){
    delete fViewer;
  }
private:
  G4VViewer* fViewer;
};

class G4ToolsSGQtZBViewer : public G4ToolsSGViewer<toolx::Qt::session,toolx::Qt::zb_viewer> {
  using parent = G4ToolsSGViewer<toolx::Qt::session,toolx::Qt::zb_viewer>;
public:
  G4ToolsSGQtZBViewer(toolx::Qt::session& a_session,G4ToolsSGSceneHandler& a_scene_handler, const G4String& a_name)
  :parent(a_session,a_scene_handler,a_name)
  ,fSGQWidget(nullptr)
  ,fDestroyCallback(0)
  {
    fDestroyCallback = new G4ToolsSGQtZBDestroyCallback(this);
  }
  virtual ~G4ToolsSGQtZBViewer() {
    delete fDestroyCallback; //it will remove the below signal/slot connection.
  }
protected:
  G4ToolsSGQtZBViewer(const G4ToolsSGQtZBViewer& a_from):parent(a_from){}
  G4ToolsSGQtZBViewer& operator=(const G4ToolsSGQtZBViewer&) {return *this;}
public:  
  virtual void Initialise() {
    if(fSGQWidget) return; //done.
    parent::Initialise();
    if(!fSGViewer) {
      G4cerr << "G4ToolsSGQtZBViewer::Initialise: ERROR: G4ToolsSGQtZBViewer has no toolx::Qt::zb_viewer." << G4endl;
      return;
    }
    fSGQWidget = fSGViewer->shell();
    if (!fSGQWidget) {
      G4cerr << "G4ToolsSGQtZBViewer::Initialise: ERROR: toolx::Qt::zb_viewer has no QWidget shell." << G4endl;
      return;
    }

   {G4UImanager* ui = G4UImanager::GetUIpointer();
    G4UIsession* session = ui->GetG4UIWindow();
     fUIQt = session? dynamic_cast<G4UIQt*>(session) :nullptr;
    if(fUIQt) {
      G4Qt* interactorManager = G4Qt::getInstance ();
      if (!interactorManager->IsExternalApp()) {
        fSGViewer->set_own_shell(false);
        fUIQt->AddTabWidget(fSGQWidget,QString(GetName()));
        QObject::connect(fSGQWidget,SIGNAL(destroyed()),fDestroyCallback,SLOT(execute()));

        QMainWindow* main_window = fUIQt->GetMainWindow();
        if(main_window) {
          main_window->show();
          interactorManager->FlushAndWaitExecution();
        }
      }
    }}

    fSGViewer->enable_keyboard_focus();
  }
  
  virtual void SetView() {
#ifdef __APPLE__
    if(fSGQWidget && fSGViewer) {
      if( (2*fSGQWidget->width() == int(fSGViewer->width()))   &&
          (2*fSGQWidget->height() == int(fSGViewer->height())) ){
        //  With Qt/Cocoa, the received size in
        // tools::sg::glarea::resizeGL is twice the QWidget::[width(),height()]!
        //  In general it does not pose problem, except when rendering 2D texts.
        // In order to have similar sizes than other platforms, we have to double
        // their pixel heights.
        fVP.SetGlobalMarkerScale(2);
      }
    }
#endif
    parent::SetView();
  }

  virtual void UpdateGUISceneTree() {
    if (fUIQt) fUIQt->UpdateSceneTree(fSceneTree);
  }

protected:
  G4UIQt* fUIQt = nullptr;
  QWidget* fSGQWidget;
  G4ToolsSGQtZBDestroyCallback* fDestroyCallback;
};

#endif
