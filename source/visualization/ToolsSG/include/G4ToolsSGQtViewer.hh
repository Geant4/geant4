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

#ifndef G4TOOLSSGQTVIEWER_HH
#define G4TOOLSSGQTVIEWER_HH

#include "G4ToolsSGViewer.hh"

#include "G4Qt.hh"
#include "G4UIQt.hh"
#include "G4UImanager.hh"

#include <toolx/Qt/sg_viewer>

#include <qmainwindow.h>

class G4ToolsSGQtDestroyCallback : public QObject {
  Q_OBJECT
public:
  G4ToolsSGQtDestroyCallback(G4VViewer* a_viewer):fViewer(a_viewer) {}
  virtual ~G4ToolsSGQtDestroyCallback() {}
private:
  G4ToolsSGQtDestroyCallback(const G4ToolsSGQtDestroyCallback&);
  G4ToolsSGQtDestroyCallback& operator=(const G4ToolsSGQtDestroyCallback&);
public slots:
  void execute(){
    //::printf("debug : G4ToolsSGQtDestroyCallback::execute\n");
    delete fViewer;
  }
private:
  G4VViewer* fViewer;
};

class G4ToolsSGQtViewer : public G4ToolsSGViewer<toolx::Qt::session,toolx::Qt::sg_viewer> {
  typedef G4ToolsSGViewer<toolx::Qt::session,toolx::Qt::sg_viewer> parent;
public:
  G4ToolsSGQtViewer(toolx::Qt::session& a_session,G4ToolsSGSceneHandler& a_scene_handler, const G4String& a_name)
  :parent(a_session,a_scene_handler,a_name)
  ,fSGQWidget(nullptr)
  ,fDestroyCallback(0)
  {
    fDestroyCallback = new G4ToolsSGQtDestroyCallback(this);
  }
  virtual ~G4ToolsSGQtViewer() {
    delete fDestroyCallback; //it will remove the below signal/slot connection.
  }
protected:
  G4ToolsSGQtViewer(const G4ToolsSGQtViewer& a_from):parent(a_from){}
  G4ToolsSGQtViewer& operator=(const G4ToolsSGQtViewer&) {return *this;}
public:  
  virtual void Initialise() {
    if(fSGQWidget) return; //done.
    parent::Initialise();
    if(!fSGViewer) {
      G4cerr << "G4ToolsSGQtViewer::Initialise: ERROR: G4ToolsSGQtViewer has no toolx::Qt::sg_viewer." << G4endl;
      return;
    }
    fSGQWidget = fSGViewer->shell();
    if (!fSGQWidget) {
      G4cerr << "G4ToolsSGQtViewer::Initialise: ERROR: toolx::Qt::sg_viewer has no QWidget shell." << G4endl;
      return;
    }

   {G4UImanager* ui = G4UImanager::GetUIpointer();
    G4UIsession* session = ui->GetG4UIWindow();
    G4UIQt* uiQt = session? dynamic_cast<G4UIQt*>(session) :nullptr;
    if(uiQt) {      
      G4Qt* interactorManager = G4Qt::getInstance ();
      if (!interactorManager->IsExternalApp()) {
        fSGViewer->set_own_shell(false);
        uiQt->AddTabWidget(fSGQWidget,QString(GetName()));
        QObject::connect(fSGQWidget,SIGNAL(destroyed()),fDestroyCallback,SLOT(execute()));

        QMainWindow* main_window = uiQt->GetMainWindow();
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
  
protected:
  QWidget* fSGQWidget;
  G4ToolsSGQtDestroyCallback* fDestroyCallback;
};

#endif
