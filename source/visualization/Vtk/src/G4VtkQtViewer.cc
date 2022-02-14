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

#include "G4VSceneHandler.hh"

#include "G4VtkQtViewer.hh"
#include "G4VtkQtSceneHandler.hh"

#include "G4UImanager.hh"
#include "G4UIQt.hh"
#include "G4Qt.hh"

#include <array>
#include "vtkVersion.h"
#include "vtkNew.h"
#include "vtkNamedColors.h"
#include "vtkCylinderSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkGenericOpenGLRenderWindow.h"

G4VtkQtViewer::G4VtkQtViewer (G4VSceneHandler& sceneHandler, const G4String& name) :
  G4VtkViewer(sceneHandler, name) {
}

G4VtkQtViewer::~G4VtkQtViewer() {}

void G4VtkQtViewer::Initialise()
{
  CreateMainWindow(this, QString(GetName()));

  // Specific GL render window and interactor for Qt
  _renderWindow          = vtkGenericOpenGLRenderWindow::New();
  renderWindowInteractor = vtkRenderWindowInteractor::New();

  _renderWindow->AddRenderer(renderer);
#if VTK_MAJOR_VERSION == 8
  this->SetRenderWindow(_renderWindow);
#else
  this->setRenderWindow(_renderWindow);
#endif

  // Set callback to match VTK parameters to Geant4
  geant4Callback->SetGeant4ViewParameters(&fVP);
  renderer->AddObserver(vtkCommand::EndEvent, geant4Callback);
}

void G4VtkQtViewer::CreateMainWindow(QVTKOpenGLNativeWidget *vtkWidget,
                                     const QString& name) {
  // G4Qt* interactorManager = G4Qt::getInstance ();
  G4UImanager* UI = G4UImanager::GetUIpointer();
  fUiQt = static_cast<G4UIQt*> (UI->GetG4UIWindow());
  fUiQt->AddTabWidget((QWidget*)vtkWidget,name);
}

void G4VtkQtViewer::FinishView()
{
  G4VtkViewer::FinishView();
  // force a widget repaint as there should already
  // be a rendered buffer when visualiser starts up
  // paintGL();
  // repaint();
}