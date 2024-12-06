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

#include "G4VtkQtViewer.hh"

#include "G4LogicalVolume.hh"
#include "G4Qt.hh"
#include "G4UIQt.hh"
#include "G4UImanager.hh"
#include "G4VSceneHandler.hh"
#include "G4VtkInteractorStyle.hh"
#include "G4VtkQtSceneHandler.hh"
#include "G4VtkUtility.hh"

#include <qcoreapplication.h>
#include <qlabel.h>
#include <qlayout.h>
#include <qlineedit.h>
#include <qopenglcontext.h>
#include <qpushbutton.h>
#include <qthread.h>
#include <qtreewidget.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCylinderSource.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkVersion.h>

#include <array>

G4VtkQtViewer::G4VtkQtViewer(G4VSceneHandler& sceneHandler, const G4String& name)
  : G4VtkViewer(sceneHandler, name)
{
  G4Qt::getInstance();
  //this->setFormat(QVTKOpenGLNativeWidget::defaultFormat());
}

G4VtkQtViewer::~G4VtkQtViewer()
{
}

void G4VtkQtViewer::Initialise()
{
  CreateMainWindow(this, QString(GetName()));

  // Specific GL render window and interactor for Qt
  _renderWindow = vtkGenericOpenGLRenderWindow::New();

  _renderWindow->AddRenderer(renderer);
  this->setRenderWindow(_renderWindow);

  // Set callback to match VTK parameters to Geant4
  geant4Callback->SetGeant4ViewParameters(&fVP);
  renderer->AddObserver(vtkCommand::EndEvent, geant4Callback);

  // Hidden line removal
  renderer->SetUseHiddenLineRemoval(0);

  // Shadows
  renderer->SetUseShadows(0);

  vtkSmartPointer<G4VtkInteractorStyle> style =
    vtkSmartPointer<G4VtkInteractorStyle>::New();
  this->interactor()->SetInteractorStyle(style);
}

void G4VtkQtViewer::CreateMainWindow(QVTKOpenGLNativeWidget* vtkWidget, const QString& name)
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  fUiQt = static_cast<G4UIQt*>(UI->GetG4UIWindow());
  fUiQt->AddTabWidget((QWidget*)vtkWidget, name);
  vtkWidget->setAttribute(Qt::WA_AcceptTouchEvents, false);
  fGLWidget = vtkWidget;
  createSceneTreeWidget();
}

#ifdef G4MULTITHREADED

namespace {
  G4Mutex visSubThreadMutex = G4MUTEX_INITIALIZER;
  G4Condition waitForVisSubThreadInitialized = G4CONDITION_INITIALIZER;
  G4bool visSubThreadEstablished = false;
  G4bool qObjectsSwitched = false;
  QVTKOpenGLNativeWidget* qVtkW = nullptr;
}

void G4VtkQtViewer::DoneWithMasterThread()
{
  // Called by Main Thread !

  // Initialise and check qVtkW - no need to check in subsequent functions
  qVtkW = dynamic_cast<QVTKOpenGLNativeWidget*>(fGLWidget);
  if (qVtkW == nullptr) return;

  // Done with master thread
  qVtkW->doneCurrent();

  // Set current QThread for the way back
  fQVtkContextMainThread = QThread::currentThread();
}

void G4VtkQtViewer::MovingToVisSubThread()
{
  // Still on master thread but vis thread has been launched

  if (qVtkW == nullptr) return;

  // Wait until SwitchToVisSubThread has found vis sub-thread QThread
  {
  G4AutoLock lock(&visSubThreadMutex);
  G4CONDITIONWAITLAMBDA(&waitForVisSubThreadInitialized, &lock, []{return visSubThreadEstablished;})
  }

  // Move stuff to sub-thread
  if(qVtkW->context()) qVtkW->context()->moveToThread(fQVtkContextVisSubThread);

  // Inform sub-thread
  G4AutoLock lock(&visSubThreadMutex);
  qObjectsSwitched = true;
  lock.unlock();
  G4CONDITIONBROADCAST(&waitForVisSubThreadInitialized);
}

void G4VtkQtViewer::SwitchToVisSubThread()
{
  // Called by VisSub Thread !

  if (qVtkW == nullptr) return;

  fQVtkContextVisSubThread = QThread::currentThread();

  // Let MovingToVisSubThread know we have the QThread
  {
  G4AutoLock lock(&visSubThreadMutex);
  visSubThreadEstablished = true;
  G4CONDITIONBROADCAST(&waitForVisSubThreadInitialized);
  }

  // Wait until MovingToVisSubThread has moved stuff
  {
  G4AutoLock lock(&visSubThreadMutex);
  G4CONDITIONWAITLAMBDA(&waitForVisSubThreadInitialized, &lock, []{return qObjectsSwitched;})
  }

  // make context current
  qVtkW->makeCurrent();
}

void G4VtkQtViewer::DoneWithVisSubThread()
{
  // Called by vis sub thread

  if (qVtkW == nullptr) return;

  // finish with this vis sub thread context
  qVtkW->doneCurrent();

  // and move it back to the main thread
  if(qVtkW->context()) qVtkW->context()->moveToThread(fQVtkContextMainThread);
}

void G4VtkQtViewer::SwitchToMasterThread()
{
  // Called by VisSub Thread !

  if (qVtkW == nullptr) return;

  qVtkW->makeCurrent();

  visSubThreadEstablished = false;
  qObjectsSwitched = false;
}

#endif

void G4VtkQtViewer::FinishView()
{
  auto& fVtkSceneHandler = dynamic_cast<G4VtkSceneHandler&>(fSceneHandler);
  fVtkSceneHandler.Modified();

  _renderWindow->Render();

  auto qGLW = dynamic_cast<QVTKOpenGLNativeWidget*>(fGLWidget);
  qGLW->interactor()->Initialize();
  qGLW->interactor()->Start();
}

void G4VtkQtViewer::createSceneTreeWidget() {}

void G4VtkQtViewer::createSceneTreeComponent() {}

QTreeWidgetItem*
G4VtkQtViewer::createTreeWidgetItem(const PVPath& /*fullPath*/, const QString& /*name*/,
                                    int /*copyNb*/, int /*POIndex*/, const QString& /*logicalName*/,
                                    Qt::CheckState /*state*/, QTreeWidgetItem* /*parentTreeNode*/,
                                    const G4Colour& /*color*/)
{
  QTreeWidgetItem* newItem = nullptr;
  return newItem;
}

void G4VtkQtViewer::addNonPVSceneTreeElement(const G4String& /*model*/, G4Visible& /*visible*/,
                                             int /*currentPOIndex*/)
{}

void G4VtkQtViewer::addPVSceneTreeElement(const G4String& /*model*/,
                                          G4PhysicalVolumeModel* /*pPVModel*/,
                                          int /*currentPOIndex*/)
{}

QString G4VtkQtViewer::getModelShortName(const G4String& /*model*/)
{
  QString modelShortName;
  return modelShortName;
}

bool G4VtkQtViewer::parseAndInsertInSceneTree(QTreeWidgetItem* /*parentItem*/,
                                              G4PhysicalVolumeModel* /*pPVModel*/,
                                              unsigned int /*fullPathIndex*/,
                                              const QString& /*parentRoot*/,
                                              unsigned int /*currentIndexInTreeSceneHandler*/,
                                              int /*currentPVPOIndex*/)
{
  return false;
}

void G4VtkQtViewer::EnableClipperWidget()
{
  G4VtkViewer::EnableClipperWidget();
  auto qGLW = dynamic_cast<QVTKOpenGLNativeWidget*>(fGLWidget);
  qGLW->interactor()->Initialize();
}

void G4VtkQtViewer::SetWidgetInteractor(vtkAbstractWidget* widget)
{
  auto qGLW = dynamic_cast<QVTKOpenGLNativeWidget*>(fGLWidget);
  widget->SetInteractor(qGLW->interactor());
}
