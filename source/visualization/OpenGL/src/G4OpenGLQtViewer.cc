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
//
// $Id$
//
// 
// G4OpenGLQtViewer : Class to provide Qt specific
//                     functionality for OpenGL in GEANT4
//
// 27/06/2003 : G.Barrand : implementation (at last !).

#ifdef G4VIS_BUILD_OPENGLQT_DRIVER

#include "G4OpenGLQtViewer.hh"

#include "G4OpenGLSceneHandler.hh"
#include "G4VSolid.hh"
#include "G4OpenGLQtExportDialog.hh"
#include "G4OpenGLQtMovieDialog.hh"
#include "G4Qt.hh"
#include "G4UIQt.hh"
#include "G4UImanager.hh"
#include "G4UIcommandTree.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VisCommandsGeometrySet.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4Text.hh"

#include <CLHEP/Units/SystemOfUnits.h>

#include <typeinfo>

#include <qlayout.h>
#include <qlabel.h>
#include <qdialog.h>
#include <qpushbutton.h>
#include <qprocess.h>
#include <qapplication.h>
#include <qdesktopwidget.h>

#include <qmenu.h>
#include <qimagewriter.h>

#include <qtextedit.h>
#include <qtreewidget.h>
#include <qapplication.h>
#include <qmessagebox.h>
#include <qfiledialog.h>
#include <qprinter.h>
#include <qdatetime.h>
#include <qpainter.h>
#include <qgl.h> // include <qglwidget.h>
#include <qdialog.h>
#include <qcolordialog.h>
#include <qevent.h> //include <qcontextmenuevent.h>
#include <qobject.h>
#include <qgroupbox.h>
#include <qcombobox.h>
#include <qlineedit.h>
#include <qsignalmapper.h>

//////////////////////////////////////////////////////////////////////////////
void G4OpenGLQtViewer::CreateMainWindow (
 QGLWidget* glWidget
 ,const QString& name
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{

  if(fWindow) return; //Done.

  fWindow = glWidget ;
  //  fWindow->makeCurrent();

  G4Qt* interactorManager = G4Qt::getInstance ();

  ResizeWindow(fVP.GetWindowSizeHintX(),fVP.GetWindowSizeHintY());
    
  // FIXME L.Garnier 9/11/09 Has to be check !!! 
  // Qt UI with Qt Vis
  // Qt UI with X Vis
  // X UI with Qt Vis
  // X UI with X Vis
  // Ne marche pas avec un UIBatch !! (ecran blanc)

  // return false if G4UIQt was not launch

  G4UImanager* UI = G4UImanager::GetUIpointer();
  if (UI == NULL) return;

  if (! static_cast<G4UIQt*> (UI->GetG4UIWindow())) {
    // NO UI, should be batch mode
    fBatchMode = true;
    return;
  }
  fUiQt = static_cast<G4UIQt*> (UI->GetG4UIWindow());
  
  bool isTabbedView = false;
  if ( fUiQt) {
    if (!fBatchMode) {
      if (!interactorManager->IsExternalApp()) {
        isTabbedView = fUiQt->AddTabWidget(fWindow,name,getWinWidth(),getWinHeight());
        fUISceneTreeComponentsTBWidget = fUiQt->GetSceneTreeComponentsTBWidget();
        isTabbedView = true;
      }
    }
  }
#ifdef G4DEBUG_VIS_OGL
  else {
    printf("G4OpenGLQtViewer::CreateMainWindow :: UIQt NOt found \n");
  }
#endif

  if (!isTabbedView) { // we have to do a dialog

    QWidget *myParent = getParentWidget();
#ifdef G4DEBUG_VIS_OGL
    printf("G4OpenGLQtViewer::CreateMainWindow :: getParent OK \n");
#endif
    if (myParent != NULL) {
      glWidget->setParent(myParent);  
    }
    QHBoxLayout *mainLayout = new QHBoxLayout(fGLWindow);
    
    mainLayout->setMargin(0);
    mainLayout->setSpacing(0);   
    mainLayout->addWidget(fWindow);
    if (fGLWindow->inherits("QMainWindow")) {
      fGLWindow->setWindowTitle( name);
    }
    fGLWindow->setLayout(mainLayout);

    
    //useful for MACOSX, we have to compt the menuBar height
    int offset = QApplication::desktop()->height() 
      - QApplication::desktop()->availableGeometry().height();
    
    G4int YPos= fVP.GetWindowAbsoluteLocationHintY(QApplication::desktop()->height());
    if (fVP.GetWindowAbsoluteLocationHintY(QApplication::desktop()->height())< offset) {
      YPos = offset;
    }
    fGLWindow->resize(getWinWidth(), getWinHeight());
#ifdef G4DEBUG_VIS_OGL
    printf("G4OpenGLQtViewer::CreateMainWindow :: resizing to %d %d \n",getWinWidth(), getWinHeight());
#endif
    fGLWindow->move(fVP.GetWindowAbsoluteLocationHintX(QApplication::desktop()->width()),YPos);
    fGLWindow->show();
  } else {
    fGLWindow = fWindow;
    fGLWindow->resize(getWinWidth(), getWinHeight());
  }

  if(!fWindow) return;
  
  if (!fContextMenu) 
    createPopupMenu();

}


//////////////////////////////////////////////////////////////////////////////
G4OpenGLQtViewer::G4OpenGLQtViewer (
 G4OpenGLSceneHandler& scene
)
  :G4VViewer (scene, -1)
  ,G4OpenGLViewer (scene)
  ,fWindow(0)
  ,fRecordFrameNumber(0)
  ,fContextMenu(0)
  ,fDeltaDepth(0.01)
  ,fDeltaZoom(0.05)
  ,fHoldKeyEvent(false)
  ,fHoldMoveEvent(false)
  ,fHoldRotateEvent(false)
  ,fAutoMove(false)
  ,fEncoderPath("")
  ,fTempFolderPath("")
  ,fMovieTempFolderPath("")
  ,fSaveFileName("")
  ,fParameterFileName("ppmtompeg_encode_parameter_file.par")
  ,fMovieParametersDialog(NULL)
  ,fRecordingStep(WAIT)
  ,fProcess(NULL)
  ,fNbMaxFramesPerSec(100)
  ,fNbMaxAnglePerSec(360)
  ,fLaunchSpinDelay(100)
  ,fUISceneTreeComponentsTBWidget(NULL)
  ,fNoKeyPress(true)
  ,fAltKeyPress(false)
  ,fControlKeyPress(false)
  ,fShiftKeyPress(false)
  ,fBatchMode(false)
  ,fCheckSceneTreeComponentSignalLock(false)
  ,fSceneTreeComponentTreeWidget(NULL)
  ,fOldSceneTreeComponentTreeWidget(NULL)
  ,fSceneTreeWidget(NULL)
  ,fPVRootNodeCreate(false)
  ,fHelpLine(NULL)
  ,fNbRotation(0)
  ,fTimeRotation(0)
  ,fTouchableVolumes("Touchables")
  ,fShortcutsDialog(NULL)
  ,fSceneTreeDepthSlider(NULL)
  ,fSceneTreeDepth(1)
  ,fModelShortNameItem(NULL)
  ,fMaxPOindexInserted(-1)
  ,fUiQt(NULL)
  ,signalMapperMouse(NULL)
  ,signalMapperSurface(NULL)
{

  // launch Qt if not
  if (QCoreApplication::instance () == NULL) {
    fBatchMode = true;
  }
  G4Qt::getInstance ();

  fLastPos3 = QPoint(-1,-1);    
  fLastPos2 = QPoint(-1,-1);    
  fLastPos1 = QPoint(-1,-1);    
  
  initMovieParameters();

  fLastEventTime = new QTime();
  signalMapperMouse = new QSignalMapper(this);
  signalMapperSurface = new QSignalMapper(this);

#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLQtViewer::G4OpenGLQtViewer END\n");
#endif
}

//////////////////////////////////////////////////////////////////////////////
G4OpenGLQtViewer::~G4OpenGLQtViewer (
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  // remove scene tree from layout
  // Delete all the existing buttons in the layout
  QLayoutItem *wItem;
  if (fSceneTreeWidget != NULL) {
    if (fSceneTreeWidget->layout() != NULL) {
      while ((wItem = fSceneTreeWidget->layout()->takeAt(0)) != 0) {
	delete wItem->widget();
	delete wItem;
      }
    }
    if (fUISceneTreeComponentsTBWidget != NULL) {
      fUISceneTreeComponentsTBWidget->removeTab(fUISceneTreeComponentsTBWidget->indexOf(fSceneTreeWidget));
    }
  }
  G4cout <<removeTempFolder().toStdString().c_str() <<G4endl;
}


//
//   Create a popup menu for the widget. This menu is activated by right-mouse click
//
void G4OpenGLQtViewer::createPopupMenu()    {

  fContextMenu = new QMenu("All");

  QMenu *mMouseAction = fContextMenu->addMenu("&Mouse actions");

  fMouseRotateAction = mMouseAction->addAction("Rotate", signalMapperMouse, SLOT(map()));
  fMouseMoveAction = mMouseAction->addAction("Move", signalMapperMouse, SLOT(map()));
  fMousePickAction = mMouseAction->addAction("Pick", signalMapperMouse, SLOT(map()));
  fMouseZoomOutAction = mMouseAction->addAction("Zoom out", signalMapperMouse, SLOT(map()));
  fMouseZoomInAction = mMouseAction->addAction("Zoom in", signalMapperMouse, SLOT(map()));
  QAction *shortcutsAction = mMouseAction->addAction("Show shortcuts");

  fMouseRotateAction->setCheckable(true);
  fMouseMoveAction->setCheckable(true);
  fMousePickAction->setCheckable(true);
  fMouseZoomOutAction->setCheckable(true);
  fMouseZoomInAction->setCheckable(true);
  shortcutsAction->setCheckable(false);

  connect(signalMapperMouse, SIGNAL(mapped(int)),this, SLOT(toggleMouseAction(int)));
  signalMapperMouse->setMapping(fMouseRotateAction,1);
  signalMapperMouse->setMapping(fMouseMoveAction,2);
  signalMapperMouse->setMapping(fMousePickAction,3);
  signalMapperMouse->setMapping(fMouseZoomOutAction,4);
  signalMapperMouse->setMapping(fMouseZoomInAction,5);

  QObject::connect(shortcutsAction, 
                    SIGNAL(triggered(bool)),
                    this, 
                    SLOT(showShortcuts()));

  // === Style Menu ===
  QMenu *mStyle = fContextMenu->addMenu("&Style");

  QMenu *mRepresentation = mStyle->addMenu("&Representation");
  QMenu *mProjection = mStyle->addMenu("&Projection");
  QAction *polyhedron = mRepresentation->addAction("Polyhedron");
  QAction *nurbs = mRepresentation->addAction("NURBS");

  fProjectionOrtho = mProjection->addAction("Orthographic", signalMapperSurface, SLOT(map()));
  fProjectionPerspective = mProjection->addAction("Persepective", signalMapperSurface, SLOT(map()));

  // INIT mRepresentation
  G4ViewParameters::RepStyle style;
  style = fVP.GetRepStyle();
  if (style == G4ViewParameters::polyhedron) {
    createRadioAction(polyhedron,nurbs,SLOT(toggleRepresentation(bool)),1);
  } else if (style == G4ViewParameters::nurbs) {
    createRadioAction(polyhedron,nurbs,SLOT(toggleRepresentation(bool)),2);
  } else {
    mRepresentation->clear();
  }

  // INIT mProjection
  if (fVP.GetFieldHalfAngle() == 0) {
    createRadioAction(fProjectionOrtho, fProjectionPerspective,SLOT(toggleProjection(bool)),1);
  } else {
    createRadioAction(fProjectionOrtho, fProjectionPerspective,SLOT(toggleProjection(bool)),2);
  }

  // === Drawing Menu ===
  QMenu *mDrawing = mStyle->addMenu("&Drawing");

  fDrawingWireframe = mDrawing->addAction("Wireframe", signalMapperSurface, SLOT(map()));

  fDrawingLineRemoval = mDrawing->addAction("Hidden line removal", signalMapperSurface, SLOT(map()));

  fDrawingSurfaceRemoval = mDrawing->addAction("Hidden Surface removal", signalMapperSurface, SLOT(map()));

  fDrawingLineSurfaceRemoval = mDrawing->addAction("Hidden line and surface removal", signalMapperSurface, SLOT(map()));

  fDrawingWireframe->setCheckable(true);
  fDrawingLineRemoval->setCheckable(true);
  fDrawingSurfaceRemoval->setCheckable(true);
  fDrawingLineSurfaceRemoval->setCheckable(true);

  connect(signalMapperSurface, SIGNAL(mapped(int)),this, SLOT(toggleSurfaceAction(int)));
  signalMapperSurface->setMapping(fDrawingWireframe,1);
  signalMapperSurface->setMapping(fDrawingLineRemoval,2);
  signalMapperSurface->setMapping(fDrawingSurfaceRemoval,3);
  signalMapperSurface->setMapping(fDrawingLineSurfaceRemoval,4);


  // Background Color

  QAction *backgroundColorChooser ;
  // === Action Menu ===
  backgroundColorChooser = mStyle->addAction("Background color");
  QObject ::connect(backgroundColorChooser, 
                    SIGNAL(triggered()),
                    this,
                    SLOT(actionChangeBackgroundColor()));

  // Text Color

  QAction *textColorChooser ;
  // === Action Menu ===
  textColorChooser = mStyle->addAction("Text color");
  QObject ::connect(textColorChooser, 
                    SIGNAL(triggered()),
                    this,
                    SLOT(actionChangeTextColor()));

  // Default Color

  QAction *defaultColorChooser ;
  // === Action Menu ===
  defaultColorChooser = mStyle->addAction("Default color");
  QObject ::connect(defaultColorChooser, 
                    SIGNAL(triggered()),
                    this,
                    SLOT(actionChangeDefaultColor()));


  // === Action Menu ===
  QMenu *mActions = fContextMenu->addMenu("&Actions");
  QAction *createEPS = mActions->addAction("Save as ...");
  QObject ::connect(createEPS, 
                    SIGNAL(triggered()),
                    this,
                    SLOT(actionSaveImage()));

  // === Action Menu ===
  QAction *movieParameters = mActions->addAction("Movie parameters...");
  QObject ::connect(movieParameters, 
                    SIGNAL(triggered()),
                    this,
                    SLOT(actionMovieParameters()));




  // === Special Menu ===
  QMenu *mSpecial = fContextMenu->addMenu("S&pecial");
  QMenu *mTransparency = mSpecial->addMenu("Transparency");
  QAction *transparencyOn = mTransparency->addAction("On");
  QAction *transparencyOff = mTransparency->addAction("Off");

  if (transparency_enabled == false) {
    createRadioAction(transparencyOn,transparencyOff,SLOT(toggleTransparency(bool)),2);
  } else if (transparency_enabled == true) {
    createRadioAction(transparencyOn,transparencyOff,SLOT(toggleTransparency(bool)),1);
  } else {
    mSpecial->clear();
  }


  QMenu *mAntialiasing = mSpecial->addMenu("Antialiasing");
  QAction *antialiasingOn = mAntialiasing->addAction("On");
  QAction *antialiasingOff = mAntialiasing->addAction("Off");

  if (antialiasing_enabled == false) {
    createRadioAction(antialiasingOn,antialiasingOff,SLOT(toggleAntialiasing(bool)),2);
  } else if (antialiasing_enabled == true) {
    createRadioAction(antialiasingOn,antialiasingOff,SLOT(toggleAntialiasing(bool)),1);
  } else {
    mAntialiasing->clear();
  }

  QMenu *mHaloing = mSpecial->addMenu("Haloing");
  QAction *haloingOn = mHaloing->addAction("On");
  QAction *haloingOff = mHaloing->addAction("Off");
  if (haloing_enabled == false) {
    createRadioAction(haloingOn,haloingOff,SLOT(toggleHaloing(bool)),2);
  } else if (haloing_enabled == true) {
    createRadioAction(haloingOn,haloingOff,SLOT(toggleHaloing(bool)),1);
  } else {
    mHaloing->clear();
  }

  QMenu *mAux = mSpecial->addMenu("Auxiliary edges");
  QAction *auxOn = mAux->addAction("On");
  QAction *auxOff = mAux->addAction("Off");
  if (!fVP.IsAuxEdgeVisible()) {
    createRadioAction(auxOn,auxOff,SLOT(toggleAux(bool)),2);
  } else {
    createRadioAction(auxOn,auxOff,SLOT(toggleAux(bool)),1);
  }


  QMenu *mHiddenMarkers = mSpecial->addMenu("Hidden markers");
  QAction *hiddenMarkersOn = mHiddenMarkers->addAction("On");
  QAction *hiddenMarkersOff = mHiddenMarkers->addAction("Off");
  if (fVP.IsMarkerNotHidden()) {
    createRadioAction(hiddenMarkersOn,hiddenMarkersOff,SLOT(toggleHiddenMarkers(bool)),2);
  } else {
    createRadioAction(hiddenMarkersOn,hiddenMarkersOff,SLOT(toggleHiddenMarkers(bool)),1);
  }



  QMenu *mFullScreen = mSpecial->addMenu("&Full screen");
  fFullScreenOn = mFullScreen->addAction("On");
  fFullScreenOff = mFullScreen->addAction("Off");
  createRadioAction(fFullScreenOn,fFullScreenOff,SLOT(toggleFullScreen(bool)),2);

  // INIT All
  updateToolbarAndMouseContextMenu();
}


void G4OpenGLQtViewer::G4manageContextMenuEvent(QContextMenuEvent *e)
{
  if (!fGLWindow) {
    G4cerr << "Visualization window not defined, please choose one before" << G4endl;
  } else {
  
    if (!fContextMenu) 
      createPopupMenu();
    
    // launch menu
    if ( fContextMenu ) {
      fContextMenu->exec( e->globalPos() );
      //    delete fContextMenu;
    }
  }
  e->accept();
}


/**
   Create a radio button menu. The two menu will be connected. When click on one,
   eatch state will be invert and callback method will be called.
   @param action1 first action to connect
   @param action2 second action to connect
   @param method callback method
   @param nCheck: 1 : first action will be set true. 2 : second action will be set true
*/
void G4OpenGLQtViewer::createRadioAction(QAction *action1,QAction *action2, const std::string& method,unsigned int nCheck) {

  action1->setCheckable(true);
  action2->setCheckable(true);

  if (nCheck ==1)
    action1->setChecked (true);
  else
    action2->setChecked (true);
   
  QObject ::connect(action1, SIGNAL(triggered(bool)),action2, SLOT(toggle()));
  QObject ::connect(action2, SIGNAL(triggered(bool)),action1, SLOT(toggle()));

  QObject ::connect(action1, SIGNAL(toggled(bool)),this, method.c_str());

}



/**
   Show shortcuts for this mouse action
*/
void G4OpenGLQtViewer::showShortcuts() {
  G4String text;

  text = "========= Mouse Shortcuts =========\n";
  if (fUiQt != NULL) {
    if (fUiQt->IsIconRotateSelected()) {  // rotate
      text += "Click and move mouse to rotate volume \n";
      text += "ALT + Click and move mouse to rotate volume (Toggle View/Theta-Phi Direction) \n";
      text += "CTRL + Click and zoom mouse to zoom in/out \n";
      text += "SHIFT + Click and zoommove camera point of view \n";
    } else  if (fUiQt->IsIconMoveSelected()) { //move
      text += "Move camera point of view with mouse \n";
    } else  if (fUiQt->IsIconPickSelected()) { //pick
      text += "Click and pick \n";
    }
  } else {
    text += "Click and move mouse to rotate volume \n";
    text += "ALT + Click and move mouse to rotate volume (Toggle View/Theta-Phi Direction) \n";
    text += "CTRL + Click and zoom mouse to zoom in/out \n";
    text += "SHIFT + Click and zoommove camera point of view \n";
  }
  text += "========= Move Shortcuts ========= \n";
  text += "Press left/right arrows to move volume left/right \n";
  text += "Press up/down arrows to move volume up/down \n";
  text += "Press '+'/'-' to move volume toward/forward \n";
  text += "\n";
  text += "========= Rotation (Theta/Phi) Shortcuts ========= \n";
  text += "Press SHIFT + left/right arrows to rotate volume left/right \n";
  text += "Press SHIFT + up/down arrows to rotate volume up/down \n";
  text += "\n";
  text += "========= Rotation (View Direction) Shortcuts ========= \n";
  text += "Press ALT + left/right to rotate volume around vertical direction \n";
  text += "Press ALT + up/down to rotate volume around horizontal direction \n";
  text += "\n";
  text += "========= Zoom View ========= \n";
  text += "Press CTRL + '+'/'-' to zoom into volume \n";
  text += "\n";
  text += "========= Misc ========= \n";
  text += "Press ALT +/- to slow/speed rotation/move \n";
  text += "Press H to reset view \n";
  text += "Press Esc to exit FullScreen \n";
  text += "\n";
  text += "========= Video ========= \n";
  text += "In video mode : \n";
  text += " Press SPACE to Start/Pause video recording \n";
  text += " Press RETURN to Stop video recording \n";
  text += "\n";

  G4cout << text;

  if (  fShortcutsDialog == NULL) {
    fShortcutsDialog = new QDialog();
    fShortcutsDialogInfos = new QTextEdit() ;
    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(fShortcutsDialogInfos);
    fShortcutsDialog->setLayout(mainLayout);
    fShortcutsDialog->setWindowTitle(tr("Shortcuts"));
  }
  
  fShortcutsDialogInfos->setPlainText(text.data());
  fShortcutsDialog->show();      
}



/**
   Slot activated when mouse action is toggle
   @param aAction : 1 rotate, 2 move, 3 pick, 4 zoom out, 5 zoom in
   @see G4OpenGLStoredQtViewer::DrawView
   @see G4XXXStoredViewer::CompareForKernelVisit
*/
void G4OpenGLQtViewer::toggleMouseAction(int aAction) {

  if (aAction == 1) {
    fUiQt->SetIconRotateSelected();
  } else  if (aAction == 2) {
    fUiQt->SetIconMoveSelected();
  } else  if (aAction == 3) {
    fUiQt->SetIconPickSelected();
  } else  if (aAction == 4) {
    fUiQt->SetIconZoomOutSelected();
  } else  if (aAction == 5) {
    fUiQt->SetIconZoomInSelected();
  }

  updateToolbarAndMouseContextMenu();
}


/**
   Slot activated when drawing menu is toggle
   Warning : When G4OpenGLStoredQtViewer::DrawView() method call,
   KernelVisitDecision () will be call and will set the fNeedKernelVisit
   to 1. See G4XXXStoredViewer::CompareForKernelVisit for explanations.
   It will cause a redraw of the view
   @param aAction : 1 wireframe, 2 line removal, 3 surface removal, 4 line & surface removal
   @see G4OpenGLStoredQtViewer::DrawView
   @see G4XXXStoredViewer::CompareForKernelVisit
*/
void G4OpenGLQtViewer::toggleSurfaceAction(int aAction) {

  G4ViewParameters::DrawingStyle d_style = G4ViewParameters::wireframe;
  
  if (aAction ==1) {
    d_style = G4ViewParameters::wireframe;

  } else  if (aAction ==2) {
    d_style = G4ViewParameters::hlr;

  } else  if (aAction ==3) {
    d_style = G4ViewParameters::hsr;

  } else  if (aAction ==4) {
    d_style = G4ViewParameters::hlhsr;
  }
  fVP.SetDrawingStyle(d_style);

  updateToolbarAndMouseContextMenu();
  updateQWidget();
}


/**
   SLOT Activate by a click on the representation menu
   Warning : When G4OpenGLStoredQtViewer::DrawView() method call,
   KernelVisitDecision () will be call and will set the fNeedKernelVisit
   to 1. See G4XXXStoredViewer::CompareForKernelVisit for explanations.
   It will cause a redraw of the view
   @param check : 1 polyhedron, 0 nurbs
   @see G4OpenGLStoredQtViewer::DrawView
   @see G4XXXStoredViewer::CompareForKernelVisit
*/
void G4OpenGLQtViewer::toggleRepresentation(bool check) {

  G4ViewParameters::RepStyle style;
  if (check == 1) {
    style = G4ViewParameters::polyhedron;
  } else {
    style = G4ViewParameters::nurbs;
  }
  fVP.SetRepStyle (style);

  updateToolbarAndMouseContextMenu();
  updateQWidget();
}

/**
   SLOT Activate by a click on the projection menu
   Warning : When G4OpenGLStoredQtViewer::DrawView() method call,
   KernelVisitDecision () will be call and will set the fNeedKernelVisit
   to 1. See G4XXXStoredViewer::CompareForKernelVisit for explanations.
   It will cause a redraw of the view
   @param check : 1 orthographic, 2 perspective
   @see G4OpenGLStoredQtViewer::DrawView
   @see G4XXXStoredViewer::CompareForKernelVisit
*/
void G4OpenGLQtViewer::toggleProjection(bool check) {

  if (check == 1) {
    fVP.SetOrthogonalProjection ();
  } else {
    fVP.SetPerspectiveProjection();
  }  
  updateToolbarAndMouseContextMenu();
  updateQWidget();
}


/**
   SLOT Activate by a click on the transparency menu
   @param check : 1 , 0
*/
void G4OpenGLQtViewer::toggleTransparency(bool check) {
  
  if (check) {
    transparency_enabled = true;
  } else {
    transparency_enabled = false;
  }
  SetNeedKernelVisit (true);
  updateToolbarAndMouseContextMenu();
  updateQWidget();
}

/**
   SLOT Activate by a click on the antialiasing menu
   @param check : 1 , 0
*/
void G4OpenGLQtViewer::toggleAntialiasing(bool check) {

  if (!check) {
    antialiasing_enabled = false;
    glDisable (GL_LINE_SMOOTH);
    glDisable (GL_POLYGON_SMOOTH);
  } else {
    antialiasing_enabled = true;
    glEnable (GL_LINE_SMOOTH);
    glHint (GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable (GL_POLYGON_SMOOTH);
    glHint (GL_POLYGON_SMOOTH_HINT, GL_NICEST);
  }

  updateToolbarAndMouseContextMenu();
  updateQWidget();
}

/**
   SLOT Activate by a click on the haloing menu
   @param check : 1 , 0
*/
//FIXME : I SEE NOTHING...
void G4OpenGLQtViewer::toggleHaloing(bool check) {
  if (check) {
    haloing_enabled = false;
  } else {
    haloing_enabled = true;
  }

  updateToolbarAndMouseContextMenu();
  updateQWidget();

}

/**
   SLOT Activate by a click on the auxiliaire edges menu
   @param check : 1 , 0
*/
void G4OpenGLQtViewer::toggleAux(bool check) {
  if (check) {
    fVP.SetAuxEdgeVisible(true);
  } else {
    fVP.SetAuxEdgeVisible(false);
  }
  SetNeedKernelVisit (true);
  updateToolbarAndMouseContextMenu();
  updateQWidget();
}

/**
   SLOT Activate by a click on the hidden marker menu
   @param check : 1 , 0
*/
void G4OpenGLQtViewer::toggleHiddenMarkers(bool check) {
  if (check) {
    fVP.SetMarkerHidden();
  } else {
    fVP.SetMarkerNotHidden();
  }
  //  SetNeedKernelVisit (true);
  updateToolbarAndMouseContextMenu();
  updateQWidget();
}

/**
   SLOT Activate by a click on the full screen menu
*/
void G4OpenGLQtViewer::toggleFullScreen(bool check) {
  if (check != fGLWindow->isFullScreen()) { //toggle
    fGLWindow->setWindowState(fGLWindow->windowState() ^ Qt::WindowFullScreen);
  }
}


void G4OpenGLQtViewer::savePPMToTemp() {
  if (fMovieTempFolderPath == "") {
    return;
  }
  QString fileName ="Test"+QString::number(fRecordFrameNumber)+".ppm";
  QString filePath =fMovieTempFolderPath+fileName;

  QImage image;
  image = fWindow->grabFrameBuffer();
  bool res = false;
  
  res = image.save(filePath,0);
  if (res == false) { 
    resetRecording();
    setRecordingInfos("Can't save tmp file "+filePath);
    return;
  }
  
  setRecordingInfos("File "+fileName+" saved");
  fRecordFrameNumber++;
}



void G4OpenGLQtViewer::actionSaveImage() {
  QString filters;
  QList<QByteArray> formats =  QImageWriter::supportedImageFormats ();
  for (int i = 0; i < formats.size(); ++i) {
    filters +=formats.at(i) + ";;";
  }
  filters += "eps;;";
  filters += "ps;;";
  filters += "pdf";
  QString* selectedFormat = new QString();
  std::string name;
  name =  QFileDialog::getSaveFileName ( fGLWindow,
                                         tr("Save as ..."),
                                         ".",
                                         filters,
                                         selectedFormat ).toStdString().c_str(); 
  // bmp jpg jpeg png ppm xbm xpm
  if (name.empty()) {
    return;
  }
  name += "." + selectedFormat->toStdString();
  QString format = selectedFormat->toLower();
  setPrintFilename(name.c_str(),0);
  G4OpenGLQtExportDialog* exportDialog= new G4OpenGLQtExportDialog(fGLWindow,format,fWindow->height(),fWindow->width());
  if(  exportDialog->exec()) {

    QImage image;
    bool res = false;
    if ((exportDialog->getWidth() !=fWindow->width()) ||
        (exportDialog->getHeight() !=fWindow->height())) {
      setPrintSize(exportDialog->getWidth(),exportDialog->getHeight());
      if ((format != QString("eps")) && (format != QString("ps"))) {
        G4cerr << "Export->Change Size : This function is not implemented, to export in another size, please resize your frame to what you need" << G4endl;
      
        //    rescaleImage(exportDialog->getWidth(),exportDialog->getHeight());// re-scale image
        //      QGLWidget* glResized = fWindow;

        // FIXME :
        // L.Garnier : I've try to implement change size function, but the problem is 
        // the renderPixmap function call the QGLWidget to resize and it doesn't draw
        // the content of this widget... It only draw the background.

        //      fWindow->renderPixmap (exportDialog->getWidth()*2,exportDialog->getHeight()*2,true );

        //      QPixmap pixmap = fWindow->renderPixmap ();
      
        //      image = pixmap->toImage();
        //      glResized->resize(exportDialog->getWidth()*2,exportDialog->getHeight()*2);
        //      image = glResized->grabFrameBuffer();
      }      
    } else {
      image = fWindow->grabFrameBuffer();
    }    
    if (format == QString("eps")) {
      fVectoredPs = exportDialog->getVectorEPS();
      printEPS();
    } else if (format == "ps") {
      fVectoredPs = true;
      printEPS();
    } else if (format == "pdf") {

      res = printPDF(name,exportDialog->getNbColor(),image);

    } else if ((format == "tif") ||
               (format == "tiff") ||
               (format == "jpg") ||
               (format == "jpeg") ||
               (format == "png") ||
               (format == "pbm") ||
               (format == "pgm") ||
               (format == "ppm") ||
               (format == "bmp") ||
               (format == "xbm") ||
               (format == "xpm")) {
      res = image.save(QString(name.c_str()),0,exportDialog->getSliderValue());
    } else {
      G4cerr << "This version of G4UI Could not generate the selected format" << G4endl;
    }
    if ((format == QString("eps")) && (format == QString("ps"))) {
      if (res == false) {
        G4cerr << "Error while saving file... "<<name.c_str()<< G4endl;
      } else {
        G4cout << "File "<<name.c_str()<<" has been saved " << G4endl;
      }
    }
    
  } else { // cancel selected
    return;
  }
  
}


void G4OpenGLQtViewer::actionChangeBackgroundColor() {

  //   //I need to revisit the kernel if the background colour changes and
  //   //hidden line removal is enabled, because hlr drawing utilises the
  //   //background colour in its drawing...
  //   // (Note added by JA 13/9/2005) Background now handled in view
  //   // parameters.  A kernel visit is triggered on change of background.

#if QT_VERSION < 0x040500
  bool a;
  const QColor color = QColor(QColorDialog::getRgba (QColor(Qt::black).rgba(),&a,fGLWindow));
#else
  const QColor color =
    QColorDialog::getColor(Qt::black,
                           fGLWindow,
                           " Get background color and transparency",
                           QColorDialog::ShowAlphaChannel);
#endif
  if (color.isValid()) {
    G4Colour colour(((G4double)color.red())/255,
                    ((G4double)color.green())/255,
                    ((G4double)color.blue())/255,
                    ((G4double)color.alpha())/255);
    fVP.SetBackgroundColour(colour);

    updateToolbarAndMouseContextMenu();
    updateQWidget();
  }
}

void G4OpenGLQtViewer::actionChangeTextColor() {

#if QT_VERSION < 0x040500
  bool a;
  const QColor color = QColor(QColorDialog::getRgba (QColor(Qt::yellow).rgba(),&a,fGLWindow));
#else
  const QColor& color =
    QColorDialog::getColor(Qt::yellow,
                           fGLWindow,
                           " Get text color and transparency",
                           QColorDialog::ShowAlphaChannel);
#endif
  if (color.isValid()) {
    G4Colour colour(((G4double)color.red())/255,
                    ((G4double)color.green())/255,
                    ((G4double)color.blue())/255,
                    ((G4double)color.alpha())/255);

    fVP.SetDefaultTextColour(colour);

    updateToolbarAndMouseContextMenu();
    updateQWidget();
  }
}

void G4OpenGLQtViewer::actionChangeDefaultColor() {

#if QT_VERSION < 0x040500
  bool a;
  const QColor color = QColor(QColorDialog::getRgba (QColor(Qt::white).rgba(),&a,fGLWindow));
#else
  const QColor& color =
    QColorDialog::getColor(Qt::white,
                           fGLWindow,
                           " Get default color and transparency",
                           QColorDialog::ShowAlphaChannel);
#endif
  if (color.isValid()) {
    G4Colour colour(((G4double)color.red())/255,
                    ((G4double)color.green())/255,
                    ((G4double)color.blue())/255,
                    ((G4double)color.alpha())/255);

    fVP.SetDefaultColour(colour);

    updateToolbarAndMouseContextMenu();
    updateQWidget();
  }
}


void G4OpenGLQtViewer::actionMovieParameters() {
  showMovieParametersDialog();
}


void G4OpenGLQtViewer::showMovieParametersDialog() {
  if (!fMovieParametersDialog) {
    fMovieParametersDialog= new G4OpenGLQtMovieDialog(this,fGLWindow);
    displayRecordingStatus();
    fMovieParametersDialog->checkEncoderSwParameters();
    fMovieParametersDialog->checkSaveFileNameParameters();
    fMovieParametersDialog->checkTempFolderParameters();
    if (getEncoderPath() == "") {
      setRecordingInfos("ppmtompeg is needed to encode in video format. It is available here: http://netpbm.sourceforge.net ");
    }
  }
  fMovieParametersDialog->show();
}


/*
// http://www.google.com/codesearch?hl=en&q=+jpg+Qt+quality+QDialog+show:FZkUoth8oiw:TONpW2mR-_c:tyTfrKMO-xI&sa=N&cd=2&ct=rc&cs_p=http://soft.proindependent.com/src/qtiplot-0.8.9.zip&cs_f=qtiplot-0.8.9/qtiplot/src/application.cpp#a0

void Graph::exportToSVG(const QString& fname)
{
// enable workaround for Qt3 misalignments
QwtPainter::setSVGMode(true);
QPicture picture;
QPainter p(&picture);
d_plot->print(&p, d_plot->rect());
p.end();

picture.save(fname, "svg");
}
*/



void G4OpenGLQtViewer::FinishView()
{
  glFlush ();

  // L. Garnier 10/2009 : Not necessary and cause problems on mac OS X 10.6
  //  fWindow->swapBuffers ();
}

/**
   Save the current mouse press point
   @param p mouse click point
*/
void G4OpenGLQtViewer::G4MousePressEvent(QMouseEvent *evnt)
{
  if ((evnt->buttons() & Qt::LeftButton) && (! (evnt->modifiers() & Qt::ControlModifier ))){
    fWindow->setMouseTracking(true);
    fAutoMove = false; // stop automove
    fLastPos1 = evnt->pos();
    fLastPos2 = fLastPos1;
    fLastPos3 = fLastPos2;
    fLastEventTime->start();
    if (fUiQt != NULL) {

      if (fUiQt->IsIconPickSelected()){  // pick
        fVP.SetPicking(true);
        Pick(evnt->pos().x(),evnt->pos().y());
        fVP.SetPicking(false);
 
      } else if (fUiQt->IsIconZoomInSelected()) {  // zoomIn
        // Move click point to center of OGL

        float deltaX = ((float)getWinWidth()/2-evnt->pos().x());
        float deltaY = ((float)getWinHeight()/2-evnt->pos().y());

        G4double coefTrans = 0;
        coefTrans = ((G4double)getSceneNearWidth())/((G4double)getWinWidth());
        if (getWinHeight() <getWinWidth()) {
          coefTrans = ((G4double)getSceneNearWidth())/((G4double)getWinHeight());
        }
        fVP.IncrementPan(-deltaX*coefTrans,deltaY*coefTrans,0);
        fVP.SetZoomFactor(1.5 * fVP.GetZoomFactor());
  
        updateQWidget();
        
      } else if (fUiQt->IsIconZoomOutSelected()) {  // zoomOut
        // Move click point to center of OGL
        moveScene(((float)getWinWidth()/2-evnt->pos().x()),((float)getWinHeight()/2-evnt->pos().y()),0,true);

        fVP.SetZoomFactor(0.75 * fVP.GetZoomFactor());
        updateQWidget();

      }
    }
  }
}

/**
 */
void G4OpenGLQtViewer::G4MouseReleaseEvent()
{
  fSpinningDelay = fLastEventTime->elapsed();
  QPoint delta = (fLastPos3-fLastPos1);
  if ((delta.x() == 0) && (delta.y() == 0)) {
    return;
  }
  if (fSpinningDelay < fLaunchSpinDelay ) {
    fAutoMove = true;
    QTime lastMoveTime;
    lastMoveTime.start();
    // try to addapt speed move/rotate looking to drawing speed
    float correctionFactor = 5;
    while (fAutoMove) {
      if ( lastMoveTime.elapsed() >= (int)(1000/fNbMaxFramesPerSec)) {
        float lTime = 1000/lastMoveTime.elapsed();
        if (((((float)delta.x())/correctionFactor)*lTime > fNbMaxAnglePerSec) ||
            ((((float)delta.x())/correctionFactor)*lTime < -fNbMaxAnglePerSec) ) {
          correctionFactor = (float)delta.x()*(lTime/fNbMaxAnglePerSec);
          if (delta.x() <0 ) {
            correctionFactor = -correctionFactor;
          }
        }
        if (((((float)delta.y())/correctionFactor)*lTime > fNbMaxAnglePerSec) ||
            ((((float)delta.y())/correctionFactor)*lTime < -fNbMaxAnglePerSec) ) {
          correctionFactor = (float)delta.y()*(lTime/fNbMaxAnglePerSec);
          if (delta.y() <0 ) {
            correctionFactor = -correctionFactor;
          }
        }
                
        // Check Qt Versions for META Keys
                
        // Click and move mouse to rotate volume
        // ALT + Click and move mouse to rotate volume (View Direction)
        // SHIFT + Click and move camera point of view
        // CTRL + Click and zoom mouse to zoom in/out

        lastMoveTime.start();

        bool rotate = false;
        bool move = false;

        if (fUiQt != NULL) {
          if (fUiQt->IsIconRotateSelected()) {  // rotate
            rotate = true;
          } else if (fUiQt->IsIconMoveSelected()) {  // move
            move = true;
          }
        } else {
          rotate = true;
        }
        if (rotate) {  // rotate
          if (fNoKeyPress) {
            rotateQtScene(((float)delta.x())/correctionFactor,((float)delta.y())/correctionFactor);
          } else if (fAltKeyPress) {
            rotateQtSceneToggle(((float)delta.x())/correctionFactor,((float)delta.y())/correctionFactor);
          }
#ifdef G4DEBUG_VIS_OGL
          fNbRotation ++;
          fTimeRotation += lastMoveTime.elapsed();
          printf("G4OpenGLQtViewer %f \n",fTimeRotation/(float)fNbRotation);
#endif
          
        } else if (move) {  // move
          moveScene(-((float)delta.x())/correctionFactor,-((float)delta.y())/correctionFactor,0,true);
        }
      }
      ((QApplication*)G4Qt::getInstance ())->processEvents();
    }
  }
  fWindow->setMouseTracking(false);

}


void G4OpenGLQtViewer::G4MouseDoubleClickEvent()
{
  fWindow->setMouseTracking(true);
}


/**
   @param pos_x mouse x position
   @param pos_y mouse y position
   @param mButtons mouse button active
   @param mAutoMove true: apply this move till another evnt came, false :one time move
*/

void G4OpenGLQtViewer::G4MouseMoveEvent(QMouseEvent *evnt)
{

  Qt::MouseButtons mButtons = evnt->buttons();

  updateKeyModifierState(evnt->modifiers());

  if (fAutoMove) {
    return;
  }

  fLastPos3 = fLastPos2;
  fLastPos2 = fLastPos1;
  fLastPos1 = QPoint(evnt->x(), evnt->y());

  int deltaX = fLastPos2.x()-fLastPos1.x();
  int deltaY = fLastPos2.y()-fLastPos1.y();

  bool rotate = false;
  bool move = false;
  if (fUiQt != NULL) {
    if (fUiQt->IsIconRotateSelected()) {  // rotate
      rotate = true;
    } else if (fUiQt->IsIconMoveSelected()) {  // move
      move = true;
    }
  } else {
    rotate = true;
  }
  if (rotate) {  // rotate
    if (mButtons & Qt::LeftButton) {
      if (fNoKeyPress) {
        rotateQtScene(((float)deltaX),((float)deltaY));
      } else if (fAltKeyPress) {
        rotateQtSceneToggle(((float)deltaX),((float)deltaY));
      } else if (fShiftKeyPress) {
        unsigned int sizeWin;
        sizeWin = getWinWidth();
        if (getWinHeight() < getWinWidth()) {
          sizeWin = getWinHeight();
        }

        // L.Garnier : 08/2010 100 is the good value, but don't ask me why !
        float factor = ((float)100/(float)sizeWin) ;
        moveScene(-(float)deltaX*factor,-(float)deltaY*factor,0,false);
      } else if (fControlKeyPress) {
        fVP.SetZoomFactor(fVP.GetZoomFactor()*(1+((float)deltaY))); 
      }
    }
  } else if (move) {  // move
    if (mButtons & Qt::LeftButton) {
      moveScene(-(float)deltaX,-(float)deltaY,0,true);
    }
  }

  fLastEventTime->start();
}


/**
   Move the scene of dx, dy, dz values.
   @param dx delta mouse x position
   @param dy delta mouse y position
   @param mouseMove : true if event comes from a mouse move, false if event comes from key action
*/

void G4OpenGLQtViewer::moveScene(float dx,float dy, float dz,bool mouseMove)
{
  if (fHoldMoveEvent)
    return;
  fHoldMoveEvent = true;

  G4double coefTrans = 0;
  GLdouble coefDepth = 0;
  if(mouseMove) {
    coefTrans = ((G4double)getSceneNearWidth())/((G4double)getWinWidth());
    if (getWinHeight() <getWinWidth()) {
      coefTrans = ((G4double)getSceneNearWidth())/((G4double)getWinHeight());
    }
  } else {
    coefTrans = getSceneNearWidth()*fPan_sens;
    coefDepth = getSceneDepth()*fDeltaDepth;
  }
  fVP.IncrementPan(-dx*coefTrans,dy*coefTrans,dz*coefDepth);
  
  updateQWidget();
  if (fAutoMove)
    ((QApplication*)G4Qt::getInstance ())->processEvents();
  
  fHoldMoveEvent = false;
}


/**
   @param dx delta mouse x position
   @param dy delta mouse y position
*/

void G4OpenGLQtViewer::rotateQtScene(float dx, float dy)
{
  if (fHoldRotateEvent)
    return;
  fHoldRotateEvent = true;
  
  rotateScene(dx,dy);

  updateQWidget();
  
  fHoldRotateEvent = false;
}

/**
   @param dx delta mouse x position
   @param dy delta mouse y position
*/

void G4OpenGLQtViewer::rotateQtSceneToggle(float dx, float dy)
{
  if (fHoldRotateEvent)
    return;
  fHoldRotateEvent = true;
  
  rotateSceneToggle(dx,dy);

  updateQWidget();
  
  fHoldRotateEvent = false;
}





/** This is the benning of a rescale function. It does nothing for the moment
    @param aWidth : new width
    @param aHeight : new height
*/
void G4OpenGLQtViewer::rescaleImage(
 int /* aWidth */
,int /* aHeight */
){
  //  GLfloat* feedback_buffer;
  //  GLint returned;
  //  FILE* file;
  
  //   feedback_buffer = new GLfloat[size];
  //   glFeedbackBuffer (size, GL_3D_COLOR, feedback_buffer);
  //   glRenderMode (GL_FEEDBACK);
  
  //   DrawView();
  //   returned = glRenderMode (GL_RENDER);

}



/**
   Generate Postscript or PDF form image
   @param aFilename : name of file
   @param aInColor : numbers of colors : 1->BW 2->RGB
   @param aImage : Image to print
*/
bool G4OpenGLQtViewer::printPDF (
 const std::string aFilename
 ,int aInColor
 ,QImage aImage
)
{

  QPrinter printer;
  //  printer.setPageSize(pageSize);

  // FIXME : L. Garnier 4/12/07
  // This is not working, it does nothing. Image is staying in color mode
  // So I have desactivate the B/W button in GUI
  if ((!aImage.isGrayscale ()) &&(aInColor ==1 )) {
    aImage = aImage.convertToFormat ( aImage.format(), Qt::MonoOnly);
  }


  if (aFilename.substr(aFilename.size()-3) == ".ps") {
#if QT_VERSION > 0x040200
  #if QT_VERSION < 0x050000
    printer.setOutputFormat(QPrinter::PostScriptFormat);
  #else
    G4cout << "Output in postscript not implemented in this Qt version" << G4endl;
  #endif
#else
    G4cout << "Output in postscript not implemented in this Qt version" << G4endl;
#endif
  } else {
#if QT_VERSION > 0x040100
    printer.setOutputFormat(QPrinter::PdfFormat);
#endif
  }
#if QT_VERSION > 0x040100
  printer.setOutputFileName(QString(aFilename.c_str()));
#endif
  //  printer.setFullPage ( true);
  QPainter paint(&printer);
  paint.drawImage (0,0,aImage);
  paint.end();
  return true;
}


void G4OpenGLQtViewer::G4wheelEvent (QWheelEvent * evnt) 
{
  fVP.SetZoomFactor(fVP.GetZoomFactor()+(fVP.GetZoomFactor()*(evnt->delta())/1200)); 
  updateQWidget();
}


void G4OpenGLQtViewer::G4keyPressEvent (QKeyEvent * evnt) 
{
  if (fHoldKeyEvent)
    return;

  fHoldKeyEvent = true;

  
  // with no modifiers
  updateKeyModifierState(evnt->modifiers());
  if ((fNoKeyPress) || (evnt->modifiers() == Qt::KeypadModifier )) {
    if (evnt->key() == Qt::Key_Down) { // go down
      moveScene(0,1,0,false);
    }
    else if (evnt->key() == Qt::Key_Up) {  // go up
      moveScene(0,-1,0,false);
    }
    if (evnt->key() == Qt::Key_Left) { // go left
      moveScene(-1,0,0,false);
    }
    else if (evnt->key() == Qt::Key_Right) { // go right
      moveScene(1,0,0,false);
    }
    if (evnt->key() == Qt::Key_Minus) { // go backward
      moveScene(0,0,1,false);
    }
    else if (evnt->key() == Qt::Key_Plus) { // go forward
      moveScene(0,0,-1,false);
    }

    // escaped from full screen
    if (evnt->key() == Qt::Key_Escape) {
      toggleFullScreen(false);
    }
  }    
  // several case here : If return is pressed, in every case -> display the movie parameters dialog
  // If one parameter is wrong -> put it in red (only save filenam could be wrong..)
  // If encoder not found-> does nothing.Only display a message in status box
  // If all ok-> generate parameter file
  // If ok -> put encoder button enabled
  
  if ((evnt->key() == Qt::Key_Return) || (evnt->key() == Qt::Key_Enter)){ // end of video
    stopVideo();
  }
  if (evnt->key() == Qt::Key_Space){ // start/pause of video
    startPauseVideo();
  }
  
  // H : Return Home view
  if (evnt->key() == Qt::Key_H){ // go Home
    ResetView();

    updateQWidget();
  }

  // Shift Modifier
  if (fShiftKeyPress) {
    if (evnt->key() == Qt::Key_Down) { // rotate phi
      rotateQtScene(0,-fRot_sens);
    }
    else if (evnt->key() == Qt::Key_Up) { // rotate phi
      rotateQtScene(0,fRot_sens);
    }
    if (evnt->key() == Qt::Key_Left) { // rotate theta
      rotateQtScene(fRot_sens,0);
    }
    else if (evnt->key() == Qt::Key_Right) { // rotate theta
      rotateQtScene(-fRot_sens,0);
    }
    if (evnt->key() == Qt::Key_Plus) { // go forward  ("Plus" imply 
      // "Shift" on Mac French keyboard
      moveScene(0,0,-1,false);
    }

    // Alt Modifier
  }
  if ((fAltKeyPress)) {
    if (evnt->key() == Qt::Key_Down) { // rotate phi
      rotateQtSceneToggle(0,-fRot_sens);
    }
    else if (evnt->key() == Qt::Key_Up) { // rotate phi
      rotateQtSceneToggle(0,fRot_sens);
    }
    if (evnt->key() == Qt::Key_Left) { // rotate theta
      rotateQtSceneToggle(fRot_sens,0);
    }
    else if (evnt->key() == Qt::Key_Right) { // rotate theta
      rotateQtSceneToggle(-fRot_sens,0);
    }

    // Rotatio +/-
    if (evnt->key() == Qt::Key_Plus) {
      fRot_sens = fRot_sens/0.7;
      G4cout << "Auto-rotation set to : " << fRot_sens << G4endl;
    }
    else if (evnt->key() == Qt::Key_Minus) {
      fRot_sens = fRot_sens*0.7;
      G4cout << "Auto-rotation set to : " << fRot_sens << G4endl;
    }

    // Control Modifier OR Command on MAC
  }
  if ((fControlKeyPress)) {
    if (evnt->key() == Qt::Key_Plus) {
      fVP.SetZoomFactor(fVP.GetZoomFactor()*(1+fDeltaZoom)); 
      updateQWidget();
    }
    else if (evnt->key() == Qt::Key_Minus) {
      fVP.SetZoomFactor(fVP.GetZoomFactor()*(1-fDeltaZoom)); 
      updateQWidget();
    }
  }  
  
  fHoldKeyEvent = false;
}
  

void  G4OpenGLQtViewer::updateKeyModifierState(const Qt::KeyboardModifiers& modifier) {
  // Check Qt Versions for META Keys
    
  fNoKeyPress = true;
  fAltKeyPress = false;
  fShiftKeyPress = false;
  fControlKeyPress = false;
  
  if (modifier & Qt::AltModifier ) {
    fAltKeyPress = true;
    fNoKeyPress = false;
  }
  if (modifier & Qt::ShiftModifier ) {
    fShiftKeyPress = true;
    fNoKeyPress = false;
  }
  if (modifier & Qt::ControlModifier ) {
    fControlKeyPress = true;
    fNoKeyPress = false;
  }
}


/** Stop the video. Check all parameters and enable encoder button if all is ok.
 */
void G4OpenGLQtViewer::stopVideo() {

  // if encoder parameter is wrong, display parameters dialog and return
  if (!fMovieParametersDialog) {
    showMovieParametersDialog();
  }
  setRecordingStatus(STOP);

  if (fRecordFrameNumber >0) {
    // check parameters if they were modified (Re APPLY them...)
    if (!(fMovieParametersDialog->checkEncoderSwParameters())) {
      setRecordingStatus(BAD_ENCODER);
    }  else if (!(fMovieParametersDialog->checkSaveFileNameParameters())) {
      setRecordingStatus(BAD_OUTPUT);
    }
  } else {
    resetRecording();
    setRecordingInfos("No frame to encode.");
  }
}

/** Stop the video. Check all parameters and enable encoder button if all is ok.
 */
void G4OpenGLQtViewer::saveVideo() {

  // if encoder parameter is wrong, display parameters dialog and return
  if (!fMovieParametersDialog) {
    showMovieParametersDialog();
  }

  fMovieParametersDialog->checkEncoderSwParameters();
  fMovieParametersDialog->checkSaveFileNameParameters();
  
  if (fRecordingStep == STOP) {
    setRecordingStatus(SAVE);
    generateMpegEncoderParameters();
    encodeVideo();
  }
}


/** Start/Pause the video..
 */
void G4OpenGLQtViewer::startPauseVideo() {
   
  // first time, if temp parameter is wrong, display parameters dialog and return

  if ( fRecordingStep == WAIT) {
    if ( fRecordFrameNumber == 0) {
      if (getTempFolderPath() == "") { // BAD_OUTPUT
        showMovieParametersDialog();
        setRecordingInfos("You should specified the temp folder in order to make movie");
        return;
      } else  {
        // remove temp folder if it was create
        QString tmp = removeTempFolder();
        if (tmp !="") {
          setRecordingInfos(tmp);
          return;
        }
        tmp = createTempFolder();
        if (tmp != "") {
          setRecordingInfos("Can't create temp folder."+tmp);
          return;
        }
      }
    }
  }
  if (fRecordingStep == WAIT) {
    setRecordingStatus(START); 
  } else if (fRecordingStep == START) {
    setRecordingStatus(PAUSE);
  } else if (fRecordingStep == PAUSE) {
    setRecordingStatus(CONTINUE);
  } else if (fRecordingStep == CONTINUE) {
    setRecordingStatus(PAUSE);
  }
}

void G4OpenGLQtViewer::setRecordingStatus(RECORDING_STEP step) {

  fRecordingStep = step;
  displayRecordingStatus();
}


void G4OpenGLQtViewer::displayRecordingStatus() {
  
  QString txtStatus = "";
  if (fRecordingStep == WAIT) {
    txtStatus  = "Waiting to start...";
    fRecordFrameNumber = 0; // reset the frame number
  } else if (fRecordingStep == START) {
    txtStatus  = "Start Recording...";
  } else if (fRecordingStep == PAUSE) {
    txtStatus  = "Pause Recording...";
  } else if (fRecordingStep == CONTINUE) {
    txtStatus  = "Continue Recording...";
  } else if (fRecordingStep == STOP) {
    txtStatus  = "Stop Recording...";
  } else if (fRecordingStep == READY_TO_ENCODE) {
    txtStatus  = "Ready to Encode...";
  } else if (fRecordingStep == ENCODING) {
    txtStatus  = "Encoding...";
  } else if (fRecordingStep == FAILED) {
    txtStatus  = "Failed to encode...";
  } else if ((fRecordingStep == BAD_ENCODER)
             || (fRecordingStep == BAD_OUTPUT)
             || (fRecordingStep == BAD_TMP)) {
    txtStatus  = "Correct above errors first";
  } else if (fRecordingStep == SUCCESS) {
    txtStatus  = "File encoded successfully";
  } else {
  }

  if (fMovieParametersDialog) {
    fMovieParametersDialog->setRecordingStatus(txtStatus);
  } else {
    G4cout << txtStatus.toStdString().c_str() << G4endl;
  }
  setRecordingInfos("");
}


void G4OpenGLQtViewer::setRecordingInfos(const QString& txt) {
  if (fMovieParametersDialog) {
    fMovieParametersDialog->setRecordingInfos(txt);
  } else {
    G4cout << txt.toStdString().c_str() << G4endl;
  }
}

/** Init the movie parameters. Temp dir and encoder path
 */
void G4OpenGLQtViewer::initMovieParameters() {
  //init encoder
  
  //look for encoderPath
  fProcess = new QProcess();
     
  QObject ::connect(fProcess,SIGNAL(finished ( int)),
                    this,SLOT(processLookForFinished()));
  fProcess->setReadChannelMode(QProcess::MergedChannels);
  fProcess->start ("which ppmtompeg");
  
}

/** @return encoder path or "" if it does not exist
 */
QString G4OpenGLQtViewer::getEncoderPath() {
  return fEncoderPath;
}
 

/**
 * set the new encoder path
 * @return "" if correct. The error otherwise
 */
QString G4OpenGLQtViewer::setEncoderPath(QString path) {
  if (path == "") {
    return "File does not exist";
  }

  path =  QDir::cleanPath(path);
  QFileInfo *f = new QFileInfo(path);
  if (!f->exists()) {
    return "File does not exist";
  } else if (f->isDir()) {
    return "This is a directory";
  } else if (!f->isExecutable()) {
    return "File exist but is not executable";
  } else if (!f->isFile()) {
    return "This is not a file";
  }
  fEncoderPath = path;

  if (fRecordingStep == BAD_ENCODER) {
    setRecordingStatus(STOP);
  } 
  return "";
}


bool G4OpenGLQtViewer::isRecording(){
  if ((fRecordingStep == START) || (fRecordingStep == CONTINUE)) {
    return true;
  }
  return false;
}

bool G4OpenGLQtViewer::isPaused(){
  if (fRecordingStep == PAUSE) {
    return true;
  }
  return false;
}

bool G4OpenGLQtViewer::isEncoding(){
  if (fRecordingStep == ENCODING) {
    return true;
  }
  return false;
}

bool G4OpenGLQtViewer::isWaiting(){
  if (fRecordingStep == WAIT) {
    return true;
  }
  return false;
}

bool G4OpenGLQtViewer::isStopped(){
  if (fRecordingStep == STOP) {
    return true;
  }
  return false;
}

bool G4OpenGLQtViewer::isFailed(){
  if (fRecordingStep == FAILED) {
    return true;
  }
  return false;
}

bool G4OpenGLQtViewer::isSuccess(){
  if (fRecordingStep == SUCCESS) {
    return true;
  }
  return false;
}

bool G4OpenGLQtViewer::isBadEncoder(){
  if (fRecordingStep == BAD_ENCODER) {
    return true;
  }
  return false;
}
bool G4OpenGLQtViewer::isBadTmp(){
  if (fRecordingStep == BAD_TMP) {
    return true;
  }
  return false;
}
bool G4OpenGLQtViewer::isBadOutput(){
  if (fRecordingStep == BAD_OUTPUT) {
    return true;
  }
  return false;
}

void G4OpenGLQtViewer::setBadEncoder(){
  fRecordingStep = BAD_ENCODER;
  displayRecordingStatus();
}
void G4OpenGLQtViewer::setBadTmp(){
  fRecordingStep = BAD_TMP;
  displayRecordingStatus();
}
void G4OpenGLQtViewer::setBadOutput(){
  fRecordingStep = BAD_OUTPUT;
  displayRecordingStatus();
}

void G4OpenGLQtViewer::setWaiting(){
  fRecordingStep = WAIT;
  displayRecordingStatus();
}


bool G4OpenGLQtViewer::isReadyToEncode(){
  if (fRecordingStep == READY_TO_ENCODE) {
    return true;
  }
  return false;
}

void G4OpenGLQtViewer::resetRecording() {
  setRecordingStatus(WAIT);
}

/**
 * set the temp folder path
 * @return "" if correct. The error otherwise
 */
QString G4OpenGLQtViewer::setTempFolderPath(QString path) {

  if (path == "") {
    return "Path does not exist";
  }
  path =  QDir::cleanPath(path);
  QFileInfo *d = new QFileInfo(path);
  if (!d->exists()) {
    return "Path does not exist";
  } else if (!d->isDir()) {
    return "This is not a directory";
  } else if (!d->isReadable()) {
    return path +" is read protected";
  } else if (!d->isWritable()) {
    return path +" is write protected";
  }
  
  if (fRecordingStep == BAD_TMP) {
    setRecordingStatus(WAIT); 
  }
  fTempFolderPath = path;
  return "";
}

/** @return the temp folder path or "" if it does not exist
 */
QString G4OpenGLQtViewer::getTempFolderPath() {
  return fTempFolderPath;
}
 
/**
 * set the save file name path
 * @return "" if correct. The error otherwise
 */
QString G4OpenGLQtViewer::setSaveFileName(QString path) {

  if (path == "") {
    return "Path does not exist";
  }
  
  QFileInfo *file = new QFileInfo(path);
  QDir dir = file->dir();
  path =  QDir::cleanPath(path);
  if (file->exists()) {
    return "File already exist, please choose a new one";
  } else if (!dir.exists()) {
    return "Dir does not exist";
  } else if (!dir.isReadable()) {
    return path +" is read protected";
  }
  
  if (fRecordingStep == BAD_OUTPUT) {
    setRecordingStatus(STOP); 
  }
  fSaveFileName = path;
  return "";
}

/** @return the save file path
 */
QString G4OpenGLQtViewer::getSaveFileName() {
  return fSaveFileName ;
}

/** Create a Qt_temp folder in the temp folder given
 * The temp folder will be like this /tmp/QtMovie_12-02-2008_12_12_58/
 * @return "" if success. Error message if not.
 */
QString G4OpenGLQtViewer::createTempFolder() {
  fMovieTempFolderPath = "";
  //check
  QString tmp = setTempFolderPath(fTempFolderPath);
  if (tmp != "") {
    return tmp;
  }
  QString sep = QString(QDir::separator());
  QString path = sep+"QtMovie_"+QDateTime::currentDateTime ().toString("dd-MM-yyyy_hh-mm-ss")+sep; 
  QDir *d = new QDir(QDir::cleanPath(fTempFolderPath));
  // check if it is already present
  if (d->exists(path)) {
    return "Folder "+path+" already exists.Please remove it first";
  }
  if (d->mkdir(fTempFolderPath+path)) {
    fMovieTempFolderPath = fTempFolderPath+path;
    return "";
  }
  return "Can't create "+fTempFolderPath+path;
}

/** Remove the Qt_temp folder in the temp folder
 */
QString G4OpenGLQtViewer::removeTempFolder() {
  // remove files in Qt_temp folder
  if (fMovieTempFolderPath == "") {
    return "";
  }
  QDir *d = new QDir(QDir::cleanPath(fMovieTempFolderPath));
  if (!d->exists()) {
    return "";  // already remove
  }

  d->setFilter( QDir::Files );
  QStringList subDirList = d->entryList();
  int res = true;
  QString error = "";
  for (QStringList::ConstIterator it = subDirList.begin() ;(it != subDirList.end()) ; it++) {
    const QString currentFile = *it;
    if (!d->remove(currentFile)) {
      res = false;
      QString file = fMovieTempFolderPath+currentFile;
      error +="Removing file failed : "+file;
    } else {
    }
  }
  if (res) {
    if (d->rmdir(fMovieTempFolderPath)) {
      fMovieTempFolderPath = "";
      return "";
    } else {
      return "Dir "+fMovieTempFolderPath+" should be empty, but could not remove it";
    }

  }
  return "Could not remove "+fMovieTempFolderPath+" because of the following errors :"+error;
}



bool G4OpenGLQtViewer::hasPendingEvents () {
  return ((QApplication*)G4Qt::getInstance ())->hasPendingEvents ();
}

bool G4OpenGLQtViewer::generateMpegEncoderParameters () {

  // save the parameter file
  FILE* fp;
  fp = fopen (QString(fMovieTempFolderPath+fParameterFileName).toStdString().c_str(), "w");

  if (fp == NULL) {
    setRecordingInfos("Generation of parameter file failed");
    return false;
  }

  fprintf (fp,"# Pattern affects speed, quality and compression. See the User's Guide\n");
  fprintf (fp,"# for more info.\n");
  fprintf (fp,"\n");
  fprintf (fp,"PATTERN		IBBPBBPBBPBBPBBP\n");
  fprintf (fp,"OUTPUT		%s\n",getSaveFileName().toStdString().c_str());
  fprintf (fp,"\n");
  fprintf (fp,"# You must specify the type of the input files.  The choices are:\n");
  fprintf (fp,"#    YUV, PPM, JMOVIE, Y, JPEG, PNM\n");
  fprintf (fp,"#	(must be upper case)\n");
  fprintf (fp,"#\n");
  fprintf (fp,"BASE_FILE_FORMAT	PPM\n");
  fprintf (fp,"\n");
  fprintf (fp,"#\n");
  fprintf (fp,"# if YUV format (or using parallel version), must provide width and height\n");
  fprintf (fp,"# YUV_SIZE	widthxheight\n");
  fprintf (fp,"# this option is ignored if BASE_FILE_FORMAT is not YUV and you're running\n");
  fprintf (fp,"# on just one machine\n");
  fprintf (fp,"#\n");
  fprintf (fp,"YUV_SIZE	352x240\n");
  fprintf (fp,"\n");
  fprintf (fp,"# If you are using YUV, there are different supported file formats.\n");
  fprintf (fp,"# EYUV or UCB are the same as previous versions of this encoder.\n");
  fprintf (fp,"# (All the Y's, then U's then V's, in 4:2:0 subsampling.)\n");
  fprintf (fp,"# Other formats, such as Abekas, Phillips, or a general format are\n");
  fprintf (fp,"# permissible, the general format is a string of Y's, U's, and V's\n");
  fprintf (fp,"# to specify the file order.\n");
  fprintf (fp,"\n");
  fprintf (fp,"INPUT_FORMAT UCB\n");
  fprintf (fp,"\n");
  fprintf (fp,"# the conversion statement\n");
  fprintf (fp,"#\n");
  fprintf (fp,"# Each occurrence of '*' will be replaced by the input file\n");
  fprintf (fp,"#\n");
  fprintf (fp,"# e.g., if you have a bunch of GIF files, then this might be:\n");
  fprintf (fp,"#	INPUT_CONVERT	giftoppm *\n");
  fprintf (fp,"#\n");
  fprintf (fp,"# e.g., if you have a bunch of files like a.Y a.U a.V, etc., then:\n");
  fprintf (fp,"#	INPUT_CONVERT	cat *.Y *.U *.V\n");
  fprintf (fp,"#\n");
  fprintf (fp,"# e.g., if you are grabbing from laser disc you might have something like\n");
  fprintf (fp,"#	INPUT_CONVERT	goto frame *; grabppm\n");
  fprintf (fp,"# 'INPUT_CONVERT *' means the files are already in the base file format\n");
  fprintf (fp,"#\n");
  fprintf (fp,"INPUT_CONVERT	* \n");
  fprintf (fp,"\n");
  fprintf (fp,"# number of frames in a GOP.\n");
  fprintf (fp,"#\n");
  fprintf (fp,"# since each GOP must have at least one I-frame, the encoder will find the\n");
  fprintf (fp,"# the first I-frame after GOP_SIZE frames to start the next GOP\n");
  fprintf (fp,"#\n");
  fprintf (fp,"# later, will add more flexible GOP signalling\n");
  fprintf (fp,"#\n");
  fprintf (fp,"GOP_SIZE	16\n");
  fprintf (fp,"\n");
  fprintf (fp,"# number of slices in a frame\n");
  fprintf (fp,"#\n");
  fprintf (fp,"# 1 is a good number.  another possibility is the number of macroblock rows\n");
  fprintf (fp,"# (which is the height divided by 16)\n");
  fprintf (fp,"#\n");
  fprintf (fp,"SLICES_PER_FRAME	1\n");
  fprintf (fp,"\n");
  fprintf (fp,"# directory to get all input files from (makes this file easier to read)\n");
  fprintf (fp,"INPUT_DIR	%s\n",fMovieTempFolderPath.toStdString().c_str());
  fprintf (fp,"\n");
  fprintf (fp,"# There are a bunch of ways to specify the input files.\n");
  fprintf (fp,"# from a simple one-per-line listing, to the following \n");
  fprintf (fp,"# way of numbering them.  See the manual for more information.\n");
  fprintf (fp,"INPUT\n");
  fprintf (fp,"# '*' is replaced by the numbers 01, 02, 03, 04\n");
  fprintf (fp,"# if I instead do [01-11], it would be 01, 02, ..., 09, 10, 11\n");
  fprintf (fp,"# if I instead do [1-11], it would be 1, 2, 3, ..., 9, 10, 11\n");
  fprintf (fp,"# if I instead do [1-11+3], it would be 1, 4, 7, 10\n");
  fprintf (fp,"# the program assumes none of your input files has a name ending in ']'\n");
  fprintf (fp,"# if you do, too bad!!!\n");
  fprintf (fp,"#\n");
  fprintf (fp,"#\n");
  fprintf (fp,"Test*.ppm	[0-%d]\n",fRecordFrameNumber-1);
  fprintf (fp,"# can have more files here if you want...there is no limit on the number\n");
  fprintf (fp,"# of files\n");
  fprintf (fp,"END_INPUT\n");
  fprintf (fp,"\n");
  fprintf (fp,"\n");
  fprintf (fp,"\n");
  fprintf (fp,"# Many of the remaining options have to do with the motion search and qscale\n");
  fprintf (fp,"\n");
  fprintf (fp,"# FULL or HALF -- must be upper case\n");
  fprintf (fp,"# Should be FULL for computer generated images\n");
  fprintf (fp,"PIXEL		FULL\n");
  fprintf (fp,"\n");
  fprintf (fp,"# means +/- this many pixels for both P and B frame searches\n");
  fprintf (fp,"# specify two numbers if you wish to serc different ranges in the two.\n");
  fprintf (fp,"RANGE		10\n");
  fprintf (fp,"\n");
  fprintf (fp,"# The two search algorithm parameters below mostly affect speed,\n");
  fprintf (fp,"# with some affect on compression and almost none on quality.\n");
  fprintf (fp,"\n");
  fprintf (fp,"# this must be one of {EXHAUSTIVE, SUBSAMPLE, LOGARITHMIC}\n");
  fprintf (fp,"PSEARCH_ALG	LOGARITHMIC\n");
  fprintf (fp,"\n");
  fprintf (fp,"# this must be one of {SIMPLE, CROSS2, EXHAUSTIVE}\n");
  fprintf (fp,"#\n");
  fprintf (fp,"# note that EXHAUSTIVE is really, really, really slow\n");
  fprintf (fp,"#\n");
  fprintf (fp,"BSEARCH_ALG	SIMPLE\n");
  fprintf (fp,"\n");
  fprintf (fp,"#\n");
  fprintf (fp,"# these specify the q-scale for I, P, and B frames\n");
  fprintf (fp,"# (values must be between 1 and 31)\n");
  fprintf (fp,"# These are the Qscale values for the entire frame in variable bit-rate\n");
  fprintf (fp,"# mode, and starting points (but not important) for constant bit rate\n");
  fprintf (fp,"#\n");
  fprintf (fp,"\n");
  fprintf (fp,"# Qscale (Quantization scale) affects quality and compression,\n");
  fprintf (fp,"# but has very little effect on speed.\n");
  fprintf (fp,"\n");
  fprintf (fp,"IQSCALE		4\n");
  fprintf (fp,"PQSCALE		5\n");
  fprintf (fp,"BQSCALE		12\n");
  fprintf (fp,"\n");
  fprintf (fp,"# this must be ORIGINAL or DECODED\n");
  fprintf (fp,"REFERENCE_FRAME	ORIGINAL\n");
  fprintf (fp,"\n");
  fprintf (fp,"# for parallel parameters see parallel.param in the exmaples subdirectory\n");
  fprintf (fp,"\n");
  fprintf (fp,"# if you want constant bit-rate mode, specify it as follows (number is bits/sec):\n");
  fprintf (fp,"#BIT_RATE  1000000\n");
  fprintf (fp,"\n");
  fprintf (fp,"# To specify the buffer size (327680 is default, measused in bits, for 16bit words)\n");
  fprintf (fp,"BUFFER_SIZE 327680\n");
  fprintf (fp,"\n");
  fprintf (fp,"# The frame rate is the number of frames/second (legal values:\n");
  fprintf (fp,"# 23.976, 24, 25, 29.97, 30, 50 ,59.94, 60\n");
  fprintf (fp,"FRAME_RATE 30\n");
  fprintf (fp,"\n");
  fprintf (fp,"# There are many more options, see the users manual for examples....\n");
  fprintf (fp,"# ASPECT_RATIO, USER_DATA, GAMMA, IQTABLE, etc.\n");
  fprintf (fp,"\n");
  fprintf (fp,"\n");
  fclose (fp);

  setRecordingInfos("Parameter file "+fParameterFileName+" generated in "+fMovieTempFolderPath);
  setRecordingStatus(READY_TO_ENCODE);
  return true;
}

void G4OpenGLQtViewer::encodeVideo()
{
  if ((getEncoderPath() != "") && (getSaveFileName() != "")) {
    setRecordingStatus(ENCODING);
    
    fProcess = new QProcess();
#if QT_VERSION > 0x040100
    QObject ::connect(fProcess,SIGNAL(finished ( int,QProcess::ExitStatus)),
                      this,SLOT(processEncodeFinished()));
    QObject ::connect(fProcess,SIGNAL(readyReadStandardOutput ()),
                      this,SLOT(processEncodeStdout()));
#else
    QObject ::connect(fProcess,SIGNAL(finished ( int)),
                      this,SLOT(processEncodeFinished()));
    QObject ::connect(fProcess,SIGNAL(readyReadStandardOutput ()),
                      this,SLOT(processEncodeStdout()));
#endif
    fProcess->setReadChannelMode(QProcess::MergedChannels);
    fProcess->start (fEncoderPath, QStringList(fMovieTempFolderPath+fParameterFileName));
  }
}


// FIXME : does not work on Qt3
void G4OpenGLQtViewer::processEncodeStdout()
{
  QString tmp = fProcess->readAllStandardOutput ().data();
  int start = tmp.lastIndexOf("ESTIMATED TIME");
  tmp = tmp.mid(start,tmp.indexOf("\n",start)-start);
  setRecordingInfos(tmp);
}


void G4OpenGLQtViewer::processEncodeFinished()
{

  QString txt = "";
  txt = getProcessErrorMsg();
  if (txt == "") {
    setRecordingStatus(SUCCESS);
  } else {
    setRecordingStatus(FAILED);
  }
  //  setRecordingInfos(txt+removeTempFolder());
}


void G4OpenGLQtViewer::processLookForFinished() 
{

  QString txt = getProcessErrorMsg();
  if (txt != "") {
    fEncoderPath = "";
  } else {
    fEncoderPath = QString(fProcess->readAllStandardOutput ().data()).trimmed();
    // if not found, return "not found"
    if (fEncoderPath.contains(" ")) {
      fEncoderPath = "";
    } else if (!fEncoderPath.contains("ppmtompeg")) {
      fEncoderPath = "";
    }
    setEncoderPath(fEncoderPath);
  }
  // init temp folder
  setTempFolderPath(QDir::temp ().absolutePath ());
}


QString G4OpenGLQtViewer::getProcessErrorMsg()
{
  QString txt = "";
  if (fProcess->exitCode() != 0) {
    switch (fProcess->error()) {
    case QProcess::FailedToStart:
      txt = "The process failed to start. Either the invoked program is missing, or you may have insufficient permissions to invoke the program.\n";
      break;
    case QProcess::Crashed:
      txt = "The process crashed some time after starting successfully.\n";
      break;
    case QProcess::Timedout:
      txt = "The last waitFor...() function timed out. The state of QProcess is unchanged, and you can try calling waitFor...() again.\n";
      break;
    case QProcess::WriteError:
      txt = "An error occurred when attempting to write to the process. For example, the process may not be running, or it may have closed its input channel.\n";
      break;
    case QProcess::ReadError:
      txt = "An error occurred when attempting to read from the process. For example, the process may not be running.\n";
      break;
    case QProcess::UnknownError:
      txt = "An unknown error occurred. This is the default return value of error().\n";
      break;
    }
  }
  return txt;
}




QWidget *G4OpenGLQtViewer::getParentWidget() 
{
  // launch Qt if not
  G4Qt* interactorManager = G4Qt::getInstance ();
  // G4UImanager* UI = 
  // G4UImanager::GetUIpointer();
  
  bool found = false;
  
  // create window
  if (((QApplication*)interactorManager->GetMainInteractor())) {
    // look for the main window
    QWidgetList wl = QApplication::allWidgets();
    QWidget *widget = NULL;
    for (int i=0; i < wl.size(); i++) {
      widget = wl.at(i);
      if ((found== false) && (widget->inherits("QMainWindow"))) {
        fGLWindow = new QDialog(widget,Qt::WindowTitleHint | Qt::WindowSystemMenuHint | Qt::WindowMinMaxButtonsHint);
        found = true;
      }
    }
    
    if (found==false) {
#ifdef G4DEBUG_VIS_OGL
      printf("G4OpenGLQtViewer::CreateMainWindow case Qapp exist, but not found\n");
#endif
      fGLWindow = new QDialog();
    }
  } else {
#ifdef G4DEBUG_VIS_OGL
    printf("G4OpenGLQtViewer::CreateMainWindow case Qapp exist\n");
#endif
    fGLWindow = new QDialog();
#ifdef G4DEBUG_VIS_OGL
    printf("G4OpenGLQtViewer::GetParentWidget fGLWindow\n");
#endif
  }
  if (found) {
    return fGLWindow;
  } else {
    return NULL;
  }
}


void G4OpenGLQtViewer::initSceneTreeComponent(){

  fSceneTreeWidget = new QWidget();
  fSceneTreeWidget->setLayout (new QVBoxLayout());

  if (fUISceneTreeComponentsTBWidget != NULL) {
    fUISceneTreeComponentsTBWidget->addTab(fSceneTreeWidget,QString(GetName().data()));
  }

  QLayout* layoutSceneTreeComponentsTBWidget = fSceneTreeWidget->layout();


  fSceneTreeComponentTreeWidget = new QTreeWidget();
  fSceneTreeComponentTreeWidget->setSelectionMode(QAbstractItemView::SingleSelection);
  fSceneTreeComponentTreeWidget->setHeaderLabel ("Scene tree : "+QString(GetName().data()));
  fSceneTreeComponentTreeWidget->setColumnHidden (1,true);  // copy number
  fSceneTreeComponentTreeWidget->setColumnHidden (2,true);  // PO index
  fSceneTreeComponentTreeWidget->setColumnHidden (3,true);  // Informations
  //   data(0) : POindex
  //   data(1) : copy number
  //   data(2) : g4color

  layoutSceneTreeComponentsTBWidget->addWidget(fSceneTreeComponentTreeWidget);

  connect(fSceneTreeComponentTreeWidget,SIGNAL(itemChanged(QTreeWidgetItem*, int)),SLOT(sceneTreeComponentItemChanged(QTreeWidgetItem*, int)));
  connect(fSceneTreeComponentTreeWidget,SIGNAL(itemSelectionChanged ()),SLOT(sceneTreeComponentSelected()));
  connect(fSceneTreeComponentTreeWidget,SIGNAL(itemDoubleClicked ( QTreeWidgetItem*, int)),SLOT(changeColorAndTransparency( QTreeWidgetItem*, int)));


  // Depth slider
  QLabel *depth = new QLabel();
  depth->setText("Depth :");

  QWidget* depthWidget = new QWidget();
  QGroupBox *groupBox = new QGroupBox(tr("Touchables slider"),depthWidget);
  QHBoxLayout *groupBoxLayout = new QHBoxLayout(groupBox);

  QLabel *zero = new QLabel();
  zero->setText("Show all");          
  QLabel *one = new QLabel();
  one->setText("Hide all");
  fSceneTreeDepthSlider = new QSlider ( Qt::Horizontal, groupBox);
  fSceneTreeDepthSlider->setMaximum (1000);
  fSceneTreeDepthSlider->setMinimum (0);
  fSceneTreeDepthSlider->setTickPosition(QSlider::TicksAbove);
  groupBoxLayout->addWidget(zero);
  groupBoxLayout->addWidget(fSceneTreeDepthSlider);
  groupBoxLayout->addWidget(one);

  groupBox->setLayout(groupBoxLayout);
  layoutSceneTreeComponentsTBWidget->addWidget(groupBox);

  connect( fSceneTreeDepthSlider, SIGNAL( valueChanged(int) ), this, SLOT( changeDepthInSceneTree(int) ) );

  // Search line
  QWidget *helpWidget = new QWidget();
  QHBoxLayout *helpLayout = new QHBoxLayout();
  QPushButton* select = new QPushButton("select item(s)");
  fHelpLine = new QLineEdit();
  helpLayout->addWidget(new QLabel("Search :",helpWidget));
  helpLayout->addWidget(fHelpLine);
  helpLayout->addWidget(select);
  helpWidget->setLayout(helpLayout);
  layoutSceneTreeComponentsTBWidget->addWidget(helpWidget);

  connect( fHelpLine, SIGNAL( returnPressed () ), this, SLOT(changeSearchSelection()));
  connect( select, SIGNAL( clicked () ), this, SLOT(changeSearchSelection()));


  fTreeItemModels.clear();

  fPVRootNodeCreate = false;

  fMaxPOindexInserted = -1;
}


// set the component to check/unchecked, also go into its child
// and set the same status to all his childs
void G4OpenGLQtViewer::setCheckComponent(QTreeWidgetItem* item,bool check)
{
  if (item != NULL) {
    if (check) {
      item->setCheckState(0,Qt::Checked);
    } else {
      item->setCheckState(0,Qt::Unchecked);
    }
      updatePositivePoIndexSceneTreeWidgetQuickMap(item->data(0,Qt::UserRole).toInt(),item);
    int nChildCount = item->childCount();
    for (int i = 0; i < nChildCount; i++) {
      setCheckComponent(item->child(i),check);
    }
  }
}


void G4OpenGLQtViewer::DrawText(const G4Text& g4text)
{
  if (isGl2psWriting()) {

    G4OpenGLViewer::DrawText(g4text);

  } else {

    if (!fWindow) return;

    G4VSceneHandler::MarkerSizeType sizeType;
    G4double size = fSceneHandler.GetMarkerSize(g4text,sizeType);

    QFont font = QFont();
    font.setPointSizeF(size);

    G4Point3D position = g4text.GetPosition();

    const G4String& textString = g4text.GetText();
    const char* textCString = textString.c_str();
  
    // Calculate move for centre and right adjustment
    QFontMetrics* f = new QFontMetrics (font);
    G4double span = f->width(textCString);
    G4double xmove = 0., ymove = 0.;
    switch (g4text.GetLayout()) {
    case G4Text::left: break;
    case G4Text::centre: xmove -= span / 2.; break;
    case G4Text::right: xmove -= span;
    }
    
    //Add offsets
    xmove += g4text.GetXOffset();
    ymove += g4text.GetYOffset();

    fWindow->renderText
      ((position.x()+(2*xmove)/getWinWidth()),
       (position.y()+(2*ymove)/getWinHeight()),
       position.z(),
       textCString,
       font);

  }
}


void G4OpenGLQtViewer::ResetView () {
  G4OpenGLViewer::ResetView();
  fDeltaDepth = 0.01;
  fDeltaZoom = 0.05;
}




void G4OpenGLQtViewer::addPVSceneTreeElement(const G4String& model, G4PhysicalVolumeModel* pPVModel, int currentPOIndex) {

  const QString& modelShortName = getModelShortName(model);

  if (modelShortName == "") {
    return ;
  }
  // try to init it
  if (fSceneTreeComponentTreeWidget == NULL) {
    initSceneTreeComponent();
  }

  // if no UI
  if (fSceneTreeComponentTreeWidget == NULL) {
    return;
  }

  fSceneTreeComponentTreeWidget->blockSignals(true);

  // Create the "volume" node if not
  //  if (fSceneTreeComponentTreeWidget->topLevelItemCount () == 0) {
  if (!fPVRootNodeCreate) {
    const G4VisAttributes* visAttr = GetApplicableVisAttributes(pPVModel->GetFullPVPath().at(0).GetPhysicalVolume()->GetLogicalVolume()->GetVisAttributes());
    const G4Colour& color = visAttr->GetColour();
    
    fModelShortNameItem = createTreeWidgetItem(pPVModel->GetFullPVPath(),
                                               modelShortName,
                                               0, // currentPVCopyNb
                                               -1, // currentPVPOIndex
                                               "",
                                               Qt::Checked,
                                               NULL,
                                               color);
    fPVRootNodeCreate = true;
  }

  bool added = parseAndInsertInSceneTree(fModelShortNameItem,pPVModel,0,modelShortName,0,currentPOIndex);
  if (!added) {
  }
  
  fSceneTreeComponentTreeWidget->blockSignals(false);
  
}


/**
   if treeNode is NULL, then add this treeNode to the TreeWidget
   @return the inserted item
*/
QTreeWidgetItem* G4OpenGLQtViewer::createTreeWidgetItem(
 const PVPath& fullPath
 ,const QString& name
 ,int copyNb
 ,int POIndex
 ,const QString& logicalName
 ,Qt::CheckState state
 ,QTreeWidgetItem * parentTreeNode
 ,const G4Colour& color
) {

  // Set depth
  if (fullPath.size() > fSceneTreeDepth) {
    fSceneTreeDepth = fullPath.size();
    // Change slider value
    if (fSceneTreeDepthSlider) {
      fSceneTreeDepthSlider->setTickInterval(1000/(fSceneTreeDepth+1));
    }
  }
  QTreeWidgetItem * newItem = NULL;
  if (parentTreeNode == NULL) {
    newItem = new QTreeWidgetItem(fSceneTreeComponentTreeWidget);
  } else {
    newItem = new QTreeWidgetItem(parentTreeNode);
  }


  newItem->setText(0,name);
  newItem->setData(1,Qt::UserRole,copyNb);
  newItem->setText(2,QString::number(POIndex));
  newItem->setData(0, Qt::UserRole, POIndex);
  newItem->setText(3,logicalName);
  newItem->setFlags(newItem->flags()|Qt::ItemIsUserCheckable);
  newItem->setCheckState(0,state);
  updatePositivePoIndexSceneTreeWidgetQuickMap(POIndex,newItem);

  changeQColorForTreeWidgetItem(newItem,QColor((int)(color.GetRed()*255),
                                               (int)(color.GetGreen()*255),
                                               (int)(color.GetBlue()*255),
                                               (int)(color.GetAlpha()*255)));

  // If invisible
  if ((state == Qt::Unchecked) && (POIndex == -1)) {
    newItem->setForeground (0, QBrush( Qt::gray) );
    
    // Set a tootip
    newItem->setToolTip (0,QString(
                                   "This node exists in the geometry but has not been\n")+
                         "drawn, perhaps because it has been set invisible. It \n"+
                         "cannot be made visible with a click on the button.\n"+
                         "To see it, change the visibility, for example, with \n"+
                         "/vis/geometry/set/visibility " + logicalName + " 0 true\n"+
                         "and rebuild the view with /vis/viewer/rebuild.\n"+
                         "Click here will only show/hide all child components");
  }

  // special case: if alpha=0, it is a totally transparent objet,
  // then, do not redraw it
  if (color.GetAlpha() == 0) {
    state = Qt::Unchecked;
    newItem->setCheckState(0,state);
    updatePositivePoIndexSceneTreeWidgetQuickMap(POIndex,newItem);
  }

  fTreeItemModels.insert(std::pair <int, PVPath > (POIndex,fullPath) );

  // Check last status of this item and change if necessary 
  // open/close/hidden/visible/selected
  changeOpenCloseVisibleHiddenSelectedColorSceneTreeElement(newItem);
  return newItem;
}


//
//   Recursive function.
//   Try to insert the given item :
//   - If not present and last item of the path: insert it and mark it CHECK
//   - If not present and NOT last item of the path: insert it and mark it UNCHECKED
//   - If already present and name/PO/Transformation identical, then it is a transparent
//     object : Change the PO number and transparency
//   - If already present and PO different, then it is an unvisible item : Have to
//     set it visible
//   - else : Create a new element
//   @return true if inserted, false if already present
//
bool G4OpenGLQtViewer::parseAndInsertInSceneTree(
 QTreeWidgetItem * parentItem
 ,G4PhysicalVolumeModel* pPVModel
 ,unsigned int fullPathIndex
 ,const QString& parentRoot
 ,unsigned int currentIndexInTreeSceneHandler
 ,int currentPVPOIndex
) {

  if (parentItem == NULL) {
    return false;
  }

  const PVPath& fullPath = pPVModel->GetFullPVPath();

  std::ostringstream oss;
  oss << fullPath.at(fullPathIndex).GetCopyNo();
  std::string currentPVName = G4String(fullPath.at(fullPathIndex).GetPhysicalVolume()->GetName()+" ["+oss.str()+"]").data();

  int currentPVCopyNb = fullPath.at(fullPathIndex).GetCopyNo();

  const G4VisAttributes* visAttr = GetApplicableVisAttributes(fullPath.at(fullPathIndex).GetPhysicalVolume()->GetLogicalVolume()->GetVisAttributes());
  const G4Colour& color = visAttr->GetColour();

  // look in all children in order to get if their is already a
  // child corresponding:
  // - if so, go into this child
  // - if not : create it as invisible

  // Realy quick check if the POindex is already there
  QTreeWidgetItem* subItem = NULL;
  QList<QTreeWidgetItem *> parentItemList;
  
    
  // first of all, very quick check if it was not the same as last one

  // Check only if it is a transparent object
  // If it is the last item and it is not transparent -> nothing to look for, 
  // simply add it.
  if ((currentIndexInTreeSceneHandler == (fullPath.size()-1)) && ((color.GetAlpha() == 1.))) {
  } else {
    QString lookForString = QString(currentPVName.c_str());
    for (int i = 0;i < parentItem->childCount(); i++ ) {
      if (parentItem->child(i)->text(0) == lookForString) {
        parentItemList.push_back(parentItem->child(i));
      }
    }
  }
   
  for (int i = 0; i < parentItemList.size(); ++i) {
    const std::string& parentItemName = parentItemList.at(i)->text(0).toStdString();
    int parentItemCopyNb = parentItemList.at(i)->data(1,Qt::UserRole).toInt();
    int parentItemPOIndex = parentItemList.at(i)->data(0,Qt::UserRole).toInt();
        
    // if already inside
    // -> return true
    // special case, do not have to deal with hierarchy except for PhysicalVolume

    
    /* Physical Volume AND copy number equal AND  name equal */
    if (((parentRoot == fTouchableVolumes) &&  (currentPVCopyNb == parentItemCopyNb)
         && (currentPVName == parentItemName))        ||
        /* NOT a Physical Volume AND  copy number equal */
        ((parentRoot != fTouchableVolumes) && (currentPVCopyNb == parentItemCopyNb) 
         /*AND  name equal AND  PO index equal*/
         && (currentPVName == parentItemName) && (currentPVPOIndex == parentItemPOIndex) )) {
      
      // then check for the Transform3D
      bool sameTransform = true;
      if (parentItemPOIndex >= 0) {
        const PVPath& fullPathTmp = fTreeItemModels[parentItemPOIndex];
        if (fullPathTmp.size() > 0) {
          if (fullPathTmp.at(fullPathTmp.size()-1).GetTransform () == pPVModel->GetTransformation ()) {
            sameTransform = true;
          } else {
            sameTransform = false;
          }
        }
      }

      // Same transformation, then try to change the PO index
      if (sameTransform == true) {
        // already exist in the tree, is it a transparent object ?
        // If so, then have to change the PO index ONLY if it is the last
        // and then change the state ONLY if POIndex has change
        // If not, then go deaper
        
        // last element
        if (currentIndexInTreeSceneHandler == (fullPath.size()-1)) {
          
          parentItemList.at(i)->setText(2,QString::number(currentPVPOIndex));
          parentItemList.at(i)->setData(0, Qt::UserRole,currentPVPOIndex);
          
          fTreeItemModels.insert(std::pair <int, PVPath >(currentPVPOIndex,fullPath) );
          
          // Then remove tooltip and special font
          QFont f = QFont();
          parentItemList.at(i)->setFont (0,f);
          
          // set foreground
          parentItemList.at(i)->setForeground (0,QBrush());
          
          // Set a tootip
          parentItemList.at(i)->setToolTip (0,"");
          
          changeQColorForTreeWidgetItem(parentItemList.at(i),QColor((int)(color.GetRed()*255),
                                                                    (int)(color.GetGreen()*255),
                                                                    (int)(color.GetBlue()*255),
                                                                    (int)(color.GetAlpha()*255)));
          
          // set check only if there is something to display
          if (color.GetAlpha() > 0) {
            parentItemList.at(i)->setCheckState(0,Qt::Checked);
            updatePositivePoIndexSceneTreeWidgetQuickMap(currentPVPOIndex,parentItemList.at(i));
          }
          return false;
        } else {
          subItem = parentItemList.at(i);        
        }
        
        // Exists but not the end of path, then forget get it
      } else if (currentIndexInTreeSceneHandler < (fullPath.size()-1))  {
        subItem = parentItemList.at(i);
      }
    }

  } // end for
  
  // if it the last, then add it and set it checked
  if (currentIndexInTreeSceneHandler == (fullPath.size()-1)) {
    /* subItem =*/ createTreeWidgetItem(fullPath,
                                   QString(currentPVName.c_str()),
                                   currentPVCopyNb,
                                   currentPVPOIndex,
                                   QString(fullPath.at(fullPathIndex).GetPhysicalVolume()->GetLogicalVolume()->GetName().data()),
                                   Qt::Checked,
                                   parentItem,
                                   color);

    if (currentPVPOIndex > fMaxPOindexInserted) {
      fMaxPOindexInserted = currentPVPOIndex;
    }

  } else {
    
    // if no child found, then this child is create and marked as invisible, then go inside
    if (subItem == NULL) {
      
      if (currentIndexInTreeSceneHandler < (fullPath.size()-1))  {
        subItem = createTreeWidgetItem(fullPath,
                                       QString(currentPVName.c_str()),
                                       currentPVCopyNb,
                                       -1,
                                       QString(fullPath.at(fullPathIndex).GetPhysicalVolume()->GetLogicalVolume()->GetName().data()),
                                       Qt::Unchecked,
                                       parentItem,
                                       color);
      }
    }
    
    return parseAndInsertInSceneTree(subItem,pPVModel,fullPathIndex+1,parentRoot,currentIndexInTreeSceneHandler+1,currentPVPOIndex);
  }
  return true;
}


void G4OpenGLQtViewer::changeOpenCloseVisibleHiddenSelectedColorSceneTreeElement(
 QTreeWidgetItem* subItem
)
{
  // Check if object with the same POIndex is the same in old tree
  QTreeWidgetItem* oldItem = NULL;

  QTreeWidgetItem* foundItem = getOldTreeWidgetItem(subItem->data(0,Qt::UserRole).toInt());

  if (foundItem != NULL) {
    if (isSameSceneTreeElement(foundItem,subItem)) {
      oldItem = foundItem;
    }
  }
  if (foundItem == NULL) {  // PO should have change, parse all

    // POindex > 0
    std::map <int, QTreeWidgetItem*>::const_iterator i;
    i = fOldPositivePoIndexSceneTreeWidgetQuickMap.begin();
    while (i != fOldPositivePoIndexSceneTreeWidgetQuickMap.end()) {
      if (isSameSceneTreeElement(i->second,subItem)) {
        oldItem = i->second;
        i = fOldPositivePoIndexSceneTreeWidgetQuickMap.end();
      } else {
        i++;
      }
    }
    // POindex == 0 ?
    if (oldItem == NULL) {
      unsigned int a = 0;
      while (a < fOldNullPoIndexSceneTreeWidgetQuickVector.size()) {
        if (isSameSceneTreeElement(fOldNullPoIndexSceneTreeWidgetQuickVector[a],subItem)) {
          oldItem = fOldNullPoIndexSceneTreeWidgetQuickVector[a];
          a = fOldNullPoIndexSceneTreeWidgetQuickVector.size();
        } else {
          a++;
        }
      }
    }
  }

  // if found : retore old state
  if (oldItem != NULL) {
    subItem->setFlags(oldItem->flags());   // flags
    subItem->setCheckState(0,oldItem->checkState(0)); // check state
    subItem->setSelected(oldItem->isSelected());  // selected
    subItem->setExpanded(oldItem->isExpanded ());  // expand

    // change color
    // when we call this function, the color in the item is the one of vis Attr

    std::map <int, QTreeWidgetItem* >::iterator it;
    
    // getOldPO
    int oldPOIndex = oldItem->data(0,Qt::UserRole).toInt();
    it = fOldPositivePoIndexSceneTreeWidgetQuickMap.find(oldPOIndex);
    QColor color;

    // get old Vis Attr Color
    std::map <int, QColor >::iterator itVis;
    itVis = fOldVisAttrColorMap.find(oldPOIndex);

    QColor oldVisAttrColor;
    const QColor& newVisAttrColor = subItem->data(2,Qt::UserRole).value<QColor>();

    bool visAttrChange = false;
    // if old vis attr color found
    if (itVis != fOldVisAttrColorMap.end()) {
      oldVisAttrColor = itVis->second;
      if (oldVisAttrColor != newVisAttrColor) {
        visAttrChange = true;
      }
    } else {
      visAttrChange = true;
    }
    
    if (visAttrChange) {
      fOldVisAttrColorMap.insert(std::pair <int, QColor > (subItem->data(0,Qt::UserRole).toInt(),newVisAttrColor) );

    } else { // if no changes, get old PO value
      // if old PO found
      if (it != fOldPositivePoIndexSceneTreeWidgetQuickMap.end()) {
        color = (it->second)->data(2,Qt::UserRole).value<QColor>(); 
      } else {
        color = oldItem->data(2,Qt::UserRole).value<QColor>(); 
      }
#ifdef G4DEBUG_VIS_OGL
      printf("====color name:%s\n",color.name().toStdString().c_str());
#endif
      changeQColorForTreeWidgetItem(subItem,color);
    }
  }

  return;
}



// Check if both items are identical.
// For that, check name, copy number, transformation
// special case for "non Touchables", do not check the PO index, check only the name
bool G4OpenGLQtViewer::isSameSceneTreeElement(
 QTreeWidgetItem* parentOldItem
 ,QTreeWidgetItem* parentNewItem
) {

  int newPO = -1;
  int oldPO = -1;
  
  int newCpNumber = -1;
  int oldCpNumber = -1;

  bool firstWhile = true;

  while ((parentOldItem != NULL) && (parentNewItem != NULL)) {

    // check transform, optimize getting data(..,..) that consume lot of time
    if (!firstWhile) {
      oldPO = parentOldItem->data(0,Qt::UserRole).toInt();
      newPO = parentNewItem->data(0,Qt::UserRole).toInt();
    }
    firstWhile = false;
    
    if ((oldPO >= 0) &&
        (newPO >= 0)) {
      const PVPath& oldFullPath = fOldTreeItemModels[oldPO];
      const PVPath& newFullPath = fTreeItemModels[newPO];
      if ((oldFullPath.size() > 0) &&
          (newFullPath.size() > 0)) {
        if (oldFullPath.size() != newFullPath.size()) {
          return false;
        }
        if (oldFullPath.at(oldFullPath.size()-1).GetTransform () == newFullPath.at(newFullPath.size()-1).GetTransform ()) {
          newCpNumber = newFullPath.at(newFullPath.size()-1).GetCopyNo();
          oldCpNumber = oldFullPath.at(oldFullPath.size()-1).GetCopyNo();
          // ok
        } else {
          return false;
        }
      }
    }
    
    // Check copy Number
    if (oldCpNumber == -1) {
      oldCpNumber = parentOldItem->data(1,Qt::UserRole).toInt();
    }
    if (newCpNumber == -1) {
      newCpNumber = parentNewItem->data(1,Qt::UserRole).toInt();
    }
    if ((oldCpNumber != newCpNumber) ||
        // Check name
        (parentOldItem->text(0) != parentNewItem->text(0)) ) { 
      // try to optimize
      return false;  
    } else if ((parentOldItem->text(0) != parentNewItem->text(0)) || // Logical Name
               (parentOldItem->text(3) != parentNewItem->text(3))) {   // Check logical name
      return false;
    } else {
      parentOldItem = parentOldItem->parent();
      parentNewItem = parentNewItem->parent();
    }
  } // end while

  return true;
}


void G4OpenGLQtViewer::addNonPVSceneTreeElement(
 const G4String& model
 ,int currentPOIndex
 ,const std::string& modelDescription
 ,const G4Visible& visible
) {

  QString modelShortName = getModelShortName(model);
  G4Colour color;

  // Special case for text
  try {
    const G4Text& g4Text = dynamic_cast<const G4Text&>(visible);
    color = fSceneHandler.GetTextColour(g4Text);
  }
  catch (std::bad_cast) {
    color = fSceneHandler.GetColour(visible);
  }
  if (modelShortName == "") {
    return ;
  }
  // try to init it
  if (fSceneTreeComponentTreeWidget == NULL) {
    initSceneTreeComponent();
  }

  // if no UI
  if (fSceneTreeComponentTreeWidget == NULL) {
    return;
  }

  fSceneTreeComponentTreeWidget->blockSignals(true);

  // Create the "Model" node if not

  QList<QTreeWidgetItem *> resItem;
  resItem =  fSceneTreeComponentTreeWidget->findItems (modelShortName, Qt::MatchExactly, 0 );
  QTreeWidgetItem * currentItem = NULL;
  const PVPath tmpFullPath;

  if (resItem.empty()) {
    currentItem = createTreeWidgetItem(tmpFullPath,
                                       modelShortName,
                                       0, // currentPVCopyNb
                                       -1, // currentPVPOIndex
                                       "",
                                       Qt::Checked,
                                       NULL,
                                       color);
  } else {
    currentItem = resItem.first();
  }

  // Is this volume already in the tree AND PO is not the same?
  const QList<QTreeWidgetItem *>&
  resItems =  fSceneTreeComponentTreeWidget->findItems (QString(modelDescription.c_str()), Qt::MatchFixedString| Qt::MatchCaseSensitive|Qt::MatchRecursive, 0 );
  
  bool alreadyPresent = false;
  for (int i = 0; i < resItems.size(); ++i) {
    if (currentPOIndex == resItems.at(i)->data(0,Qt::UserRole).toInt()) {
      alreadyPresent = true;
    }
  }
  if (!alreadyPresent) {
    createTreeWidgetItem(tmpFullPath,
                         QString(modelDescription.c_str()),
                         0, // currentPVCopyNb
                         currentPOIndex,
                         "",
                         Qt::Checked,
                         currentItem,
                         color);
  }
  fSceneTreeComponentTreeWidget->blockSignals(false);
  
}


/**
   Get the short name for a given label
*/
QString G4OpenGLQtViewer::getModelShortName(const G4String& model) {

  QString modelShortName = model.data();
  if (modelShortName.mid(0,modelShortName.indexOf(" ")) == "G4PhysicalVolumeModel") {
    modelShortName = fTouchableVolumes;
  } else {
    if (modelShortName.mid(0,2) == "G4") {
      modelShortName = modelShortName.mid(2);
    }
    if (modelShortName.indexOf("Model") != -1) {
      modelShortName = modelShortName.mid(0,modelShortName.indexOf("Model"));
    }
  }
  return modelShortName;
}



bool G4OpenGLQtViewer::isTouchableVisible(int POindex){

  // If no scene tree (Immediate viewer)
  if (fSceneTreeComponentTreeWidget == NULL) {
    return false;
  }
  
  // should be the next one
  if (fLastSceneTreeWidgetAskFor != fPositivePoIndexSceneTreeWidgetQuickMap.end()) {
    fLastSceneTreeWidgetAskFor++;
  }

  QTreeWidgetItem* item = getTreeWidgetItem(POindex);
  
  if (item != NULL) {
    if ( item->checkState(0) == Qt::Checked) {
      return true;
    }
  }  
  return false;
}


bool G4OpenGLQtViewer::parseAndCheckVisibility(QTreeWidgetItem * treeNode,int POindex){
  bool isFound = false;
  for (int i = 0; i < treeNode->childCount() ; ++i) {

    if (treeNode->child(i)->data(0,Qt::UserRole).toInt() == POindex) {
      if (treeNode->child(i)->checkState(0) == Qt::Checked) {
        return true;
      }
    }
    isFound = parseAndCheckVisibility(treeNode->child(i),POindex);
    if (isFound) {
      return true;
    }
  } // end for
  return false;
}


void G4OpenGLQtViewer::sceneTreeComponentItemChanged(QTreeWidgetItem* item, int) {

  if (fCheckSceneTreeComponentSignalLock == false) {
    fCheckSceneTreeComponentSignalLock = true;
    G4bool checked = false;
    if (item->checkState(0) == Qt::Checked) {
      checked = true;
    }
    setCheckComponent(item,checked);
    updateQWidget();

    // change vis attributes to set new visibility
    G4int iPO = item->data(0,Qt::UserRole).toInt();
    if (iPO >= 0 && fTreeItemModels.find(iPO) != fTreeItemModels.end()) {
      const PVPath& fullPath = fTreeItemModels[iPO];
      // If a physical volume
      if (fullPath.size()) {

        // Instantiate a working copy of a G4VisAttributes object...
        G4VisAttributes workingVisAtts;
        // and set the visibility.
        workingVisAtts.SetVisibility(checked);

        // Add a vis atts modifier to the view parameters...
        fVP.AddVisAttributesModifier
        (G4ModelingParameters::VisAttributesModifier
         (workingVisAtts,
          G4ModelingParameters::VASVisibility,
          fullPath));
        // G4ModelingParameters::VASVisibility tells G4PhysicalVolumeModel that
        // it is the visibility that should be picked out and merged with the
        // touchable's normal vis attributes.
      }
    }

    fCheckSceneTreeComponentSignalLock = false;
  }
}


void G4OpenGLQtViewer::sceneTreeComponentSelected() {
  /*
  if (fSceneTreeComponentTreeWidget) {
    if (! fSceneTreeComponentTreeWidget->selectedItems().isEmpty ()) {
      QTreeWidgetItem* item = fSceneTreeComponentTreeWidget->selectedItems ().first();
      if (item) {

        const PVPath& fullPath = fTreeItemModels[item->data(0,Qt::UserRole).toInt()];
        if (fullPath.size() > 0) {
	  const G4VisAttributes* visAttr2 = GetApplicableVisAttributes(fullPath.at(fullPath.size()-1).GetPhysicalVolume()->GetLogicalVolume()->GetVisAttributes());
	}
      }
    }
  }
  */
}


void G4OpenGLQtViewer::changeDepthInSceneTree (int val){

  // If no scene tree (Immediate viewer)
  if (fSceneTreeComponentTreeWidget == NULL) {
    return;
  }

  // max depth :   fSceneTreeDepth
  // val is between 0 and 1
  // 0  .1  .2  .3  .4  .5  .6  .7  .8  .9  1
  // 1      1.4          2                  
  // 1         2         3         4

  // Get the depth :
  double depth = 1 + ((double)val)/1000 * ((double)fSceneTreeDepth+1);
  
  // lock update on scene tree items
  fCheckSceneTreeComponentSignalLock = true;
  
  // Disable redraw each time !
  G4bool currentAutoRefresh = fVP.IsAutoRefresh();
  fVP.SetAutoRefresh(false);

  for (int b=0;b<fSceneTreeComponentTreeWidget->topLevelItemCount();b++) {
    changeDepthOnSceneTreeItem(depth,1.,fSceneTreeComponentTreeWidget->topLevelItem(b));
  }

  // Enable redraw !
  fVP.SetAutoRefresh(currentAutoRefresh);
  updateQWidget();

  // unlock update on scene tree items
  fCheckSceneTreeComponentSignalLock = false;

}


void G4OpenGLQtViewer::changeColorAndTransparency(QTreeWidgetItem* item,int) {

  if (item == NULL) {
    return;
  }
  const QColor& old = QColor(item->data(2,Qt::UserRole).value<QColor>());

#if QT_VERSION < 0x040500
  bool a;
  const QColor& color = QColor(QColorDialog::getRgba (old.rgba(),&a,fSceneTreeComponentTreeWidget));
#else
  const QColor& color = QColorDialog::getColor(old,
                                        fSceneTreeComponentTreeWidget,
                                        " Get color and transparency",
                                        QColorDialog::ShowAlphaChannel);
#endif
  
  if (color.isValid()) {

    // change vis attributes to set new colour
    G4int iPO = item->data(0,Qt::UserRole).toInt();
    if (iPO >= 0 && fTreeItemModels.find(iPO) != fTreeItemModels.end()) {
      const PVPath& fullPath = fTreeItemModels[iPO];
      // If a physical volume
      if (fullPath.size()) {
      
        // Instantiate a working copy of a G4VisAttributes object...
        G4VisAttributes workingVisAtts;
        // and set the colour.
        G4Colour g4c(((G4double)color.red())/255,
                     ((G4double)color.green())/255,
                     ((G4double)color.blue())/255,
                     ((G4double)color.alpha())/255);
        workingVisAtts.SetColour(g4c);
        
        // Add a vis atts modifier to the view parameters...
        fVP.AddVisAttributesModifier
        (G4ModelingParameters::VisAttributesModifier
         (workingVisAtts,
          G4ModelingParameters::VASColour,
          fullPath));
        // G4ModelingParameters::VASColour tells G4PhysicalVolumeModel that it is
        // the colour that should be picked out and merged with the touchable's
        // normal vis attributes.
      }
    }
    
    // set scene tree parameters
    changeQColorForTreeWidgetItem(item,color);
  }
}


G4Colour G4OpenGLQtViewer::getColorForPoIndex(int poIndex) {

  QTreeWidgetItem* item = getTreeWidgetItem(poIndex);

  if (item != NULL) {
    
    const QColor& color = item->data(2,Qt::UserRole).value<QColor>();
    G4Colour g4c(((G4double)color.red())/255,
                 ((G4double)color.green())/255,
                 ((G4double)color.blue())/255,
                 ((G4double)color.alpha())/255);

    return g4c;
  }
  return G4Colour();
}


void G4OpenGLQtViewer::changeSearchSelection()
{
  const QString& searchText = fHelpLine->text();
  if (fSceneTreeComponentTreeWidget == NULL) {
    return;
  }

  // unselect all
  for (int a=0; a<fSceneTreeComponentTreeWidget->topLevelItemCount(); a++) {
    fSceneTreeComponentTreeWidget->topLevelItem(a)->setExpanded(false);
    fSceneTreeComponentTreeWidget->topLevelItem(a)->setSelected(false);
    clearSceneTreeSelection(fSceneTreeComponentTreeWidget->topLevelItem(a));
  }

  QList<QTreeWidgetItem *> itemList = fSceneTreeComponentTreeWidget->findItems (searchText,Qt::MatchContains | Qt::MatchRecursive,0);

  for (int i = 0; i < itemList.size(); ++i) {
    QTreeWidgetItem* expandParentItem = itemList.at(i);
    while (expandParentItem->parent() != NULL) {
      expandParentItem->parent()->setExpanded(true);
      expandParentItem = expandParentItem->parent();
    }
    itemList.at(i)->setSelected(true);
  }    

}


void G4OpenGLQtViewer::clearSceneTreeSelection(QTreeWidgetItem* item) {
  for (int a=0; a<item->childCount(); a++) {
    item->child(a)->setSelected(false);
    item->child(a)->setExpanded(false);
    clearSceneTreeSelection(item->child(a));
  }

}


bool G4OpenGLQtViewer::isPVVolume(QTreeWidgetItem* item) {
  QTreeWidgetItem* sParent = item;
  while (sParent->parent() != NULL) {
    sParent = sParent->parent();
  }
  if (sParent->text(0) != fTouchableVolumes) {
    return false;
  }
  // item is the "Touchable" node
  if (item->text(0) == fTouchableVolumes) {
    return false;
  }
  return true;
}


void G4OpenGLQtViewer::changeDepthOnSceneTreeItem(
 double lookForDepth
 ,double currentDepth
 ,QTreeWidgetItem* item
) {
  double transparencyLevel = 0.;
  
  // look for a 2.2 depth and we are at level 3 
  // -> Set all theses items to Opaque
  // ONLY if it is a PV volume !
  if (isPVVolume(item)) {
    if ((lookForDepth-currentDepth) < 0) {
      item->setCheckState(0,Qt::Checked);
      updatePositivePoIndexSceneTreeWidgetQuickMap(item->data(0,Qt::UserRole).toInt(),item);
      transparencyLevel = 1;
    } else if ((lookForDepth-currentDepth) > 1 ){
      item->setCheckState(0,Qt::Unchecked);
      updatePositivePoIndexSceneTreeWidgetQuickMap(item->data(0,Qt::UserRole).toInt(),item);
      transparencyLevel = 0;
    } else {
      item->setCheckState(0,Qt::Checked);
      updatePositivePoIndexSceneTreeWidgetQuickMap(item->data(0,Qt::UserRole).toInt(),item);
      transparencyLevel = 1-(lookForDepth-currentDepth);
    }
  }
  
  if (item->data(0,Qt::UserRole).toInt() >= 0) {
    const G4Colour& color = getColorForPoIndex(item->data(0,Qt::UserRole).toInt());
    
    // We are less depper (ex:tree depth:2) than lookForDepth (ex:3.1)
    // -> Have to hide this level ONLY if it was not hidden before
    
    // Not on a top level item case
    // Do not set if it was already set

    // Should add them all the time in case of an older command has change transparency
    // before. Should be checked in changeDepthInSceneTree for duplicated commands
    // Do not change transparency if not visible by humain (and avoid precision value
    // problems..)
    if (((color.GetAlpha()-transparencyLevel) >  0.000001) ||
        ((color.GetAlpha()-transparencyLevel) < -0.000001))  {
      if ((item->text(3) != "")) {
        // FIXME : Should not test this here because of transparent
        // volume that will came after and with a different alpha level
        // Good thing to do is to check and suppress doubles in changeDepthInSceneTree
        // and then check if last (transparents volumes) has to change alpha

        changeQColorForTreeWidgetItem(item,QColor((int)(color.GetRed()*255),
                                                  (int)(color.GetGreen()*255),
                                                  (int)(color.GetBlue()*255),
                                                  (int)(transparencyLevel*255)));
      }
    }
  }

  for (int b=0;b< item->childCount();b++) {
    changeDepthOnSceneTreeItem(lookForDepth,currentDepth+1,item->child(b));
  }
}


void G4OpenGLQtViewer::clearTreeWidget(){
  // be careful about calling this twice

  if (fSceneTreeComponentTreeWidget) {

    if (fOldSceneTreeComponentTreeWidget == NULL) {
      fOldSceneTreeComponentTreeWidget = new QTreeWidget();
    }
    // clear old scene tree
    int tmp = fOldSceneTreeComponentTreeWidget->topLevelItemCount();
    while (tmp > 0) {
      delete fOldSceneTreeComponentTreeWidget->takeTopLevelItem (0);
      tmp = fOldSceneTreeComponentTreeWidget->topLevelItemCount();
    }

    if (fSceneTreeComponentTreeWidget->topLevelItemCount () > 0) {

      fPVRootNodeCreate = false;

      // reset all old
      fOldPositivePoIndexSceneTreeWidgetQuickMap.clear();
      fOldNullPoIndexSceneTreeWidgetQuickVector.clear();
      fOldTreeItemModels.clear();

      // Clone everything  
      for (int b =0; b <fSceneTreeComponentTreeWidget->topLevelItemCount();b++) {
        // All tree widgets are in :
        // then we could get the old POindex and get 
        // .visible/Hidden
        // .Check/Uncheck
        // .selected
        // .colour status from std::map

        // clone top level items
        int poIndex = fSceneTreeComponentTreeWidget->topLevelItem(b)->data(0,Qt::UserRole).toInt();
        if (poIndex != -1) {
          fOldPositivePoIndexSceneTreeWidgetQuickMap.insert(std::pair <int, QTreeWidgetItem*> (poIndex,cloneWidgetItem(fSceneTreeComponentTreeWidget->topLevelItem(b))));
        } else {
          fOldNullPoIndexSceneTreeWidgetQuickVector.push_back(cloneWidgetItem(fSceneTreeComponentTreeWidget->topLevelItem(b)));
        }

        // clone leaves
        cloneSceneTree(fSceneTreeComponentTreeWidget->topLevelItem(b));
      }
      // delete all elements

      fOldTreeItemModels.insert(fTreeItemModels.begin(), fTreeItemModels.end());  

      // all is copy, then clear scene tree
      int tmp2 = fSceneTreeComponentTreeWidget->topLevelItemCount();
      while (tmp2 > 0) {
        delete fSceneTreeComponentTreeWidget->takeTopLevelItem (0);
        tmp2 = fSceneTreeComponentTreeWidget->topLevelItemCount();
      }
      fPositivePoIndexSceneTreeWidgetQuickMap.clear();

      // put correct value in paramaters
      fOldLastSceneTreeWidgetAskFor = fOldPositivePoIndexSceneTreeWidgetQuickMap.begin();
      fSceneTreeDepth = 1;
      fModelShortNameItem = NULL;
      fMaxPOindexInserted = -1;

    }
  }
}


/**
   Clone : 
   - Open/close
   - Visible/hidden
   - Selected
*/
QTreeWidgetItem * G4OpenGLQtViewer::cloneWidgetItem(QTreeWidgetItem* item) {

  QTreeWidgetItem* cloneItem = new QTreeWidgetItem();

  // Clone what is create createTreeWidgetItem step

  cloneItem->setText(0,item->text(0));
  cloneItem->setData(1,Qt::UserRole,item->data(1,Qt::UserRole).toInt());
  cloneItem->setText(2,item->text(2));
  cloneItem->setData(0, Qt::UserRole,item->data(0,Qt::UserRole).toInt());
  cloneItem->setText(3,item->text(3));
  cloneItem->setFlags(item->flags());
  cloneItem->setToolTip(0,item->toolTip(0));        
  cloneItem->setCheckState(0,item->checkState(0));
  cloneItem->setSelected(item->isSelected()); 
  cloneItem->setExpanded(item->isExpanded ());

  cloneItem->setData(2,Qt::UserRole,item->data(2,Qt::UserRole).value<QColor>());

  return cloneItem;
}


/**
   Clone the current tree in order to get a snapshot of old version
*/
void G4OpenGLQtViewer::cloneSceneTree(
 QTreeWidgetItem* rootItem
) {
  
  for (int b=0;b< rootItem->childCount();b++) {

    QTreeWidgetItem *child = rootItem->child(b);

    // clone top level items
    int poIndex = child->data(0,Qt::UserRole).toInt();
    if (poIndex != -1) {
      fOldPositivePoIndexSceneTreeWidgetQuickMap.insert(std::pair <int, QTreeWidgetItem*> (poIndex,cloneWidgetItem(child)));
    } else {
      fOldNullPoIndexSceneTreeWidgetQuickVector.push_back(cloneWidgetItem(child));
    }
    cloneSceneTree(child);
  }
}


/**
   Update the quick scene tree visibility map (used by parseAndCheckVisibility)
*/
 void G4OpenGLQtViewer::updatePositivePoIndexSceneTreeWidgetQuickMap(int POindex,QTreeWidgetItem* item) {

   // Check state
   std::map <int, QTreeWidgetItem*>::iterator i;
   i = fPositivePoIndexSceneTreeWidgetQuickMap.find(POindex);

   if (i == fPositivePoIndexSceneTreeWidgetQuickMap.end()) {
     fPositivePoIndexSceneTreeWidgetQuickMap.insert(std::pair <int, QTreeWidgetItem*> (POindex,item) );
     fLastSceneTreeWidgetAskFor = fPositivePoIndexSceneTreeWidgetQuickMap.end();
   } else {
     i->second = item;
   }
 }



void G4OpenGLQtViewer::changeQColorForTreeWidgetItem(QTreeWidgetItem* item,const QColor& qc) {

  int POIndex = item->data(0,Qt::UserRole).toInt();
  updatePositivePoIndexSceneTreeWidgetQuickMap(POIndex,item );

  QPixmap pixmap = QPixmap(QSize(16, 16));
  if (item->data(0,Qt::UserRole).toInt() != -1) {
    pixmap.fill (qc);
  } else {
    pixmap.fill (QColor(255,255,255,255));
  }
  QPainter painter(&pixmap);
  painter.setPen(Qt::black);
  painter.drawRect(0,0,15,15); // Draw contour

  item->setIcon(0,pixmap);
  item->setData(2,Qt::UserRole,qc);
}



/**
   @return the corresponding item if existing.
   Look into fPositivePoIndexSceneTreeWidgetQuickMap
 */
QTreeWidgetItem* G4OpenGLQtViewer::getTreeWidgetItem(int POindex){

  // -1 is not a visible item
  if (POindex == -1) {
    return NULL;
  }

  if (fPositivePoIndexSceneTreeWidgetQuickMap.size() == 0){
    return NULL;
  }

  if (POindex == fLastSceneTreeWidgetAskFor->first) {
    if (fLastSceneTreeWidgetAskFor->second != NULL) {
      return fLastSceneTreeWidgetAskFor->second;
    }
  }
  
  // if not, use the "find" algorithm
  fLastSceneTreeWidgetAskFor = fPositivePoIndexSceneTreeWidgetQuickMap.find(POindex);
  
  if (fLastSceneTreeWidgetAskFor != fPositivePoIndexSceneTreeWidgetQuickMap.end()) {
    return fLastSceneTreeWidgetAskFor->second;
  }
  return NULL;
}

/**
   @return the corresponding item if existing in the old tree
   Look into fOldPositivePoIndexSceneTreeWidgetQuickMap
 */
QTreeWidgetItem* G4OpenGLQtViewer::getOldTreeWidgetItem(int POindex){


  // -1 is not a visible item
  if (POindex == -1) {
    return NULL;
  }

  if (fOldPositivePoIndexSceneTreeWidgetQuickMap.size() == 0){
    return NULL;
  }

  // Should be call only once by item addition
  fOldLastSceneTreeWidgetAskFor++;

  if (POindex == fOldLastSceneTreeWidgetAskFor->first) {
    if (fOldLastSceneTreeWidgetAskFor->second != NULL) {
      if (fOldLastSceneTreeWidgetAskFor != fOldPositivePoIndexSceneTreeWidgetQuickMap.end()) {
        return fOldLastSceneTreeWidgetAskFor->second;
      }
    }
  }
  
  // if not, use the "find" algorithm
  fOldLastSceneTreeWidgetAskFor = fOldPositivePoIndexSceneTreeWidgetQuickMap.find(POindex);
  
  if (fOldLastSceneTreeWidgetAskFor != fOldPositivePoIndexSceneTreeWidgetQuickMap.end()) {
    return fOldLastSceneTreeWidgetAskFor->second;
  }
  return NULL;
}



/**
   Should replace actual tree by the one in this class
   and update tree
*/
void G4OpenGLQtViewer::displaySceneTreeComponent() {
  // no UI
  if (fUISceneTreeComponentsTBWidget == NULL) {
    return; 
  }
  if (fSceneTreeComponentTreeWidget == NULL) {
    return; 
  }

  // sort tree items
  fSceneTreeComponentTreeWidget->sortItems (0, Qt::AscendingOrder );

  return;
}


/**
   Update the toolbar Icons/Mouse context menu
   - Change ortho/perspective
   - Change surface style
   - Change cursor style
 */
void G4OpenGLQtViewer::updateToolbarAndMouseContextMenu(){
   if (fBatchMode) {
     return;
   }

  G4ViewParameters::DrawingStyle
  d_style = fVP.GetDrawingStyle();
  
  // Surface style
  if (d_style == G4ViewParameters::wireframe) {
    if (fUiQt) fUiQt->SetIconWireframeSelected();
    if (fContextMenu) {
      fDrawingWireframe->setChecked(true);
      fDrawingLineRemoval->setChecked(false);
      fDrawingSurfaceRemoval->setChecked(false);
      fDrawingLineSurfaceRemoval->setChecked(false);
    }    
  } else if (d_style == G4ViewParameters::hlr) {
    if (fUiQt) fUiQt->SetIconHLRSelected();
    if (fContextMenu) {
      fDrawingLineRemoval->setChecked(true);
      fDrawingWireframe->setChecked(false);
      fDrawingSurfaceRemoval->setChecked(false);
      fDrawingLineSurfaceRemoval->setChecked(false);
    }    
  } else if (d_style == G4ViewParameters::hsr) {
    if (fUiQt) fUiQt->SetIconSolidSelected();
    if (fContextMenu) {
      fDrawingSurfaceRemoval->setChecked(true);
      fDrawingWireframe->setChecked(false);
      fDrawingLineRemoval->setChecked(false);
      fDrawingLineSurfaceRemoval->setChecked(false);
    }
  } else if (d_style == G4ViewParameters::hlhsr) {
    if (fUiQt) fUiQt->SetIconHLHSRSelected();
    if (fContextMenu) {
      fDrawingLineSurfaceRemoval->setChecked(true);
      fDrawingWireframe->setChecked(false);
      fDrawingLineRemoval->setChecked(false);
      fDrawingSurfaceRemoval->setChecked(false);
      fDrawingLineSurfaceRemoval->setChecked(false);
    }
  }
  
  
  // projection style
  G4double d_proj = fVP.GetFieldHalfAngle () ;
  if (d_proj == 0.) { // ortho
    if (fUiQt) fUiQt->SetIconOrthoSelected();
    if (fContextMenu) {
      fProjectionOrtho->setChecked(true);
      fProjectionPerspective->setChecked(false);
    }    
  } else {
    if (fUiQt) fUiQt->SetIconPerspectiveSelected();
      if (fContextMenu) {
        fProjectionPerspective->setChecked(true);
        fProjectionOrtho->setChecked(false);
      }
  }
  
  
  // mouse style : They are controlled by UI !
  if (fUiQt && fContextMenu) {
    if (fUiQt->IsIconPickSelected()) {
      fMousePickAction->setChecked(true);
      fMouseZoomOutAction->setChecked(false);
      fMouseZoomInAction->setChecked(false);
      fMouseRotateAction->setChecked(false);
      fMouseMoveAction->setChecked(false);
    } else if (fUiQt->IsIconZoomOutSelected()) {
      fMouseZoomOutAction->setChecked(true);
      fMousePickAction->setChecked(false);
      fMouseZoomInAction->setChecked(false);
      fMouseRotateAction->setChecked(false);
      fMouseMoveAction->setChecked(false);
    } else if (fUiQt->IsIconZoomInSelected()) {
      fMouseZoomInAction->setChecked(true);
      fMousePickAction->setChecked(false);
      fMouseZoomOutAction->setChecked(false);
      fMouseRotateAction->setChecked(false);
      fMouseMoveAction->setChecked(false);
    } else if (fUiQt->IsIconRotateSelected()) {
      fMouseRotateAction->setChecked(true);
      fMousePickAction->setChecked(false);
      fMouseZoomOutAction->setChecked(false);
      fMouseZoomInAction->setChecked(false);
      fMouseMoveAction->setChecked(false);
    } else if (fUiQt->IsIconMoveSelected()) {
      fMouseMoveAction->setChecked(true);
      fMousePickAction->setChecked(false);
      fMouseZoomOutAction->setChecked(false);
      fMouseZoomInAction->setChecked(false);
      fMouseRotateAction->setChecked(false);
    }
  }
}


/*
  
void MultiLayer::exportToSVG(const QString& fname)
{
QPicture picture;
QPainter p(&picture);
for (int i=0;i<(int)graphsList->count();i++)
{
Graph *gr=(Graph *)graphsList->at(i);
Plot *myPlot= (Plot *)gr->plotWidget();
      
QPoint pos=gr->pos();
      
int width=int(myPlot->frameGeometry().width());
int height=int(myPlot->frameGeometry().height());
      
myPlot->print(&p, QRect(pos,QSize(width,height)));
}
  
p.end();
picture.save(fname, "svg");
}
*/
#endif
