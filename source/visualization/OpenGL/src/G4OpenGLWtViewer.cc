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
// $Id: G4OpenGLWtViewer.cc 86360 2014-11-10 08:34:16Z gcosmo $
//
// 
// G4OpenGLWtViewer : Class to provide Wt specific
//                     functionality for OpenGL in GEANT4
//
// 27/06/2003 : G.Barrand : implementation (at last !).

#ifdef G4VIS_BUILD_OPENGLWT_DRIVER

#include "G4OpenGLWtViewer.hh"
#include "G4VViewer.hh"
#include "G4VSceneHandler.hh"
#include "G4OpenGLSceneHandler.hh"

#include "G4ios.hh"
#include "G4VisExtent.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Point3D.hh"
#include "G4Normal3D.hh"
#include "G4Scene.hh"
//#include "G4OpenGLWtExportDialog.hh"
//#include "G4OpenGLWtMovieDialog.hh"
#include "G4UnitsTable.hh"
#include "G4Wt.hh"
#include "G4UIWt.hh"
#include "G4UImanager.hh"
#include "G4UIcommandTree.hh"
#include <Wt/WHBoxLayout>
#include <Wt/WApplication>
#include <Wt/WTime>




//////////////////////////////////////////////////////////////////////////////
void G4OpenGLWtViewer::CreateMainWindow (
  Wt::WGLWidget* glWidget
 ,Wt::WString name
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLWtViewer::CreateMainWindow \n");
#endif

  if(fWindow) return; //Done.

  fWindow = glWidget ;
  //  fWindow->makeCurrent();

//  G4Wt* interactorManager = G4Wt::getInstance ();
  // return false if G4UIWt was not launch
  
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if (UI == NULL) return;
  
  if (! static_cast<G4UIWt*> (UI->GetG4UIWindow())) {
    // NO UI, should be batch mode
    fBatchMode = true;
    return;
  }
  fUiWt = static_cast<G4UIWt*> (UI->GetG4UIWindow());
  
  bool isTabbedView = false;
  if ( fUiWt) {
    if (!fBatchMode) {
//      if (!interactorManager->IsExternalApp()) {

      // resize window to get the good size at the beginning
      ResizeWindow(fVP.GetWindowSizeHintX(),fVP.GetWindowSizeHintY());
      
      isTabbedView = fUiWt->AddTabWidget(fWindow->parent(),name,getWinWidth(),getWinHeight());
      // change color
      fWindow->parent()->decorationStyle().setBackgroundColor (Wt::WColor("blue"));

      // Have to resize !
#ifdef G4DEBUG_VIS_OGL
      printf("G4OpenGLWtViewer::CreateMainWindow :: resize :%d %d\n",getWinWidth(),getWinHeight());
#endif
      fWindow->resize(getWinWidth(),getWinHeight());

      fUISceneTreeComponentsTBWidget = fUiWt->GetSceneTreeComponentsTBWidget();
      fWindow->resize(fWindow->parent()->width(),fWindow->parent()->height());
      isTabbedView = true;
      //      }
    }
  }
#ifdef G4DEBUG_VIS_OGL
  else {
    printf("G4OpenGLWtViewer::CreateMainWindow :: UIWt NOt found \n");
  }
#endif
  

/* if (!isTabbedView) { // we have to do a dialog
    
    Wt::WWidget *myParent = getParentWidget();
#ifdef G4DEBUG_VIS_OGL
    printf("G4OpenGLWtViewer::CreateMainWindow :: getParent OK \n");
#endif
    if (myParent != NULL) {
      glWidget->setParent(myParent);
    }
    Wt::WHBoxLayout *mainLayout = new Wt::WHBoxLayout(fGLWindow);
    
    mainLayout->setMargin(0);
    mainLayout->setSpacing(0);
    mainLayout->addWidget(fWindow);
    if (fGLWindow->inherits("Wt::WContainerWidget")) {
      fGLWindow->setWindowTitle( name);
    }
    fGLWindow->setLayout(mainLayout);
    
 */   
/*
 //useful for MACOSX, we have to compt the menuBar height
    int offset = QApplication::desktop()->height()
    - QApplication::desktop()->availableGeometry().height();
    
    G4int YPos= fVP.GetWindowAbsoluteLocationHintY(QApplication::desktop()->height());
    if (fVP.GetWindowAbsoluteLocationHintY(QApplication::desktop()->height())< offset) {
      YPos = offset;
    }
#ifdef G4DEBUG_VIS_OGL
    printf("G4OpenGLQtViewer::CreateMainWindow :: resizing to %d %d \n",getWinWidth(), getWinHeight());
#endif
    fGLWindow->move(fVP.GetWindowAbsoluteLocationHintX(QApplication::desktop()->width()),YPos);
*/
//    fGLWindow->show();
//  } else {

    fGLWindow = fWindow;
//  }
  
  if(!fWindow) return;
  
#ifdef _A_FINIR_FIXME
  if (!fContextMenu)
    createPopupMenu();
#endif

}

/**  Close the dialog and set the pointer to NULL
 */
// void G4OpenGLWtViewer::dialogClosed() {
//   //  fGLWindow = NULL;
// }

//////////////////////////////////////////////////////////////////////////////
G4OpenGLWtViewer::G4OpenGLWtViewer (
                                    G4OpenGLSceneHandler& scene
                                    )
  :G4VViewer (scene, -1)
  ,G4OpenGLViewer (scene)
  ,fWindow(0)
  ,fRecordFrameNumber(0)
 //#ifdef _A_FINIR_FIXME  ,fContextMenu(0)
  ,fMouseAction(STYLE1)
  ,fDeltaRotation(1)
  ,fDeltaSceneTranslation(0.01)
  ,fDeltaDepth(0.01)
  ,fDeltaZoom(0.05)
  ,fDeltaMove(0.05)
  ,fHoldKeyEvent(false)
  ,fHoldMoveEvent(false)
  ,fHoldRotateEvent(false)
  ,fAutoMove(false)
  ,fEncoderPath("")
  ,fTempFolderPath("")
  ,fMovieTempFolderPath("")
  ,fSaveFileName("")
  ,fParameterFileName("mpeg_encode_parameter_file.par")
  ,fMovieParametersDialog(NULL)
  ,fRecordingStep(WAIT)
  ,fProcess(NULL)
  ,fNbMaxFramesPerSec(100)
  ,fNbMaxAnglePerSec(360)
  ,fLaunchSpinDelay(100)
  ,fXRot(0)
  ,fYRot(0)
  ,fNoKeyPress(true)
  ,fAltKeyPress(false)
  ,fControlKeyPress(false)
  ,fShiftKeyPress(false)
  ,fBatchMode(false)
  ,fUiWt(NULL)
{

  // launch Wt if not
  G4Wt::getInstance ();

// FIXME : all stuff with G4VIS_BUILD_OPENGL_ES_DRIVER
    //  G4OpenGLViewer::SetView(this);

#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLWtViewer::Create \n");
#endif
  fLastPos3 = Wt::WPoint(-1,-1);    
  fLastPos2 = Wt::WPoint(-1,-1);    
  fLastPos1 = Wt::WPoint(-1,-1);    
  
  mMatrix.setToIdentity();

#ifdef _A_FINIR_FIXME
  initMovieParameters();
#endif

  fLastEventTime = new Wt::WTime();

#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLWtViewer::G4OpenGLWtViewer END\n");
#endif

}





void G4OpenGLWtViewer::resizeGL(int width, int height)
  {
#ifdef G4DEBUG_VIS_OGL
    printf("G4OpenGLWtViewer resizeGL %d %d\n",width,height);
#endif

  }




//////////////////////////////////////////////////////////////////////////////
G4OpenGLWtViewer::~G4OpenGLWtViewer (
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
#ifdef _A_FINIR_FIXME
  G4cout <<removeTempFolder().toUTF8().c_str() <<G4endl;
#endif
}





#ifdef _A_FINIR_FIXME
/**
   Create a popup menu for the widget. This menu is activated by right-mouse click
*/
void G4OpenGLWtViewer::createPopupMenu()    {

  fContextMenu = new WMenu("All");

  WMenu *mMouseAction = fContextMenu->addMenu("&Mouse actions");

  fRotateAction = mMouseAction->addAction("Rotate");
  fMoveAction = mMouseAction->addAction("Move");
  fPickAction = mMouseAction->addAction("Pick");
  WAction *shortcutsAction = mMouseAction->addAction("Show shortcuts");

  fRotateAction->setCheckable(true);
  fMoveAction->setCheckable(false);
  fPickAction->setCheckable(false);
  shortcutsAction->setCheckable(false);

  fRotateAction->setChecked(true);
  fMoveAction->setChecked(false);
  fPickAction->setChecked(false);
  shortcutsAction->setChecked(false);

  WObject ::connect(fRotateAction, 
                    SIGNAL(triggered(bool)),
                    this, 
                    SLOT(actionMouseRotate()));

  WObject ::connect(fMoveAction, 
                    SIGNAL(triggered(bool)),
                    this, 
                    SLOT(actionMouseMove()));

  WObject ::connect(fPickAction, 
                    SIGNAL(triggered(bool)),
                    this, 
                    SLOT(actionMousePick()));

  WObject ::connect(shortcutsAction, 
                    SIGNAL(triggered(bool)),
                    this, 
                    SLOT(showShortcuts()));

  // === Style Menu ===
  WMenu *mStyle = fContextMenu->addMenu("&Style");

  WMenu *mRepresentation = mStyle->addMenu("&Representation");
  WMenu *mProjection = mStyle->addMenu("&Projection");
  WAction *polyhedron = mRepresentation->addAction("Polyhedron");
  WAction *nurbs = mRepresentation->addAction("NURBS");

  WAction *ortho = mProjection->addAction("Orthographic");
  WAction *perspective = mProjection->addAction("Persepective");

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
    createRadioAction(ortho, perspective,SLOT(toggleProjection(bool)),1);
  } else {
    createRadioAction(ortho, perspective,SLOT(toggleProjection(bool)),2);
  }

  // === Drawing Menu ===
  WMenu *mDrawing = mStyle->addMenu("&Drawing");

  fDrawingWireframe = mDrawing->addAction("Wireframe");
  fDrawingWireframe->setCheckable(true);

  fDrawingLineRemoval = mDrawing->addAction("Hidden line removal");
  fDrawingLineRemoval->setCheckable(true);

  fDrawingSurfaceRemoval = mDrawing->addAction("Hidden Surface removal");
  fDrawingSurfaceRemoval->setCheckable(true);

  fDrawingLineSurfaceRemoval = mDrawing->addAction("Hidden line and surface removal");
  fDrawingLineSurfaceRemoval->setCheckable(true);

  // INIT Drawing
  G4ViewParameters::DrawingStyle d_style;
  d_style = fVP.GetDrawingStyle();
  
  if (d_style == G4ViewParameters::wireframe) {
    fDrawingWireframe->setChecked(true);
  } else if (d_style == G4ViewParameters::hlr) {
    fDrawingLineRemoval->setChecked(true);
  } else if (d_style == G4ViewParameters::hsr) {
    fDrawingSurfaceRemoval->setChecked(true);
  } else if (d_style == G4ViewParameters::hlhsr) {
    fDrawingLineSurfaceRemoval->setChecked(true);
  } else {
    mDrawing->clear();
  }
  WObject ::connect(fDrawingWireframe, 
                    SIGNAL(triggered(bool)),
                    this, 
                    SLOT(actionDrawingWireframe()));
  WObject ::connect(fDrawingLineRemoval, 
                    SIGNAL(triggered(bool)),
                    this, 
                    SLOT(actionDrawingLineRemoval()));
  WObject ::connect(fDrawingSurfaceRemoval, 
                    SIGNAL(triggered(bool)),
                    this, 
                    SLOT(actionDrawingSurfaceRemoval()));
  WObject ::connect(fDrawingLineSurfaceRemoval, 
                    SIGNAL(triggered(bool)),
                    this, 
                    SLOT(actionDrawingLineSurfaceRemoval()));

  // Background Color

  WAction *backgroundColorChooser ;

  // === Action Menu ===
  backgroundColorChooser = mStyle->addAction("Background color");
  WObject ::connect(backgroundColorChooser, 
                    SIGNAL(triggered()),
                    this,
                    SLOT(actionChangeBackgroundColor()));

  // Text Color

  WAction *textColorChooser ;
  // === Action Menu ===
  textColorChooser = mStyle->addAction("Text color");
  WObject ::connect(textColorChooser, 
                    SIGNAL(triggered()),
                    this,
                    SLOT(actionChangeTextColor()));

  // Default Color

  WAction *defaultColorChooser ;
  // === Action Menu ===
  defaultColorChooser = mStyle->addAction("Default color");
  WObject ::connect(defaultColorChooser, 
                    SIGNAL(triggered()),
                    this,
                    SLOT(actionChangeDefaultColor()));


  // === Action Menu ===
  WMenu *mActions = fContextMenu->addMenu("&Actions");
  WAction *createEPS = mActions->addAction("Save as ...");
  WObject ::connect(createEPS, 
                    SIGNAL(triggered()),
                    this,
                    SLOT(actionSaveImage()));

  // === Action Menu ===
  WAction *movieParameters = mActions->addAction("Movie parameters...");
  WObject ::connect(movieParameters, 
                    SIGNAL(triggered()),
                    this,
                    SLOT(actionMovieParameters()));




  // === Special Menu ===
  WMenu *mSpecial = fContextMenu->addMenu("S&pecial");
  WMenu *mTransparency = mSpecial->addMenu("Transparency");
  WAction *transparencyOn = mTransparency->addAction("On");
  WAction *transparencyOff = mTransparency->addAction("Off");

  if (transparency_enabled == false) {
    createRadioAction(transparencyOn,transparencyOff,SLOT(toggleTransparency(bool)),2);
  } else if (transparency_enabled == true) {
    createRadioAction(transparencyOn,transparencyOff,SLOT(toggleTransparency(bool)),1);
  } else {
    mSpecial->clear();
  }


  WMenu *mAntialiasing = mSpecial->addMenu("Antialiasing");
  WAction *antialiasingOn = mAntialiasing->addAction("On");
  WAction *antialiasingOff = mAntialiasing->addAction("Off");

  if (antialiasing_enabled == false) {
    createRadioAction(antialiasingOn,antialiasingOff,SLOT(toggleAntialiasing(bool)),2);
  } else if (antialiasing_enabled == true) {
    createRadioAction(antialiasingOn,antialiasingOff,SLOT(toggleAntialiasing(bool)),1);
  } else {
    mAntialiasing->clear();
  }

  WMenu *mHaloing = mSpecial->addMenu("Haloing");
  WAction *haloingOn = mHaloing->addAction("On");
  WAction *haloingOff = mHaloing->addAction("Off");

  if (haloing_enabled == false) {
    createRadioAction(haloingOn,haloingOff,SLOT(toggleHaloing(bool)),2);
  } else if (haloing_enabled == true) {
    createRadioAction(haloingOn,haloingOff,SLOT(toggleHaloing(bool)),1);
  } else {
    mHaloing->clear();
  }

  WMenu *mAux = mSpecial->addMenu("Auxiliary edges");
  WAction *auxOn = mAux->addAction("On");
  WAction *auxOff = mAux->addAction("Off");
  if (!fVP.IsAuxEdgeVisible()) {
    createRadioAction(auxOn,auxOff,SLOT(toggleAux(bool)),1);
  } else {
    createRadioAction(auxOn,auxOff,SLOT(toggleAux(bool)),2);
  }



  WMenu *mFullScreen = mSpecial->addMenu("&Full screen");
  fFullScreenOn = mFullScreen->addAction("On");
  fFullScreenOff = mFullScreen->addAction("Off");
  createRadioAction(fFullScreenOn,fFullScreenOff,SLOT(toggleFullScreen(bool)),2);

}


void G4OpenGLWtViewer::G4manageContextMenuEvent(WContextMenuEvent *e)
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
void G4OpenGLWtViewer::createRadioAction(WAction *action1,WAction *action2, const std::string& method,unsigned int nCheck) {

  action1->setCheckable(true);
  action2->setCheckable(true);

  if (nCheck ==1)
    action1->setChecked (true);
  else
    action2->setChecked (true);
   
  WObject ::connect(action1, SIGNAL(triggered(bool)),action2, SLOT(toggle()));
  WObject ::connect(action2, SIGNAL(triggered(bool)),action1, SLOT(toggle()));

  WObject ::connect(action1, SIGNAL(toggled(bool)),this, method.c_str());

}

/**
   Slot activate when mouseAction->rotate menu is set 
 */
void G4OpenGLWtViewer::actionMouseRotate() {
  emit toggleMouseAction(STYLE1);
}


/**
   Slot activate when mouseAction->rotate menu is set 
 */
void G4OpenGLWtViewer::actionMouseMove() {
  emit toggleMouseAction(STYLE2);
}


/**
   Slot activate when mouseAction->zoom menu is set 
 */
void G4OpenGLWtViewer::actionMousePick() {
  emit toggleMouseAction(STYLE3);
}


/**
   Slot activate when drawing->wireframe menu is set 
 */
void G4OpenGLWtViewer::actionDrawingWireframe() {
  emit toggleDrawingAction(1);
}

/**
   Slot activate when drawing->line removal menu is set 
 */
void G4OpenGLWtViewer::actionDrawingLineRemoval() {
  emit toggleDrawingAction(2);
}

/**
   Slot activate when drawing->surface removal menu is set 
 */
void G4OpenGLWtViewer::actionDrawingSurfaceRemoval() {
  emit toggleDrawingAction(3);
}

/**
   Slot activate when drawing->wireframe menu is set 
 */
void G4OpenGLWtViewer::actionDrawingLineSurfaceRemoval() {
  emit toggleDrawingAction(4);
}


/**
   Slot activated when mouse action is toggle
   @param aAction : STYLE1, STYLE2, STYLE3
 */
void G4OpenGLWtViewer::toggleMouseAction(mouseActions aAction) {
  
  if ((aAction == STYLE1) || //initialize all
      (aAction == STYLE2) ||
      (aAction == STYLE3))  {
    fRotateAction->setChecked (false);
    fMoveAction->setChecked (false);
    fPickAction->setChecked (false);
    fVP.SetPicking(false);
    fMouseAction = aAction;
  }
  // rotate
  if (aAction == STYLE1) {  // rotate
    showShortcuts();
    fRotateAction->setChecked (true);
  } else  if (aAction == STYLE2) { //move
    fMoveAction->setChecked (true);
  } else  if (aAction == STYLE3) { //pick
    fPickAction->setChecked (true);
    fVP.SetPicking(true);
  }
}

#endif
/**
   Show shortcuts for this mouse action
 */
void G4OpenGLWtViewer::showShortcuts() {
  G4cout << "========= Mouse Shortcuts =========" << G4endl;
  if (fMouseAction == STYLE1) {  // rotate
    G4cout << "Click and move mouse to rotate volume " << G4endl;
    G4cout << "ALT + Click and move mouse to rotate volume (View Direction)" << G4endl;
    G4cout << "CTRL + Click and zoom mouse to zoom in/out" << G4endl;
    G4cout << "SHIFT + Click and zoommove camera point of view" << G4endl;
  } else  if (fMouseAction == STYLE2) { //move
    G4cout << "Move camera point of view with mouse" << G4endl;
  } else  if (fMouseAction == STYLE3) { //pick
    G4cout << "Click and pick " << G4endl;
  }
  G4cout << "========= Move Shortcuts =========" << G4endl;
  G4cout << "Press left/right arrows to move volume left/right" << G4endl;
  G4cout << "Press up/down arrows to move volume up/down" << G4endl;
  G4cout << "Press '+'/'-' to move volume toward/forward" << G4endl;
  G4cout <<  G4endl;
  G4cout << "========= Rotation (Theta/Phi) Shortcuts =========" << G4endl;
  G4cout << "Press SHIFT + left/right arrows to rotate volume left/right" << G4endl;
  G4cout << "Press SHIFT + up/down arrows to rotate volume up/down" << G4endl;
  G4cout <<  G4endl;
  G4cout << "========= Rotation (View Direction) Shortcuts =========" << G4endl;
  G4cout << "Press ALT + left/right to rotate volume around vertical direction" << G4endl;
  G4cout << "Press ALT + up/down to rotate volume around horizontal direction" << G4endl;
  G4cout <<  G4endl;
  G4cout << "========= Zoom View =========" << G4endl;
  G4cout << "Press CTRL + '+'/'-' to zoom into volume" << G4endl;
  G4cout <<  G4endl;
  G4cout << "========= Misc =========" << G4endl;
  G4cout << "Press ALT +/- to slow/speed rotation/move" << G4endl;
  G4cout << "Press H to reset view" << G4endl;
  G4cout << "Press Esc to exit FullScreen" << G4endl;
  G4cout <<  G4endl;
  G4cout << "========= Video =========" << G4endl;
  G4cout << "In video mode : " << G4endl;
  G4cout << " Press SPACE to Start/Pause video recording " << G4endl;
  G4cout << " Press RETURN to Stop video recording " << G4endl;
  G4cout <<  G4endl;
}



#ifdef _A_FINIR_FIXME

/**
   Slot activated when drawing menu is toggle
   Warning : When G4OpenGLStoredWtViewer::DrawView() method call,
   KernelVisitDecision () will be call and will set the fNeedKernelVisit
   to 1. See G4XXXStoredViewer::CompareForKernelVisit for explanations.
   It will cause a redraw of the view
   @param aAction : 1 wireframe, 2 line removal, 3 surface removal, 4 line & surface removal
   @see G4OpenGLStoredWtViewer::DrawView
   @see G4XXXStoredViewer::CompareForKernelVisit
 */
void G4OpenGLWtViewer::toggleDrawingAction(int aAction) {

  G4ViewParameters::DrawingStyle d_style = G4ViewParameters::wireframe;
  

  // initialize
  if ((aAction >0) && (aAction <5)) {
    fDrawingWireframe->setChecked (false);
    fDrawingLineRemoval->setChecked (false);
    fDrawingSurfaceRemoval->setChecked (false);
    fDrawingLineSurfaceRemoval->setChecked (false);
  }
  if (aAction ==1) {
    fDrawingWireframe->setChecked (true);

    d_style = G4ViewParameters::wireframe;

  } else  if (aAction ==2) {
    fDrawingLineRemoval->setChecked (true);

    d_style = G4ViewParameters::hlr;

  } else  if (aAction ==3) {
    fDrawingSurfaceRemoval->setChecked (true);

    d_style = G4ViewParameters::hsr;

  } else  if (aAction ==4) {
    fDrawingLineSurfaceRemoval->setChecked (true);
    d_style = G4ViewParameters::hlhsr;
  }
  fVP.SetDrawingStyle(d_style);

  updateWWidget();
}


/**
   SLOT Activate by a click on the representation menu
   Warning : When G4OpenGLStoredWtViewer::DrawView() method call,
   KernelVisitDecision () will be call and will set the fNeedKernelVisit
   to 1. See G4XXXStoredViewer::CompareForKernelVisit for explanations.
   It will cause a redraw of the view
   @param check : 1 polyhedron, 0 nurbs
   @see G4OpenGLStoredWtViewer::DrawView
   @see G4XXXStoredViewer::CompareForKernelVisit
*/
void G4OpenGLWtViewer::toggleRepresentation(bool check) {

  G4ViewParameters::RepStyle style;
  if (check == 1) {
    style = G4ViewParameters::polyhedron;
  } else {
    style = G4ViewParameters::nurbs;
  }
  fVP.SetRepStyle (style);

  updateWWidget();
}

/**
   SLOT Activate by a click on the projection menu
   Warning : When G4OpenGLStoredWtViewer::DrawView() method call,
   KernelVisitDecision () will be call and will set the fNeedKernelVisit
   to 1. See G4XXXStoredViewer::CompareForKernelVisit for explanations.
   It will cause a redraw of the view
   @param check : 1 orthographic, 2 perspective
   @see G4OpenGLStoredWtViewer::DrawView
   @see G4XXXStoredViewer::CompareForKernelVisit
*/
void G4OpenGLWtViewer::toggleProjection(bool check) {

  if (check == 1) {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/projection o");
  } else {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/set/projection p");
  }  
  updateWWidget();
}


/**
   SLOT Activate by a click on the transparency menu
@param check : 1 , 0
*/
void G4OpenGLWtViewer::toggleTransparency(bool check) {
  
  if (check) {
    transparency_enabled = false;
  } else {
    transparency_enabled = true;
  }
  SetNeedKernelVisit (true);
  updateWWidget();
}

/**
   SLOT Activate by a click on the antialiasing menu
@param check : 1 , 0
*/
void G4OpenGLWtViewer::toggleAntialiasing(bool check) {

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

  updateWWidget();
}

/**
   SLOT Activate by a click on the haloing menu
@param check : 1 , 0
*/
//FIXME : I SEE NOTHING...
void G4OpenGLWtViewer::toggleHaloing(bool check) {
  if (check) {
    haloing_enabled = false;
  } else {
    haloing_enabled = true;
  }

  updateWWidget();

}

/**
   SLOT Activate by a click on the auxiliaire edges menu
@param check : 1 , 0
*/
void G4OpenGLWtViewer::toggleAux(bool check) {
  if (check) {
    fVP.SetAuxEdgeVisible(true);
  } else {
    fVP.SetAuxEdgeVisible(false);
  }
  SetNeedKernelVisit (true);
  updateWWidget();
}

/**
   SLOT Activate by a click on the full screen menu
*/
void G4OpenGLWtViewer::toggleFullScreen(bool check) {
  if (check != fGLWindow->isFullScreen()) { //toggle
    fGLWindow->setWindowState(fGLWindow->windowState() ^ Wt::WindowFullScreen);
    G4cerr << "This version of Wt could not do fullScreen. Resizing the widget is the only solution available." << G4endl;
  }
}
#endif

#ifdef _A_FINIR_FIXME
void G4OpenGLWtViewer::savePPMToTemp() {
  if (fMovieTempFolderPath == "") {
    return;
  }
  Wt::WString fileName ="Test"+Wt::WString::number(fRecordFrameNumber)+".ppm";
  Wt::WString filePath =fMovieTempFolderPath+fileName;

  WImage image;
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



void G4OpenGLWtViewer::actionSaveImage() {
  Wt::WString filters;
  WList<WByteArray> formats =  WImageWriter::supportedImageFormats ();
  for (int i = 0; i < formats.size(); ++i) {
    filters +=formats.at(i) + ";;";
  }
  filters += "eps;;";
  filters += "ps;;";
  filters += "pdf";
  Wt::WString* selectedFormat = new Wt::WString();
  std::string name;
  name =  WFileDialog::getSaveFileName ( fGLWindow,
                                                    tr("Save as ..."),
                                                    ".",
                                                    filters,
                                                    selectedFormat ).toUTF8().c_str(); 
  // bmp jpg jpeg png ppm xbm xpm
  if (name.empty()) {
    return;
  }
  name += "." + selectedFormat->toUTF8();
  Wt::WString format = selectedFormat->toLower();
  setPrintFilename(name.c_str(),0);
  G4OpenGLWtExportDialog* exportDialog= new G4OpenGLWtExportDialog(fGLWindow,format,fWindow->height(),fWindow->width());
  if(  exportDialog->exec()) {

    WImage image;
    bool res = false;
    if ((exportDialog->getWidth() !=fWindow->width()) ||
        (exportDialog->getHeight() !=fWindow->height())) {
      setPrintSize(exportDialog->getWidth(),exportDialog->getHeight());
      if ((format != Wt::WString("eps")) && (format != Wt::WString("ps"))) {
      G4cerr << "Export->Change Size : This function is not implemented, to export in another size, please resize your frame to what you need" << G4endl;
      
      //    rescaleImage(exportDialog->getWidth(),exportDialog->getHeight());// re-scale image
      //      WGLWidget* glResized = fWindow;

      // FIXME :
      // L.Garnier : I've try to implement change size function, but the problem is 
      // the renderPixmap function call the WGLWidget to resize and it doesn't draw
      // the content of this widget... It only draw the background.

      //      fWindow->renderPixmap (exportDialog->getWidth()*2,exportDialog->getHeight()*2,true );

      //      WPixmap pixmap = fWindow->renderPixmap ();
      
      //      image = pixmap->toImage();
      //      glResized->resize(exportDialog->getWidth()*2,exportDialog->getHeight()*2);
      //      image = glResized->grabFrameBuffer();
      }      
    } else {
      image = fWindow->grabFrameBuffer();
    }    
    if (format == Wt::WString("eps")) {
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
      res = image.save(Wt::WString(name.c_str()),0,exportDialog->getSliderValue());
    } else {
      G4cerr << "This version of G4UI Could not generate the selected format" << G4endl;
    }
    if ((format == Wt::WString("eps")) && (format == Wt::WString("ps"))) {
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
#endif


#ifdef _A_FINIR_FIXME
void G4OpenGLWtViewer::actionChangeBackgroundColor() {

  //   //I need to revisit the kernel if the background colour changes and
  //   //hidden line removal is enabled, because hlr drawing utilises the
  //   //background colour in its drawing...
  //   // (Note added by JA 13/9/2005) Background now handled in view
  //   // parameters.  A kernel visit is triggered on change of background.

  WColor color;
  color = WColorDialog::getColor(Wt::black, fGLWindow);
  if (color.isValid()) {
    Wt::WString com = "/vis/viewer/set/background ";
    Wt::WString num;
    com += num.setNum(((float)color.red())/256)+" ";
    com += num.setNum(((float)color.green())/256)+" ";
    com += num.setNum(((float)color.blue())/256)+" ";
    G4UImanager::GetUIpointer()->ApplyCommand(com.toUTF8().c_str());
    updateWWidget();
  }
}

void G4OpenGLWtViewer::actionChangeTextColor() {

  WColor color;
  color = WColorDialog::getColor(Wt::yellow, fGLWindow);
  if (color.isValid()) {
    Wt::WString com = "/vis/viewer/set/defaultTextColour ";
    Wt::WString num;
    com += num.setNum(((float)color.red())/256)+" ";
    com += num.setNum(((float)color.green())/256)+" ";
    com += num.setNum(((float)color.blue())/256)+" ";
    G4UImanager::GetUIpointer()->ApplyCommand(com.toUTF8().c_str());
    updateWWidget();
  }
}

void G4OpenGLWtViewer::actionChangeDefaultColor() {

  WColor color;
  color = WColorDialog::getColor(Wt::white, fGLWindow);
  printf("actionChangeDefaultColor\n");
  if (color.isValid()) {
    Wt::WString com = "/vis/viewer/set/defaultColour ";
    Wt::WString num;
    com += num.setNum(((float)color.red())/256)+" ";
    com += num.setNum(((float)color.green())/256)+" ";
    com += num.setNum(((float)color.blue())/256)+" ";
    G4UImanager::GetUIpointer()->ApplyCommand(com.toUTF8().c_str());
    updateWWidget();
  }
}


void G4OpenGLWtViewer::actionMovieParameters() {
  showMovieParametersDialog();
}


void G4OpenGLWtViewer::showMovieParametersDialog() {
  if (!fMovieParametersDialog) {
    fMovieParametersDialog= new G4OpenGLWtMovieDialog(this,fGLWindow);
    displayRecordingStatus();
    fMovieParametersDialog->checkEncoderSwParameters();
    fMovieParametersDialog->checkSaveFileNameParameters();
    fMovieParametersDialog->checkTempFolderParameters();
    if (getEncoderPath() == "") {
      setRecordingInfos("mpeg_encode is needed to encode in video format. It is available here: http://bmrc.berkeley.edu/frame/research/mpeg/");
    }
  }
  fMovieParametersDialog->show();
}
#endif

/*
// http://www.google.com/codesearch?hl=en&q=+jpg+Wt+quality+WDialog+show:FZkUoth8oiw:TONpW2mR-_c:tyTfrKMO-xI&sa=N&cd=2&ct=rc&cs_p=http://soft.proindependent.com/src/qtiplot-0.8.9.zip&cs_f=qtiplot-0.8.9/qtiplot/src/application.cpp#a0

void Graph::exportToSVG(const Wt::WString& fname)
{
  // enable workaround for Wt3 misalignments
  WwtPainter::setSVGMode(true);
  WPicture picture;
  WPainter p(&picture);
  d_plot->print(&p, d_plot->rect());
  p.end();

  picture.save(fname, "svg");
}
*/



/**
   Save the current mouse press point
   @param p mouse click point
*/
void G4OpenGLWtViewer::G4MousePressEvent(Wt::WMouseEvent *event)
{
    if (event->button() & Wt::WMouseEvent::LeftButton) {
#ifdef _A_FINIR_FIXME
    fWindow->setMouseTracking(true);
#endif
    fAutoMove = false; // stop automove
    fLastPos1 = Wt::WPoint(event->widget().x,event->widget().y);
    fLastPos2 = fLastPos1;
    fLastPos3 = fLastPos2;
//    fLastEventTime->start();
    if (fMouseAction == STYLE3){  // pick
      Pick(event->widget().x,event->widget().y);
    }
  }
}

/**
*/
void G4OpenGLWtViewer::G4MouseReleaseEvent()
{
    fSpinningDelay = 1;//fLastEventTime->elapsed();
  Wt::WPoint delta = Wt::WPoint(fLastPos3.x()-fLastPos1.x(),fLastPos3.y()-fLastPos1.y());
  if ((delta.x() == 0) && (delta.y() == 0)) {
    return;
  }
  if (fSpinningDelay < fLaunchSpinDelay ) {
    fAutoMove = true;
    Wt::WTime lastMoveTime;
//    lastMoveTime.start();
    // try to addapt speed move/rotate looking to drawing speed
    float correctionFactor = 5;
    while (fAutoMove) {
//      if ( lastMoveTime.elapsed () >= (int)(1000/fNbMaxFramesPerSec)) {
        if ( 1 >= (int)(1000/fNbMaxFramesPerSec)) {
        float lTime = 1000/((float)1);
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
                
        // Check Wt Versions for META Keys
                
        // Click and move mouse to rotate volume
        // ALT + Click and move mouse to rotate volume (View Direction)
        // SHIFT + Click and move camera point of view
        // CTRL + Click and zoom mouse to zoom in/out

        if (fMouseAction == STYLE1) {  // rotate
          if (fNoKeyPress) {
            rotateWtScene(((float)delta.x())/correctionFactor,((float)delta.y())/correctionFactor);
          } else if (fAltKeyPress) {
            rotateWtSceneToggle(((float)delta.x())/correctionFactor,((float)delta.y())/correctionFactor);
          }
          
        } else if (fMouseAction == STYLE2) {  // move
          moveScene(-((float)delta.x())/correctionFactor,-((float)delta.y())/correctionFactor,0,true);
        }
//        lastMoveTime.start();
      }
#ifdef _A_FINIR_FIXME
      ((Wt::WApplication*)G4Wt::getInstance ())->processEvents();
#endif
    }
  }
#ifdef _A_FINIR_FIXME
  fWindow->setMouseTracking(false);
#endif
}


void G4OpenGLWtViewer::G4MouseDoubleClickEvent()
{
#ifdef _A_FINIR_FIXME
  fWindow->setMouseTracking(true);
#endif
}


/**
   @param pos_x mouse x position
   @param pos_y mouse y position
   @param mButtons mouse button active
   @param mAutoMove true: apply this move till another evnt came, false :one time move
*/

 void G4OpenGLWtViewer::G4MouseMoveEvent(Wt::WMouseEvent *event)
{
  
  Wt::WMouseEvent::Button mButtons = event->button();

#ifdef _A_FINIR_FIXME
  updateKeyModifierState(event->modifiers());
#endif

  if (fAutoMove) {
    return;
  }

  fLastPos3 = fLastPos2;
  fLastPos2 = fLastPos1;
  fLastPos1 = Wt::WPoint(event->widget().x, event->widget().y);

  printf("G4OpenGLWtViewer move :%d %d\n",event->widget().x, event->widget().y);
  int deltaX = fLastPos2.x()-fLastPos1.x();
  int deltaY = fLastPos2.y()-fLastPos1.y();

  if (fMouseAction == STYLE1) {  // rotate
    if (mButtons & Wt::WMouseEvent::LeftButton) {
      if (fNoKeyPress) {
        rotateWtScene(((float)deltaX),((float)deltaY));
      } else if (fAltKeyPress) {
        rotateWtSceneToggle(((float)deltaX),((float)deltaY));
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
  } else if (fMouseAction == STYLE2) {  // move
    if (mButtons & Wt::WMouseEvent::LeftButton) {
      moveScene(-deltaX,-deltaY,0,true);
    }
  }

//  fLastEventTime->start();
}


/**
   Move the scene of dx, dy, dz values.
   @param dx delta mouse x position
   @param dy delta mouse y position
   @param mouseMove : true if even comes from a mouse move, false if even comes from key action
*/

void G4OpenGLWtViewer::moveScene(float dx,float dy, float dz,bool mouseMove)
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
    coefTrans = getSceneNearWidth()*fDeltaSceneTranslation;
    coefDepth = getSceneDepth()*fDeltaDepth;
  }
  fVP.IncrementPan(-dx*coefTrans,dy*coefTrans,dz*coefDepth);
  updateWWidget();
  if (fAutoMove)
#ifdef _A_FINIR_FIXME
    ((WApplication*)G4Wt::getInstance ())->processEvents();
#endif
  
  fHoldMoveEvent = false;
}


/**
   @param dx delta mouse x position
   @param dy delta mouse y position
*/

void G4OpenGLWtViewer::rotateWtScene(float dx, float dy)
{
  if (fHoldRotateEvent)
    return;
  fHoldRotateEvent = true;
  
  if( dx != 0) {
    rotateScene(dx,0);
  }
  if( dy != 0) {
    rotateScene(0,dy);
  }
  updateWWidget();
  
  fHoldRotateEvent = false;
}

/**
   @param dx delta mouse x position
   @param dy delta mouse y position
*/

void G4OpenGLWtViewer::rotateWtSceneToggle(float dx, float dy)
{
  if (fHoldRotateEvent)
    return;
  fHoldRotateEvent = true;
  
  rotateSceneToggle(dx,dy);
  
  updateWWidget();
  
  fHoldRotateEvent = false;
}

/**
   @param dx delta mouse x position
   @param dy delta mouse y position
*/




/** This is the benning of a rescale function. It does nothing for the moment
    @param aWidth : new width
    @param aHeight : new height
*/
void G4OpenGLWtViewer::rescaleImage(
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



#ifdef _A_FINIR_FIXME
/**
   Generate Postscript or PDF form image
   @param aFilename : name of file
   @param aInColor : numbers of colors : 1->BW 2->RGB
   @param aImage : Image to print
*/
bool G4OpenGLWtViewer::printPDF (
 const std::string aFilename
,int aInColor
,WImage aImage
)
{

  WPrinter printer;
  //  printer.setPageSize(pageSize);

  // FIXME : L. Garnier 4/12/07
  // This is not working, it does nothing. Image is staying in color mode
  // So I have desactivate the B/W button in GUI
  if ((!aImage.isGrayscale ()) &&(aInColor ==1 )) {
    aImage = aImage.convertToFormat ( aImage.format(), Wt::MonoOnly);
  }


  if (aFilename.substr(aFilename.size()-3) == ".ps") {
#if WT_VERSION > 0x040200
    printer.setOutputFormat(WPrinter::PostScriptFormat);
#endif
  } else {
#if WT_VERSION > 0x040100
    printer.setOutputFormat(WPrinter::PdfFormat);
#endif
  }
#if WT_VERSION > 0x040100
  printer.setOutputFileName(Wt::WString(aFilename.c_str()));
#endif
  //  printer.setFullPage ( true);
  WPainter paint(&printer);
  paint.drawImage (0,0,aImage);
  paint.end();
  return true;
}



void G4OpenGLWtViewer::G4wheelEvent (Wt::WWheelEvent * event)
{
  fVP.SetZoomFactor(fVP.GetZoomFactor()+(fVP.GetZoomFactor()*(event->delta())/1200)); 
  updateWWidget();
}
#endif

 void G4OpenGLWtViewer::G4keyPressEvent (Wt::WKeyEvent * event) 
{
  if (fHoldKeyEvent)
    return;

  fHoldKeyEvent = true;

  
  // with no modifiers
#ifdef _A_FINIR_FIXME
  updateKeyModifierState(event->modifiers());
#endif
  if ((fNoKeyPress)) { // || (event->modifiers() == Wt::KeyboardModifier )) {
    if (event->key() == Wt::Key_Down) { // go down
      moveScene(0,1,0,false);
    }
    else if (event->key() == Wt::Key_Up) {  // go up
      moveScene(0,-1,0,false);
    }
    if (event->key() == Wt::Key_Left) { // go left
      moveScene(-1,0,0,false);
    }
    else if (event->key() == Wt::Key_Right) { // go right
      moveScene(1,0,0,false);
    }
    if (event->text() == Wt::WString("-") ) { // go backward
      moveScene(0,0,1,false);
    }
    else if (event->text() == Wt::WString("+")) { // go forward
      moveScene(0,0,-1,false);
    }

    // escaped from full screen
    if (event->key() == Wt::Key_Escape) {
#ifdef _A_FINIR_FIXME
      toggleFullScreen(false);
#endif
    }
  }    
  // several case here : If return is pressed, in every case -> display the movie parameters dialog
  // If one parameter is wrong -> put it in red (only save filenam could be wrong..)
  // If encoder not found-> does nothing.Only display a message in status box
  // If all ok-> generate parameter file
  // If ok -> put encoder button enabled
  
#ifdef _A_FINIR_FIXME
  if ( (event->key() == Wt::Key_Enter)){ // end of video
    stopVideo();
  }
  if (event->key() == Wt::Key_Space){ // start/pause of video
    startPauseVideo();
  }
#endif
  
  // H : Return Home view
  if (event->key() == Wt::Key_H){ // go Home
    fDeltaRotation = 1;
    fDeltaSceneTranslation = 0.01;
    fDeltaDepth = 0.01;
    fDeltaZoom = 0.05;
    fDeltaMove = 0.05;
    
    fVP.SetZoomFactor(1.);
    fVP.SetUpVector(G4Vector3D (0., 1., 0.));
    fVP.SetViewAndLights (G4Vector3D (0., 0., 1.));

    updateWWidget();
  }

  // Shift Modifier
  if (fShiftKeyPress) {
    if (event->key() == Wt::Key_Down) { // rotate phi
      rotateWtScene(0,-fDeltaRotation);
    }
    else if (event->key() == Wt::Key_Up) { // rotate phi
      rotateWtScene(0,fDeltaRotation);
    }
    if (event->key() == Wt::Key_Left) { // rotate theta
      rotateWtScene(fDeltaRotation,0);
    }
    else if (event->key() == Wt::Key_Right) { // rotate theta
      rotateWtScene(-fDeltaRotation,0);
    }

  // Alt Modifier
  }
  if ((fAltKeyPress)) {
    if (event->key() == Wt::Key_Down) { // rotate phi
      rotateWtSceneToggle(0,-fDeltaRotation);
    }
    else if (event->key() == Wt::Key_Up) { // rotate phi
      rotateWtSceneToggle(0,fDeltaRotation);
    }
    if (event->key() == Wt::Key_Left) { // rotate theta
      rotateWtSceneToggle(fDeltaRotation,0);
    }
    else if (event->key() == Wt::Key_Right) { // rotate theta
      rotateWtSceneToggle(-fDeltaRotation,0);
    }

    // Rotatio +/-
    if (event->text() == Wt::WString("+")) {
      fDeltaRotation = fDeltaRotation/0.7;
      G4cout << "Auto-rotation set to : " << fDeltaRotation << G4endl;
    }
    else if (event->text() == Wt::WString("-")) {
      fDeltaRotation = fDeltaRotation*0.7;
      G4cout << "Auto-rotation set to : " << fDeltaRotation << G4endl;
    }

  // Control Modifier OR Command on MAC
  }
  if ((fControlKeyPress)) {
    if (event->text() == Wt::WString("+")) {
      fVP.SetZoomFactor(fVP.GetZoomFactor()*(1+fDeltaZoom)); 
      updateWWidget();
    }
    else if (event->text() == Wt::WString("-")) {
      fVP.SetZoomFactor(fVP.GetZoomFactor()*(1-fDeltaZoom)); 
      updateWWidget();
    }
  }  
  
  fHoldKeyEvent = false;
}
  

#ifdef _A_FINIR_FIXME
void  G4OpenGLWtViewer::updateKeyModifierState(Wt::KeyboardModifiers modifier) {
  // Check Wt Versions for META Keys
    
  fNoKeyPress = true;
  fAltKeyPress = false;
  fShiftKeyPress = false;
  fControlKeyPress = false;
  
  if (modifier & Wt::AltModifier ) {
    fAltKeyPress = true;
    fNoKeyPress = false;
  }
  if (modifier & Wt::ShiftModifier ) {
    fShiftKeyPress = true;
    fNoKeyPress = false;
  }
  if (modifier & Wt::ControlModifier ) {
    fControlKeyPress = true;
    fNoKeyPress = false;
  }
}

/** Stop the video. Check all parameters and enable encoder button if all is ok.
*/
void G4OpenGLWtViewer::stopVideo() {

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
void G4OpenGLWtViewer::saveVideo() {

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
void G4OpenGLWtViewer::startPauseVideo() {
   
  // first time, if temp parameter is wrong, display parameters dialog and return

  if (( fRecordingStep == WAIT)) {
    if ( fRecordFrameNumber == 0) {
      if (getTempFolderPath() == "") { // BAD_OUTPUT
        showMovieParametersDialog();
        setRecordingInfos("You should specified the temp folder in order to make movie");
        return;
      } else  {
        // remove temp folder if it was create
        Wt::WString tmp = removeTempFolder();
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
  if ((fRecordingStep == WAIT)) {
    setRecordingStatus(START); 
  } else if (fRecordingStep == START) {
    setRecordingStatus(PAUSE);
  } else if (fRecordingStep == PAUSE) {
    setRecordingStatus(CONTINUE);
  } else if (fRecordingStep == CONTINUE) {
    setRecordingStatus(PAUSE);
  }
}

void G4OpenGLWtViewer::setRecordingStatus(RECORDING_STEP step) {

  fRecordingStep = step;
  displayRecordingStatus();
}


void G4OpenGLWtViewer::displayRecordingStatus() {
  
  Wt::WString txtStatus = "";
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
    G4cout << txtStatus.toUTF8().c_str() << G4endl;
  }
  setRecordingInfos("");
}


void G4OpenGLWtViewer::setRecordingInfos(Wt::WString txt) {
  if (fMovieParametersDialog) {
    fMovieParametersDialog->setRecordingInfos(txt);
  } else {
    G4cout << txt.toUTF8().c_str() << G4endl;
  }
}

/** Init the movie parameters. Temp dir and encoder path
*/
void G4OpenGLWtViewer::initMovieParameters() {
  //init encoder
  
   //look for encoderPath
     fProcess = new WProcess();
     
     WObject ::connect(fProcess,SIGNAL(finished ( int)),
		       this,SLOT(processLookForFinished()));
     fProcess->setReadChannelMode(WProcess::MergedChannels);
     fProcess->start ("which mpeg_encode");
  
}

/** @return encoder path or "" if it does not exist
 */
Wt::WString G4OpenGLWtViewer::getEncoderPath() {
  return fEncoderPath;
}
 

/**
 * set the new encoder path
 * @return "" if correct. The error otherwise
*/
Wt::WString G4OpenGLWtViewer::setEncoderPath(Wt::WString path) {
  if (path == "") {
    return "File does not exist";
  }

  path =  WDir::cleanPath(path);
  WFileInfo *f = new WFileInfo(path);
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

  if ((fRecordingStep == BAD_ENCODER)) {
    setRecordingStatus(STOP);
  } 
  return "";
}


bool G4OpenGLWtViewer::isRecording(){
  if ((fRecordingStep == START) || (fRecordingStep == CONTINUE)) {
    return true;
  }
  return false;
}

bool G4OpenGLWtViewer::isPaused(){
  if (fRecordingStep == PAUSE) {
    return true;
  }
  return false;
}

bool G4OpenGLWtViewer::isEncoding(){
  if (fRecordingStep == ENCODING) {
    return true;
  }
  return false;
}

bool G4OpenGLWtViewer::isWaiting(){
  if (fRecordingStep == WAIT) {
    return true;
  }
  return false;
}

bool G4OpenGLWtViewer::isStopped(){
  if (fRecordingStep == STOP) {
    return true;
  }
  return false;
}

bool G4OpenGLWtViewer::isFailed(){
  if (fRecordingStep == FAILED) {
    return true;
  }
  return false;
}

bool G4OpenGLWtViewer::isSuccess(){
  if (fRecordingStep == SUCCESS) {
    return true;
  }
  return false;
}

bool G4OpenGLWtViewer::isBadEncoder(){
  if (fRecordingStep == BAD_ENCODER) {
    return true;
  }
  return false;
}
bool G4OpenGLWtViewer::isBadTmp(){
  if (fRecordingStep == BAD_TMP) {
    return true;
  }
  return false;
}
bool G4OpenGLWtViewer::isBadOutput(){
  if (fRecordingStep == BAD_OUTPUT) {
    return true;
  }
  return false;
}

void G4OpenGLWtViewer::setBadEncoder(){
  fRecordingStep = BAD_ENCODER;
  displayRecordingStatus();
}
void G4OpenGLWtViewer::setBadTmp(){
  fRecordingStep = BAD_TMP;
  displayRecordingStatus();
}
void G4OpenGLWtViewer::setBadOutput(){
  fRecordingStep = BAD_OUTPUT;
  displayRecordingStatus();
}

void G4OpenGLWtViewer::setWaiting(){
  fRecordingStep = WAIT;
  displayRecordingStatus();
}


bool G4OpenGLWtViewer::isReadyToEncode(){
  if (fRecordingStep == READY_TO_ENCODE) {
    return true;
  }
  return false;
}

void G4OpenGLWtViewer::resetRecording() {
    setRecordingStatus(WAIT);
}

/**
 * set the temp folder path
 * @return "" if correct. The error otherwise
*/
Wt::WString G4OpenGLWtViewer::setTempFolderPath(Wt::WString path) {

  if (path == "") {
    return "Path does not exist";
  }
  path =  WDir::cleanPath(path);
  WFileInfo *d = new WFileInfo(path);
  if (!d->exists()) {
    return "Path does not exist";
  } else if (!d->isDir()) {
    return "This is not a directory";
  } else if (!d->isReadable()) {
    return path +" is read protected";
  } else if (!d->isWritable()) {
    return path +" is write protected";
  }
  
  if ((fRecordingStep == BAD_TMP)) {
    setRecordingStatus(WAIT); 
  }
  fTempFolderPath = path;
  return "";
}

/** @return the temp folder path or "" if it does not exist
 */
Wt::WString G4OpenGLWtViewer::getTempFolderPath() {
  return fTempFolderPath;
}
 
/**
 * set the save file name path
 * @return "" if correct. The error otherwise
*/
Wt::WString G4OpenGLWtViewer::setSaveFileName(Wt::WString path) {

  if (path == "") {
    return "Path does not exist";
  }
  
  WFileInfo *file = new WFileInfo(path);
  WDir dir = file->dir();
  path =  WDir::cleanPath(path);
  if (file->exists()) {
    return "File already exist, please choose a new one";
  } else if (!dir.exists()) {
    return "Dir does not exist";
  } else if (!dir.isReadable()) {
    return path +" is read protected";
  }
  
  if ((fRecordingStep == BAD_OUTPUT)) {
    setRecordingStatus(STOP); 
  }
  fSaveFileName = path;
  return "";
}

/** @return the save file path
 */
Wt::WString G4OpenGLWtViewer::getSaveFileName() {
  return fSaveFileName ;
}

/** Create a Wt_temp folder in the temp folder given
* The temp folder will be like this /tmp/WtMovie_12-02-2008_12_12_58/
* @return "" if success. Error message if not.
*/
Wt::WString G4OpenGLWtViewer::createTempFolder() {
  fMovieTempFolderPath = "";
  //check
  Wt::WString tmp = setTempFolderPath(fTempFolderPath);
  if (tmp != "") {
    return tmp;
  }
  Wt::WString sep = Wt::WString(WDir::separator());
  Wt::WString path = sep+"WtMovie_"+WDateTime::currentDateTime ().toString("dd-MM-yyyy_hh-mm-ss")+sep; 
  WDir *d = new WDir(WDir::cleanPath(fTempFolderPath));
  // check if it is already present
  if (d->exists(path)) {
    return "Folder "+path+" already exists.Please remove it first";
  }
  if (d->mkdir(fTempFolderPath+path)) {
    fMovieTempFolderPath = fTempFolderPath+path;
    return "";
  } else {
    return "Can't create "+fTempFolderPath+path;
  }
  return "-";
}

/** Remove the Wt_temp folder in the temp folder
*/
Wt::WString G4OpenGLWtViewer::removeTempFolder() {
	// remove files in Wt_temp folder
  if (fMovieTempFolderPath == "") {
    return "";
  }
  WDir *d = new WDir(WDir::cleanPath(fMovieTempFolderPath));
  if (!d->exists()) {
    return "";  // already remove
  }

  d->setFilter( WDir::Files );
  Wt::WStringList subDirList = d->entryList();
  int res = true;
  Wt::WString error = "";
  for (Wt::WStringList::ConstIterator it = subDirList.begin() ;(it != subDirList.end()) ; it++) {
    const Wt::WString currentFile = *it;
      if (!d->remove(currentFile)) {
        res = false;
        Wt::WString file = fMovieTempFolderPath+currentFile;
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



bool G4OpenGLWtViewer::hasPendingEvents () {
#ifdef _A_FINIR_FIXME
  return ((WApplication*)G4Wt::getInstance ())->hasPendingEvents ();
#endif
  return false;
}

bool G4OpenGLWtViewer::generateMpegEncoderParameters () {

		// save the parameter file
  FILE* fp;
  fp = fopen (Wt::WString(fMovieTempFolderPath+fParameterFileName).toUTF8().c_str(), "w");

  if (fp == NULL) {
    setRecordingInfos("Generation of parameter file failed");
    return false;
  }

  fprintf (fp,"# parameter file template with lots of comments to assist you\n");
  fprintf (fp,"#\n");
  fprintf (fp,"# you can use this as a template, copying it to a separate file then modifying\n");
  fprintf (fp,"# the copy\n");
  fprintf (fp,"#\n");
  fprintf (fp,"#\n");
  fprintf (fp,"# any line beginning with '#' is a comment\n");
  fprintf (fp,"#\n");
  fprintf (fp,"# no line should be longer than 255 characters\n");
  fprintf (fp,"#\n");
  fprintf (fp,"#\n");
  fprintf (fp,"# general format of each line is:\n");
  fprintf (fp,"#	  \n");
  fprintf (fp,"#\n");
  fprintf (fp,"# lines can generally be in any order\n");
  fprintf (fp,"#\n");
  fprintf (fp,"# an exception is the option 'INPUT' which must be followed by input\n");
  fprintf (fp,"# files in the order in which they must appear, followed by 'END_INPUT'\n");
  fprintf (fp,"#\n");
  fprintf (fp,"# Also, if you use the `command` method of generating input file names,\n");
  fprintf (fp,"# the command will only be executed in the INPUT_DIR if INPUT_DIR preceeds\n");
  fprintf (fp,"# the INPUT parameter.\n");
  fprintf (fp,"#\n");
  fprintf (fp,"#  MUST be in UPPER CASE\n");
  fprintf (fp,"#\n");
  fprintf (fp,"\n");
  fprintf (fp,"# Pattern affects speed, quality and compression. See the User's Guide\n");
  fprintf (fp,"# for more info.\n");
  fprintf (fp,"\n");
  fprintf (fp,"PATTERN		IBBPBBPBBPBBPBBP\n");
  fprintf (fp,"OUTPUT		%s\n",getSaveFileName().toUTF8().c_str());
  fprintf (fp,"\n");
  fprintf (fp,"# mpeg_encode really only accepts 3 different file formats, but using a\n");
  fprintf (fp,"# conversion statement it can effectively handle ANY file format\n");
  fprintf (fp,"#\n");
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
  fprintf (fp,"INPUT_DIR	%s\n",fMovieTempFolderPath.toUTF8().c_str());
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
  fprintf (fp,"# These are the Wscale values for the entire frame in variable bit-rate\n");
  fprintf (fp,"# mode, and starting points (but not important) for constant bit rate\n");
  fprintf (fp,"#\n");
  fprintf (fp,"\n");
  fprintf (fp,"# Wscale (Wuantization scale) affects quality and compression,\n");
  fprintf (fp,"# but has very little effect on speed.\n");
  fprintf (fp,"\n");
  fprintf (fp,"IWSCALE		4\n");
  fprintf (fp,"PWSCALE		5\n");
  fprintf (fp,"BWSCALE		12\n");
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
  fprintf (fp,"# ASPECT_RATIO, USER_DATA, GAMMA, IWTABLE, etc.\n");
  fprintf (fp,"\n");
  fprintf (fp,"\n");
  fclose (fp);

  setRecordingInfos("Parameter file "+fParameterFileName+" generated in "+fMovieTempFolderPath);
  setRecordingStatus(READY_TO_ENCODE);
  return true;
}

void G4OpenGLWtViewer::encodeVideo()
{
  if ((getEncoderPath() != "") && (getSaveFileName() != "")) {
    setRecordingStatus(ENCODING);
    
    fProcess = new WProcess();
    WObject ::connect(fProcess,SIGNAL(finished ( int)),
                      this,SLOT(processEncodeFinished()));
    WObject ::connect(fProcess,SIGNAL(readyReadStandardOutput ()),
                      this,SLOT(processEncodeStdout()));
    fProcess->setReadChannelMode(WProcess::MergedChannels);
    fProcess->start (fEncoderPath, Wt::WStringList(fMovieTempFolderPath+fParameterFileName));
  }
}


// FIXME : does not work on Wt3
void G4OpenGLWtViewer::processEncodeStdout()
{
  Wt::WString tmp = fProcess->readStdout ().data();
  int start = tmp.findRev("ESTIMATED TIME");
  tmp = tmp.mid(start,tmp.find("\n",start)-start);
  setRecordingInfos(tmp);
}


void G4OpenGLWtViewer::processEncodeFinished()
{

  Wt::WString txt = "";
  txt = getProcessErrorMsg();
  if (txt == "") {
    setRecordingStatus(SUCCESS);
  } else {
    setRecordingStatus(FAILED);
  }
  //  setRecordingInfos(txt+removeTempFolder());
}


void G4OpenGLWtViewer::processLookForFinished() 
 {

  Wt::WString txt = getProcessErrorMsg();
  if (txt != "") {
    fEncoderPath = "";
  } else {
    fEncoderPath = Wt::WString(fProcess->readAllStandardOutput ().data()).trimmed();
    // if not found, return "not found"
    if (fEncoderPath.contains(" ")) {
      fEncoderPath = "";
    } else if (!fEncoderPath.contains("mpeg_encode")) {
      fEncoderPath = "";
    }
    setEncoderPath(fEncoderPath);
  }
  // init temp folder
  setTempFolderPath(WDir::temp ().absolutePath ());
}


Wt::WString G4OpenGLWtViewer::getProcessErrorMsg()
{
  Wt::WString txt = "";
  if (fProcess->exitCode() != 0) {
    switch (fProcess->error()) {
    case WProcess::FailedToStart:
      txt = "The process failed to start. Either the invoked program is missing, or you may have insufficient permissions to invoke the program.\n";
      break;
    case WProcess::Crashed:
      txt = "The process crashed some time after starting successfully.\n";
      break;
    case WProcess::Timedout:
      txt = "The last waitFor...() function timed out. The state of WProcess is unchanged, and you can try calling waitFor...() again.\n";
      break;
    case WProcess::WriteError:
      txt = "An error occurred when attempting to write to the process. For example, the process may not be running, or it may have closed its input channel.\n";
      break;
    case WProcess::ReadError:
      txt = "An error occurred when attempting to read from the process. For example, the process may not be running.\n";
      break;
    case WProcess::UnknownError:
      txt = "An unknown error occurred. This is the default return value of error().\n";
      break;
    }
  }
   return txt;
}
#endif


//
// Matrix utility functions
//

/*
function mvTranslate(v) {
  multMatrix(Matrix.Translation($V([v[0], v[1], v[2]])).ensure4x4());
}


var mvMatrixStack = [];

function mvPushMatrix(m) {
  if (m) {
    mvMatrixStack.push(m.dup());
    mvMatrix = m.dup();
  } else {
    mvMatrixStack.push(mvMatrix.dup());
  }
}

function mvPopMatrix() {
  if (!mvMatrixStack.length) {
    throw("Can't pop from an empty matrix stack.");
  }
  
  mvMatrix = mvMatrixStack.pop();
  return mvMatrix;
}

function mvRotate(angle, v) {
  var inRadians = angle * Math.PI / 180.0;
  
  var m = Matrix.Rotation(inRadians, $V([v[0], v[1], v[2]])).ensure4x4();
  multMatrix(m);
}
*/



/*
  
void MultiLayer::exportToSVG(const Wt::WString& fname)
{
  WPicture picture;
  WPainter p(&picture);
  for (int i=0;i<(int)graphsList->count();i++)
    {
      Graph *gr=(Graph *)graphsList->at(i);
      Plot *myPlot= (Plot *)gr->plotWidget();
      
      Wt::WPoint pos=gr->pos();
      
      int width=int(myPlot->frameGeometry().width());
      int height=int(myPlot->frameGeometry().height());
      
      myPlot->print(&p, WRect(pos,WSize(width,height)));
    }
  
  p.end();
  picture.save(fname, "svg");
}
*/
#endif




 /*
G4UIWt::CommandEnteredCallback
G4UIWt::CommandEnteredCallback 1
G4UIWt::CommandEnteredCallback 2
G4UIWt::CommandEnteredCallback 3
G4UIWt::CommandEnteredCallback 4
G4VisCommandViewerCreate::SetNewValue Before CreateViewer
G4VisManager::CreateViewer Before CreateViewer
G4OpenGLImmediateWt::CreateViewer 
G4OpenGLImmediateWt::CreateViewer after Get Pointer
G4OpenGLImmediateWt::CreateViewer uiWt
G4UIWt::AddTabWidget 50 50
G4UIWt::AddTabWidget 4
G4UIWt::AddTabWidget 5
G4UIWt::AddTabWidget 5a
G4UIWt::AddTabWidget 5a2 69882928 69985360
G4UIWt::AddTabWidget 5b
G4UIWt::AddTabWidget 5c
G4UIWt::AddTabWidget 6
G4UIWt::AddTabWidget ADD 50 50 + 23 1880279432---------------------------------------------------
G4UIWt::AddTabWidget 7
G4UIWt::AddTabWidget 8
G4UIWt::AddTabWidget 9
G4UIWt::AddTabWidget END
G4OpenGLViewer:: Creation
G4OpenGLWtViewer::Create 
G4OpenGLWtViewer::G4OpenGLWtViewer END
G4OpenGLImmediateWtViewer INIT
G4OpenGLImmediateWt::CreateViewer lastInsert :68536784
G4OpenGLImmediateWt::CreateViewer END 
G4VisManager::CreateViewer After 1 CreateViewer
G4VisManager::CreateViewer After 1 CreateViewer
G4VisManager::CreateViewer After 2 CreateViewer
G4VisManager::CreateViewer After 3 CreateViewer
G4VisManager::CreateViewer After 4 CreateViewer
G4OpenGLImmediateWtViewer::Initialise 
G4OpenGLWtViewer::CreateMainWindow 
G4OpenGLViewer::ResizeWindow 600 600
G4OpenGLViewer::SetWinSize 600 600
G4VisManager::CreateViewer After 5 CreateViewer
G4VisManager::CreateViewer After 5 CreateViewer
G4VisManager::CreateViewer After 6 CreateViewer
G4VisManager::CreateViewer After 7 CreateViewer
G4VisManager::CreateViewer After 8 CreateViewer
G4VisManager::CreateViewer After 9 CreateViewer
G4VisManager::CreateViewer After END CreateViewer
G4VisCommandViewerCreate::SetNewValue After CreateViewer 1
G4VisCommandViewerCreate::SetNewValue After CreateViewer 2
G4VisCommandViewerCreate::SetNewValue After CreateViewer 3
G4OpenGLViewer::SetView
G4OpenGLViewer::ResizeGLView 600 600 0x104857b40
G4OpenGLViewer::ClearView
G4OpenGLViewer::ClearView set Background :0.040000 .4 .9: 1.000000
G4OpenGLViewer::ClearView flush
G4OpenGLImmediateWtViewer DrawView
G4OpenGLImmediateWtViewer updateWWidget
G4OpenGLImmediateWtViewer paintGL vvvvvvvvvvvvvvvvvvvv
  G4OpenGLWtViewer paintGL   VVVVVVVVVVVVVVVVVVVVVVVVVV
    G4OpenGLViewer::SetView
    G4OpenGLViewer::ResizeGLView 600 600 0x104857b40
    G4OpenGLViewer::ClearView
    G4OpenGLViewer::ClearView set Background :0.080000 .4 .9: 0.000000
    G4OpenGLViewer::ClearView flush
#0  G4OpenGLViewer::ClearView (this=0x1049f1540) at src/G4OpenGLViewer.cc:195
#1  0x0000000100353e1d in G4OpenGLWtViewer::paintGL (this=0x1049f0e00) at src/G4OpenGLWtViewer.cc:508
#2  0x000000010030e2bd in G4OpenGLImmediateWtViewer::paintGL (this=0x1049f0e00) at src/G4OpenGLImmediateWtViewer.cc:117
#3  0x000000010030e1b9 in G4OpenGLImmediateWtViewer::updateWWidget (this=0x1049f0e00) at src/G4OpenGLImmediateWtViewer.cc:314
#4  0x000000010030e28e in G4OpenGLImmediateWtViewer::DrawView (this=0x1049f0e00) at src/G4OpenGLImmediateWtViewer.cc:206
#5  0x00000001005fed54 in G4VisCommandViewerRefresh::SetNewValue (this=0x104240e30, newValue=@0x1045436f0) at src/G4VisCommandsViewer.cc:1189
#6  0x0000000101612e04 in ~G4String [inlined] () at src/G4UIcommand.cc:211
#7  0x0000000101612e04 in ~G4String [inlined] () at /Users/garnier/Work/geant4/source/global/management/include/G4String.hh:122
#8  0x0000000101612e04 in G4UIcommand::DoIt (this=<value temporarily unavailable, due to optimizations>, parameterList=<value temporarily unavailable, due to optimizations>) at src/G4UIcommand.cc:211

    G4OpenGLWtViewer drawScene   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      G4OpenGLWtViewer drawArrays
      G4OpenGLWtViewer drawArrays
      G4OpenGLWtViewer drawArrays
      G4OpenGLWtViewer drawArrays
      G4OpenGLWtViewer drawScene Call ComputeView
      G4OpenGLWtViewer::ComputeView 600 600   VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
        G4OpenGLWtViewer drawArrays
        G4OpenGLWtViewer::ComputeView NeedKernelVisit
        G4OpenGLWtViewer::ComputeView ProcessView
        G4VViewer::ProcessView  need ? 1
        G4VSceneHandler::ProcessScene 
        G4OpenGLWtViewer::FinishView() 
        G4VSceneHandler::ProcessScene END
        G4VViewer::ProcessView END
        G4OpenGLWtViewer::FinishView() 
      G4OpenGLWtViewer::ComputeView 600 600 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
      G4OpenGLWtViewer drawScene END Call ComputeView
    G4OpenGLWtViewer drawScene END   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    G4OpenGLImmediateWtViewer resizeGL
    G4OpenGLWtViewer resizeGL 600 600
  G4OpenGLWtViewer paintGL   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
G4OpenGLImmediateWtViewer paintGL ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
G4OpenGLImmediateWtViewer updateWWidget END
G4OpenGLImmediateWtViewer DrawView END
G4VisCommandViewerCreate::SetNewValue After CreateViewer 4
G4UIWt::CommandEnteredCallback 5
G4UIWt::CommandEnteredCallback 6
G4UIWt::CommandEnteredCallback END
G4UIWt::CommandLineSlot
G4OpenGLWtViewer initializeGL
G4OpenGLWtViewer centerpoint END
G4OpenGLWtViewer initializeGL END
G4OpenGLImmediateWtViewer updateWWidget
G4OpenGLImmediateWtViewer paintGL vvvvvvvvvvvvvvvvvvvv
  G4OpenGLWtViewer paintGL   VVVVVVVVVVVVVVVVVVVVVVVVVV
    G4OpenGLViewer::SetView
    G4OpenGLViewer::ResizeGLView 600 600 0x104857b40
    G4OpenGLViewer::ClearView
    G4OpenGLViewer::ClearView set Background :0.120000 .4 .9: 0.000000
    G4OpenGLViewer::ClearView flush
#0  G4OpenGLViewer::ClearView (this=0x1049f1540) at src/G4OpenGLViewer.cc:195
#1  0x0000000100353e1d in G4OpenGLWtViewer::paintGL (this=0x1049f0e00) at src/G4OpenGLWtViewer.cc:508
#2  0x000000010030e2bd in G4OpenGLImmediateWtViewer::paintGL (this=0x1049f0e00) at src/G4OpenGLImmediateWtViewer.cc:117
#3  0x000000010030e1b9 in G4OpenGLImmediateWtViewer::updateWWidget (this=0x1049f0e00) at src/G4OpenGLImmediateWtViewer.cc:314
#4  0x0000000100354e82 in G4OpenGLWtViewer::initializeGL (this=0x1049f0e00) at src/G4OpenGLWtViewer.cc:446
#5  0x000000010283b785 in basic_string [inlined] () at /usr/include/c++/4.2.1/bits/basic_string.h:451
#6  0x000000010283b785 in _Alloc_hider [inlined] () at /Users/garnier/Work/Devel/wt-3.2.0/src/Wt/WGLWidget.C:2067
#7  0x000000010283b785 in basic_string [inlined] () at /Users/garnier/Work/Devel/wt-3.2.0/src/Wt/WGLWidget.C:262
#8  0x000000010283b785 in std::basic_stringbuf<char, std::char_traits<char>, std::allocator<char> >::str () at /Users/garnier/Work/Devel/wt-3.2.0/src/Wt/WGLWidget.C:130
#9  std::basic_stringstream<char, std::char_traits<char>, std::allocator<char> >::str () at /usr/include/c++/4.2.1/sstream:572


    G4OpenGLWtViewer drawScene   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      G4OpenGLWtViewer drawArrays
      G4OpenGLWtViewer drawArrays
      G4OpenGLWtViewer drawArrays
      G4OpenGLWtViewer drawArrays
      G4OpenGLWtViewer drawScene Call ComputeView
      G4OpenGLWtViewer::ComputeView 600 600   VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
        G4OpenGLWtViewer drawArrays
        G4OpenGLWtViewer::ComputeView NeedKernelVisit
        G4OpenGLWtViewer::ComputeView ProcessView
        G4VViewer::ProcessView  need ? 1
        G4VSceneHandler::ProcessScene 
        G4OpenGLWtViewer::FinishView() 
        G4VSceneHandler::ProcessScene END
        G4VViewer::ProcessView END
        G4OpenGLWtViewer::FinishView() 
      G4OpenGLWtViewer::ComputeView 600 600 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
      G4OpenGLWtViewer drawScene END Call ComputeView
    G4OpenGLWtViewer drawScene END   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    G4OpenGLImmediateWtViewer resizeGL
    G4OpenGLWtViewer resizeGL 600 600
  G4OpenGLWtViewer paintGL   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
G4OpenGLImmediateWtViewer paintGL ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
G4OpenGLImmediateWtViewer updateWWidget END
G4OpenGLImmediateWtViewer resizeGL
G4OpenGLWtViewer resizeGL 600 600
G4OpenGLImmediateWtViewer paintGL vvvvvvvvvvvvvvvvvvvv
  G4OpenGLWtViewer paintGL   VVVVVVVVVVVVVVVVVVVVVVVVVV
    G4OpenGLViewer::SetView
    G4OpenGLViewer::ResizeGLView 600 600 0x104857b40
    G4OpenGLViewer::ClearView
    G4OpenGLViewer::ClearView set Background :0.160000 .4 .9: 0.000000
    G4OpenGLViewer::ClearView flush
#0  G4OpenGLViewer::ClearView (this=0x1049f1540) at src/G4OpenGLViewer.cc:195
#1  0x0000000100353e1d in G4OpenGLWtViewer::paintGL (this=0x1049f0e00) at src/G4OpenGLWtViewer.cc:508
#2  0x000000010030e2bd in G4OpenGLImmediateWtViewer::paintGL (this=0x1049f0e00) at src/G4OpenGLImmediateWtViewer.cc:117
#3  0x000000010283955b in std::basic_stringstream<char, std::char_traits<char>, std::allocator<char> >::str () at /usr/include/c++/4.2.1/sstream:527
#4  0x000000010283955b in Wt::WGLWidget::updateDom (this=<value temporarily unavailable, due to optimizations>, element=@0x1049e1800, all=true) at /Users/garnier/Work/Devel/wt-3.2.0/src/Wt/WGLWidget.C:531
#5  0x000000010283ba3c in Wt::WGLWidget::createDomElement (this=0x1049f0e00, app=<value temporarily unavailable, due to optimizations>) at /Users/garnier/Work/Devel/wt-3.2.0/src/Wt/WGLWidget.C:473
#6  0x0000000102993cfa in Wt::WWidget::createSDomElement (this=0x1049f0e00, app=0x105021e00) at /Users/garnier/Work/Devel/wt-3.2.0/src/Wt/WWidget.C:326
#7  0x00000001027ebf16 in Wt::WContainerWidget::createDomChildren (this=0x1041c1b50, parent=@0x1049e0c00, app=0x105021e00) at /Users/garnier/Work/Devel/wt-3.2.0/src/Wt/WContainerWidget.C:720
#8  0x00000001027ec5a6 in Wt::WContainerWidget::createDomElement (this=0x1041c1b50, app=0x105021e00, addChildren=true) at /Users/garnier/Work/Devel/wt-3.2.0/src/Wt/WContainerWidget.C:653
#9  0x0000000102993cfa in Wt::WWidget::createSDomElement (this=0x1041c1b50, app=0x105021e00) at /Users/garnier/Work/Devel/wt-3.2.0/src/Wt/WWidget.C:326
#10 0x00000001027ebf16 in Wt::WContainerWidget::createDomChildren (this=0x1041c1cf0, parent=@0x1049f5200, app=0x105021e00) at /Users/garnier/Work/Devel/wt-3.2.0/src/Wt/WContainerWidget.C:720
#11 0x00000001027ec5a6 in Wt::WContainerWidget::createDomElement (this=0x1041c1cf0, app=0x105021e00, addChildren=true) at /Users/garnier/Work/Devel/wt-3.2.0/src/Wt/WContainerWidget.C:653
#12 0x0000000102993cfa in Wt::WWidget::createSDomElement (this=0x1041c1cf0, app=0x105021e00) at /Users/garnier/Work/Devel/wt-3.2.0/src/Wt/WWidget.C:326
#13 0x00000001027eb5a2 in Wt::WContainerWidget::updateDomChildren (this=0x104278560, parent=@0x1049f4600, app=0x105021e00) at /Users/garnier/Work/Devel/wt-3.2.0/src/Wt/WContainerWidget.C:748
#14 0x00000001027ec62d in Wt::WContainerWidget::getDomChanges (this=0x104278560, result=@0x1045449b0, app=0x105021e00) at /Users/garnier/Work/Devel/wt-3.2.0/src/Wt/WContainerWidget.C:632
#15 0x0000000102976c93 in Wt::WWebWidget::getSDomChanges (this=0x104278560, result=@0x1045449b0, app=0x105021e00) at /Users/garnier/Work/Devel/wt-3.2.0/src/Wt/WWebWidget.C:1788


    G4OpenGLWtViewer drawScene   vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      G4OpenGLWtViewer drawArrays
      G4OpenGLWtViewer drawArrays
      G4OpenGLWtViewer drawArrays
      G4OpenGLWtViewer drawArrays
      G4OpenGLWtViewer drawScene Call ComputeView
      G4OpenGLWtViewer::ComputeView 600 600   VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
        G4OpenGLWtViewer drawArrays
        G4OpenGLWtViewer::ComputeView NeedKernelVisit
        G4OpenGLWtViewer::ComputeView ProcessView
        G4VViewer::ProcessView  need ? 1
        G4VSceneHandler::ProcessScene 
        G4OpenGLWtViewer::FinishView() 
        G4VSceneHandler::ProcessScene END
        G4VViewer::ProcessView END
        G4OpenGLWtViewer::FinishView()  
      G4OpenGLWtViewer::ComputeView 600 600 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
      G4OpenGLWtViewer drawScene END Call ComputeView
    G4OpenGLWtViewer drawScene END   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    G4OpenGLImmediateWtViewer resizeGL
    G4OpenGLWtViewer resizeGL 600 600
  G4OpenGLWtViewer paintGL   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
G4OpenGLImmediateWtViewer paintGL ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 */
