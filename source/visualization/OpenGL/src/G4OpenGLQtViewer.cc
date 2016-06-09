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
// $Id: G4OpenGLQtViewer.cc,v 1.7 2007/11/15 18:24:28 lgarnier Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// 
// G4OpenGLQtViewer : Class to provide Qt specific
//                     functionality for OpenGL in GEANT4
//
// 27/06/2003 : G.Barrand : implementation (at last !).

#ifdef G4VIS_BUILD_OPENGLQT_DRIVER

#include "G4OpenGLQtViewer.hh"

#include "G4ios.hh"
#include "G4VisExtent.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Point3D.hh"
#include "G4Normal3D.hh"
#include "G4Scene.hh"
#include "G4OpenGLQtExportDialog.hh"

#include "G4Qt.hh"
#include "G4UIsession.hh"
#include "G4UImanager.hh"
#include <qapplication.h>
#include <qlayout.h>
#include <qdialog.h>

#if QT_VERSION >= 0x040000
#include <qmenu.h>
#include <qimagewriter.h>
#else
#include <qaction.h>
#include <qwidgetlist.h>
#include <qpopupmenu.h>
#include <qimage.h>
#endif

#include <qmessagebox.h>
#include <qfiledialog.h>
#include <qprinter.h>
#include <qpainter.h>
#include <qgl.h> // include <qglwidget.h>
#include <qdialog.h>
#include <qevent.h> //include <qcontextmenuevent.h>


//////////////////////////////////////////////////////////////////////////////
/**
   Implementation of virtual method of G4VViewer
*/
void G4OpenGLQtViewer::SetView (
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
#ifdef GEANT4_QT_DEBUG
  printf("G4OpenGLQtViewer::SetView ++++++++++++++++++++\n");
#endif
  //   if(!fHDC) return;
  //   if(!fHGLRC) return;
  //   ::wglMakeCurrent(fHDC,fHGLRC);
  //  fWindow->makeCurrent();
  G4OpenGLViewer::SetView ();
#ifdef GEANT4_QT_DEBUG
  printf("G4OpenGLQtViewer::SetView --------------------\n");
#endif
}



//////////////////////////////////////////////////////////////////////////////
/**
   Implementation of virtual method of G4VViewer
*/
void G4OpenGLQtViewer::ShowView (
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
#ifdef GEANT4_QT_DEBUG
  printf("G4OpenGLQtViewer::ShowView  +++++++++++++++++++++\n");
#endif
  glFlush ();
  if (!GLWindow) {
    G4cerr << "Visualization window not defined, please choose one before" << G4endl;
  } else {
#if QT_VERSION < 0x040000
    GLWindow->setActiveWindow();
#else
    GLWindow->activateWindow();
#endif
#ifdef GEANT4_QT_DEBUG
    printf("G4OpenGLQtViewer::ShowView -----------------------\n");
#endif
  }
  //   // Empty the Windows message queue :
  //   MSG event;
  //   while ( ::PeekMessage(&event, NULL, 0, 0, PM_REMOVE) ) {
  //     ::TranslateMessage(&event);
  //     ::DispatchMessage (&event);
  //   }
}



//////////////////////////////////////////////////////////////////////////////
void G4OpenGLQtViewer::CreateGLQtContext (
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
#ifdef GEANT4_QT_DEBUG
  printf("G4OpenGLQtViewer::CreateGLQtContext \n");
#endif
}


//////////////////////////////////////////////////////////////////////////////
void G4OpenGLQtViewer::CreateMainWindow (
 QGLWidget* glWidget
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{

  if(fWindow) return; //Done.
#ifdef GEANT4_QT_DEBUG
  printf("G4OpenGLQtViewer::CreateMainWindow glWidget\n");
#endif

  // launch Qt if not
  G4Qt* interactorManager = G4Qt::getInstance ();
  //  G4UImanager* UI = G4UImanager::GetUIpointer();

  fWindow = glWidget ;
  //  fWindow->makeCurrent();

  // create window
  if (((QApplication*)interactorManager->GetMainInteractor())) {
    // look for the main window
    bool found = false;
#if QT_VERSION < 0x040000
    QWidgetList  *list = QApplication::allWidgets();
    QWidgetListIt it( *list );         // iterate over the widgets
    QWidget * widget;
    while ( (widget=it.current()) != 0 ) {  // for each widget...
      ++it;
      if ((found== false) && (widget->inherits("QMainWindow"))) {
#ifdef GEANT4_QT_DEBUG
        printf("G4OpenGLQtViewer::CreateMainWindow case Qapp exist\n");
#endif
        GLWindow = new QDialog(widget,0,FALSE,Qt::WStyle_Title | Qt::WStyle_SysMenu | Qt::WStyle_MinMax );
        found = true;
      }
    }
    delete list;                      // delete the list, not the widgets
#else
    foreach (QWidget *widget, QApplication::allWidgets()) {
      if ((found== false) && (widget->inherits("QMainWindow"))) {
#ifdef GEANT4_QT_DEBUG
        printf("G4OpenGLQtViewer::CreateMainWindow case Qapp exist\n");
#endif
        GLWindow = new QDialog(widget,Qt::WindowTitleHint | Qt::WindowSystemMenuHint | Qt::WindowMinMaxButtonsHint);
        found = true;
      }
    }
#endif

    if (found==false) {
#ifdef GEANT4_QT_DEBUG
      printf("G4OpenGLQtViewer::CreateMainWindow case Qapp exist, but not found\n");
#endif
      GLWindow = new QDialog();
    }
  } else {
#ifdef GEANT4_QT_DEBUG
    printf("G4OpenGLQtViewer::CreateMainWindow case Qapp exist\n");
#endif
    GLWindow = new QDialog();
  }

#if QT_VERSION < 0x040000
  QHBoxLayout *mainLayout = new QHBoxLayout(GLWindow);
#else
  QHBoxLayout *mainLayout = new QHBoxLayout;
#endif

  mainLayout->addWidget(fWindow);

#if QT_VERSION < 0x040000
  GLWindow->setCaption( tr( "QGl Viewer" ));
#else
  GLWindow->setLayout(mainLayout);
  GLWindow->setWindowTitle(tr("QGl Viewer"));
#endif
  GLWindow->resize(300, 300);
  GLWindow->move(900,300);
  GLWindow->show();
  
  // delete the pointer if close this
  //  GLWindow->setAttribute(Qt::WA_DeleteOnClose);

  QObject ::connect(GLWindow, 
                    SIGNAL(rejected()),
                    this, 
                    SLOT(dialogClosed()));

  WinSize_x = 400;
  WinSize_y = 400;
  if (WinSize_x < fVP.GetWindowSizeHintX ())
    WinSize_x = fVP.GetWindowSizeHintX ();
  if (WinSize_y < fVP.GetWindowSizeHintY ())
    WinSize_y = fVP.GetWindowSizeHintY ();

  if(!fWindow) return;
#ifdef GEANT4_QT_DEBUG
  printf("G4OpenGLQtViewer::CreateMainWindow glWidget END\n");
#endif

  if (!fContextMenu) 
    createPopupMenu();

}

/**  Close the dialog and set the pointer to NULL
 */
void G4OpenGLQtViewer::dialogClosed() {
  GLWindow = NULL;
}


//////////////////////////////////////////////////////////////////////////////
G4OpenGLQtViewer::G4OpenGLQtViewer (
                                    G4OpenGLSceneHandler& scene
                                    )
  :G4VViewer (scene, -1)
  ,G4OpenGLViewer (scene)
  ,fWindow(0)
  ,fContextMenu(0)
  ,fMouseAction(true)
{
#ifdef GEANT4_QT_DEBUG
  printf("G4OpenGLQtViewer::G4OpenGLQtViewer \n");
#endif
}

//////////////////////////////////////////////////////////////////////////////
G4OpenGLQtViewer::~G4OpenGLQtViewer (
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
#ifdef GEANT4_QT_DEBUG
  printf("G4OpenGLQtViewer::~G4OpenGLQtViewer \n");
#endif
  delete fContextMenu;
}


/**
   Create a popup menu for the widget. This menu is activated by right-mouse click
*/
void G4OpenGLQtViewer::createPopupMenu()    {

#if QT_VERSION < 0x040000
  fContextMenu = new QPopupMenu( GLWindow,"All" );
#else
  fContextMenu = new QMenu("All");
#endif
  Q_CHECK_PTR( fContextMenu );

#if QT_VERSION < 0x040000
  // === Mouse menu ===
  QPopupMenu *mMouseAction = new QPopupMenu( fContextMenu );
  Q_CHECK_PTR( mMouseAction );

  QAction *rotate = new QAction("&Rotate scene",CTRL+Key_N,mMouseAction);
  QAction *move =  new QAction("&Move scene",CTRL+Key_M,mMouseAction);
  rotate->addTo(mMouseAction);
  move->addTo(mMouseAction);

  fContextMenu->insertItem( "&Mouse action", mMouseAction);

#else
  // === Mouse menu ===
  QMenu *mMouseAction = fContextMenu->addMenu("&Mouse action");

  QAction *rotate = mMouseAction->addAction("&Rotate scene");
  QAction *move = mMouseAction->addAction("&Move scene");
#endif

  // INIT mMouse
  createRadioAction(rotate,move,SLOT(toggleMouseAction(bool)),1);

#if QT_VERSION < 0x040000
  // === Style Menu ===
  QPopupMenu *mStyle = new QPopupMenu(fContextMenu);

  QPopupMenu *mRepresentation = new QPopupMenu(fContextMenu);

  QAction *polyhedron = new QAction("&Polyhedron",CTRL+Key_P,mRepresentation);
  QAction *nurbs = new QAction("&NURBS",CTRL+Key_N,mRepresentation);
  polyhedron->addTo(mRepresentation);
  nurbs->addTo(mRepresentation);

  mStyle->insertItem("&Representation",mRepresentation);
  fContextMenu->insertItem("&Style",mStyle);

#else
  // === Style Menu ===
  QMenu *mStyle = fContextMenu->addMenu("&Style");

  QMenu *mRepresentation = mStyle->addMenu("&Representation");
  QAction *polyhedron = mRepresentation->addAction("Polyhedron");
  QAction *nurbs = mRepresentation->addAction("NURBS");
#endif
  // INIT mStyle
  G4ViewParameters::RepStyle style;
  style = fVP.GetRepStyle();
  if (style == G4ViewParameters::polyhedron) {
    createRadioAction(polyhedron,nurbs,SLOT(toggleRepresentation(bool)),1);
  } else if (style == G4ViewParameters::nurbs) {
    createRadioAction(polyhedron,nurbs,SLOT(toggleRepresentation(bool)),2);
  } else {
    mRepresentation->clear();
  }


#if QT_VERSION < 0x040000
  // === Drawing Menu ===
  QPopupMenu *mDrawing = new QPopupMenu(fContextMenu);
  fContextMenu->insertItem("&Drawing",mDrawing);

  fDrawingWireframe = new QPopupMenu(mDrawing);
  mDrawing->insertItem("&Wireframe",fDrawingWireframe);

  mDrawing->setCheckable(true);
  fDrawingWireframe->setCheckable(true);

  fDrawingLineRemoval = new QPopupMenu(mDrawing);
  mDrawing->insertItem("&Hidden line removal",fDrawingLineRemoval);
  fDrawingLineRemoval->setCheckable(true);

  fDrawingSurfaceRemoval = new QPopupMenu(mDrawing);
  mDrawing->insertItem("&Hidden surface removal",fDrawingSurfaceRemoval);
  fDrawingSurfaceRemoval->setCheckable(true);

  fDrawingLineSurfaceRemoval = new QPopupMenu(mDrawing);
  mDrawing->insertItem("&Hidden line and surface removal",fDrawingLineSurfaceRemoval);
  fDrawingLineSurfaceRemoval->setCheckable(true);

#else
  // === Drawing Menu ===
  QMenu *mDrawing = mStyle->addMenu("&Drawing");

  fDrawingWireframe = mDrawing->addAction("Wireframe");
  fDrawingWireframe->setCheckable(true);

  fDrawingLineRemoval = mDrawing->addAction("Hidden line removal");
  fDrawingLineRemoval->setCheckable(true);

  fDrawingSurfaceRemoval = mDrawing->addAction("Hidden Surface removal");
  fDrawingSurfaceRemoval->setCheckable(true);

  fDrawingLineSurfaceRemoval = mDrawing->addAction("Hidden line and surface removal");
  fDrawingLineSurfaceRemoval->setCheckable(true);
#endif
  // INIT Drawing
  G4ViewParameters::DrawingStyle d_style;
  d_style = fVP.GetDrawingStyle();
  
#if QT_VERSION < 0x040000
  if (d_style == G4ViewParameters::wireframe) {
    fDrawingWireframe->setItemChecked(0,true);
  } else if (d_style == G4ViewParameters::hlr) {
    fDrawingLineRemoval->setItemChecked(0,true);
  } else if (d_style == G4ViewParameters::hsr) {
    fDrawingSurfaceRemoval->setItemChecked(0,true);
  } else if (d_style == G4ViewParameters::hlhsr) {
    fDrawingLineSurfaceRemoval->setItemChecked(0,true);
  } else {
    mDrawing->clear();
  }
#else
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
#endif

  QObject ::connect(fDrawingWireframe, 
                    SIGNAL(triggered(bool)),
                    this, 
                    SLOT(actionDrawingWireframe()));
  QObject ::connect(fDrawingLineRemoval, 
                    SIGNAL(triggered(bool)),
                    this, 
                    SLOT(actionDrawingLineRemoval()));
  QObject ::connect(fDrawingSurfaceRemoval, 
                    SIGNAL(triggered(bool)),
                    this, 
                    SLOT(actionDrawingSurfaceRemoval()));
  QObject ::connect(fDrawingLineSurfaceRemoval, 
                    SIGNAL(triggered(bool)),
                    this, 
                    SLOT(actionDrawingLineSurfaceRemoval()));


#if QT_VERSION < 0x040000
  QPopupMenu *mBackground = new QPopupMenu(mStyle);
  mStyle->insertItem("&Background color",mBackground);

  QAction *white = new QAction("&White",CTRL+Key_W,mBackground);
  QAction *black =  new QAction("&Black",CTRL+Key_B,mBackground);
  white->addTo(mBackground);
  black->addTo(mBackground);

#else
  QMenu *mBackground = mStyle->addMenu("&Background color");
  QAction *white = mBackground->addAction("White");
  QAction *black = mBackground->addAction("Black");

#endif
  if (background.GetRed() == 1. &&
      background.GetGreen() == 1. &&
      background.GetBlue() == 1.) {
    createRadioAction(white,black,SLOT(toggleBackground(bool)),1);
  } else {
    createRadioAction(white,black,SLOT(toggleBackground(bool)),2);
  }


#if QT_VERSION < 0x040000
  // === Action Menu ===
  QPopupMenu *mActions = new QPopupMenu(fContextMenu);
  fContextMenu->insertItem("&Actions",mActions);

  QAction *controlPanels = new QAction("&Control panels",CTRL+Key_C,mActions);
  QAction *exitG4 =  new QAction("&Exit to G4Vis >",CTRL+Key_E,mActions);
  QAction *createEPS =  new QAction("&Save as ...",CTRL+Key_S,mActions);
  controlPanels->addTo(mActions);
  exitG4->addTo(mActions);
  createEPS->addTo(mActions);

#else
  // === Action Menu ===
  QMenu *mActions = fContextMenu->addMenu("&Actions");
  QAction *controlPanels = mActions->addAction("Control panels");
  QAction *exitG4 = mActions->addAction("Exit to G4Vis >");
  QAction *createEPS = mActions->addAction("Save as ...");
#endif

  QObject ::connect(controlPanels, 
                    SIGNAL(triggered()),
                    this, 
                    SLOT(actionControlPanels()));
  QObject ::connect(exitG4, 
                    SIGNAL(triggered()),
                    this, 
                    SLOT(actionExitG4()));
  QObject ::connect(createEPS, 
                    SIGNAL(triggered()),
                    this,
                    SLOT(actionCreateEPS()));


#if QT_VERSION < 0x040000
  // === Special Menu ===
  QPopupMenu *mSpecial = new QPopupMenu(fContextMenu);
  fContextMenu->insertItem("S&pecial",mSpecial);

  QPopupMenu *mTransparency = new QPopupMenu(mSpecial);
  mSpecial->insertItem("Transparency",mTransparency);

  QAction *transparencyOn = new QAction("&On",CTRL+Key_O,mTransparency);
  QAction *transparencyOff = new QAction("&Off",CTRL+Key_F,mTransparency);
  transparencyOn->addTo(mTransparency);
  transparencyOff->addTo(mTransparency);

#else
  // === Special Menu ===
  QMenu *mSpecial = fContextMenu->addMenu("S&pecial");
  QMenu *mTransparency = mSpecial->addMenu("Transparency");
  QAction *transparencyOn = mTransparency->addAction("On");
  QAction *transparencyOff = mTransparency->addAction("Off");
#endif

  if (transparency_enabled == false) {
    createRadioAction(transparencyOn,transparencyOff,SLOT(toggleTransparency(bool)),2);
  } else if (transparency_enabled == true) {
    createRadioAction(transparencyOn,transparencyOff,SLOT(toggleTransparency(bool)),1);
  } else {
    mSpecial->clear();
  }


#if QT_VERSION < 0x040000
  QPopupMenu *mAntialiasing = new QPopupMenu(mSpecial);
  mSpecial->insertItem("Antialiasing",mAntialiasing);

  QAction *antialiasingOn = new QAction("&On",CTRL+Key_O,mAntialiasing);
  QAction *antialiasingOff = new QAction("&Off",CTRL+Key_F,mAntialiasing);
  antialiasingOn->addTo(mAntialiasing);
  antialiasingOff->addTo(mAntialiasing);

#else
  QMenu *mAntialiasing = mSpecial->addMenu("Antialiasing");
  QAction *antialiasingOn = mAntialiasing->addAction("On");
  QAction *antialiasingOff = mAntialiasing->addAction("Off");
#endif

  if (antialiasing_enabled == false) {
    createRadioAction(antialiasingOn,antialiasingOff,SLOT(toggleAntialiasing(bool)),2);
  } else if (antialiasing_enabled == true) {
    createRadioAction(antialiasingOn,antialiasingOff,SLOT(toggleAntialiasing(bool)),1);
  } else {
    mAntialiasing->clear();
  }

#if QT_VERSION < 0x040000
  QPopupMenu *mHaloing = new QPopupMenu(mSpecial);
  mSpecial->insertItem("Haloing",mHaloing);

  QAction *haloingOn = new QAction("&On",CTRL+Key_O,mHaloing);
  QAction *haloingOff = new QAction("&Off",CTRL+Key_F,mHaloing);
  haloingOn->addTo(mHaloing);
  haloingOff->addTo(mHaloing);
#else
  QMenu *mHaloing = mSpecial->addMenu("Haloing");
  QAction *haloingOn = mHaloing->addAction("On");
  QAction *haloingOff = mHaloing->addAction("Off");
#endif
  if (haloing_enabled == false) {
    createRadioAction(haloingOn,haloingOff,SLOT(toggleHaloing(bool)),2);
  } else if (haloing_enabled == true) {
    createRadioAction(haloingOn,haloingOff,SLOT(toggleHaloing(bool)),1);
  } else {
    mHaloing->clear();
  }

#if QT_VERSION < 0x040000
  QPopupMenu *mAux = new QPopupMenu(mSpecial);
  mSpecial->insertItem("Auxiliairy edges",mAux);

  QAction *auxOn = new QAction("&On",CTRL+Key_O,mAux);
  QAction *auxOff = new QAction("&Off",CTRL+Key_F,mAux);
  auxOn->addTo(mAux);
  auxOff->addTo(mAux);

#else
  QMenu *mAux = mSpecial->addMenu("Auxiliary edges");
  QAction *auxOn = mAux->addAction("On");
  QAction *auxOff = mAux->addAction("Off");
#endif
  if (!fVP.IsAuxEdgeVisible()) {
    createRadioAction(auxOn,auxOff,SLOT(toggleAux(bool)),1);
  } else {
    createRadioAction(auxOn,auxOff,SLOT(toggleAux(bool)),2);
  }


#if QT_VERSION < 0x040000
  QPopupMenu *mFullScreen = new QPopupMenu(mSpecial);
  mSpecial->insertItem("Full screen",mFullScreen);

  QAction *fullOn = new QAction("&On",CTRL+Key_O,mFullScreen);
  QAction *fullOff = new QAction("&Off",CTRL+Key_F,mFullScreen);
  fullOn->addTo(mFullScreen);
  fullOff->addTo(mFullScreen);
#else
  QMenu *mFullScreen = mSpecial->addMenu("Full screen");
  QAction *fullOn = mFullScreen->addAction("On");
  QAction *fullOff = mFullScreen->addAction("Off");
#endif
  createRadioAction(fullOn,fullOff,SLOT(toggleFullScreen(bool)),2);

}

void G4OpenGLQtViewer::manageContextMenuEvent(QContextMenuEvent *e)
{
  if (!GLWindow) {
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
#if QT_VERSION < 0x040000
void G4OpenGLQtViewer::createRadioAction(QAction *action1,QAction *action2, const std::string& method,unsigned int nCheck) {

  if (action1->parent()->inherits("QPopupMenu")){
    ((QPopupMenu*)action1->parent())->setCheckable(true);
  }
  ((QPopupMenu*)action1->parent())->setItemChecked(0,true);
  ((QPopupMenu*)action1->parent())->setItemChecked(1,true);

  if (nCheck ==1)
    ((QPopupMenu*)action1->parent())->setItemChecked(0,true);
  else
    ((QPopupMenu*)action1->parent())->setItemChecked(1,true);
   
  QObject ::connect(action1, SIGNAL(triggered(bool)),action2, SLOT(toggle()));
  QObject ::connect(action2, SIGNAL(triggered(bool)),action1, SLOT(toggle()));

  QObject ::connect(action1, SIGNAL(toggled(bool)),this, method.c_str());
}

#else
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
#endif

/**
   Slot activate when drawing->wireframe menu is set 
 */
void G4OpenGLQtViewer::actionDrawingWireframe() {
  emit toggleDrawingAction(1);
}

/**
   Slot activate when drawing->line removal menu is set 
 */
void G4OpenGLQtViewer::actionDrawingLineRemoval() {
  emit toggleDrawingAction(2);
}

/**
   Slot activate when drawing->surface removal menu is set 
 */
void G4OpenGLQtViewer::actionDrawingSurfaceRemoval() {
  emit toggleDrawingAction(3);
}

/**
   Slot activate when drawing->wireframe menu is set 
 */
void G4OpenGLQtViewer::actionDrawingLineSurfaceRemoval() {
  emit toggleDrawingAction(4);
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
void G4OpenGLQtViewer::toggleDrawingAction(int aAction) {

  G4ViewParameters::DrawingStyle d_style;
  

  if (aAction ==1) {
#if QT_VERSION < 0x040000
    fDrawingWireframe->setItemChecked (0,true);
    fDrawingLineRemoval->setItemChecked (0,false);
    fDrawingSurfaceRemoval->setItemChecked (0,false);
    fDrawingLineSurfaceRemoval->setItemChecked (0,false);
#else
    fDrawingWireframe->setChecked (true);
    fDrawingLineRemoval->setChecked (false);
    fDrawingSurfaceRemoval->setChecked (false);
    fDrawingLineSurfaceRemoval->setChecked (false);
#endif

    d_style = G4ViewParameters::wireframe;

  } else  if (aAction ==2) {
#if QT_VERSION < 0x040000
    fDrawingWireframe->setItemChecked (0,false);
    fDrawingLineRemoval->setItemChecked (0,true);
    fDrawingSurfaceRemoval->setItemChecked (0,false);
    fDrawingLineSurfaceRemoval->setItemChecked (0,false);
#else
    fDrawingWireframe->setChecked (false);
    fDrawingLineRemoval->setChecked (true);
    fDrawingSurfaceRemoval->setChecked (false);
    fDrawingLineSurfaceRemoval->setChecked (false);
#endif

    d_style = G4ViewParameters::hlr;

  } else  if (aAction ==3) {
#if QT_VERSION < 0x040000
    fDrawingWireframe->setItemChecked (0,false);
    fDrawingLineRemoval->setItemChecked (0,false);
    fDrawingSurfaceRemoval->setItemChecked (0,true);
    fDrawingLineSurfaceRemoval->setItemChecked (0,false);
#else
    fDrawingWireframe->setChecked (false);
    fDrawingLineRemoval->setChecked (false);
    fDrawingSurfaceRemoval->setChecked (true);
    fDrawingLineSurfaceRemoval->setChecked (false);
#endif

    d_style = G4ViewParameters::hsr;

  } else  if (aAction ==4) {
#if QT_VERSION < 0x040000
    fDrawingWireframe->setItemChecked (0,false);
    fDrawingLineRemoval->setItemChecked (0,false);
    fDrawingSurfaceRemoval->setItemChecked (0,false);
    fDrawingLineSurfaceRemoval->setItemChecked (0,true);
#else
    fDrawingWireframe->setChecked (false);
    fDrawingLineRemoval->setChecked (false);
    fDrawingSurfaceRemoval->setChecked (false);
    fDrawingLineSurfaceRemoval->setChecked (true);
#endif
    d_style = G4ViewParameters::hlhsr;
  }
  fVP.SetDrawingStyle(d_style);

  updateQWidget();
#ifdef GEANT4_QT_DEBUG
  printf("G4OpenGLQtViewer::toggleDrawingAction\n");
#endif
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

#ifdef GEANT4_QT_DEBUG
  printf("G4OpenGLQtViewer::toggleRepresentation 3%d\n",check);
#endif
  updateQWidget();
#ifdef GEANT4_QT_DEBUG
  printf("G4OpenGLQtViewer::toggleRepresentation 4%d\n",check);
#endif
}

/**
   SLOT Activate by a click on the background menu
@param check : 1 white, 0 black
*/
void G4OpenGLQtViewer::toggleBackground(bool check) {

  //   //I need to revisit the kernel if the background colour changes and
  //   //hidden line removal is enabled, because hlr drawing utilises the
  //   //background colour in its drawing...
  //   // (Note added by JA 13/9/2005) Background now handled in view
  //   // parameters.  A kernel visit is triggered on change of background.
  if (check == 1) {
    ((G4ViewParameters&)this->GetViewParameters()).
      SetBackgroundColour(G4Colour(1.,1.,1.));  // White
  } else {
    ((G4ViewParameters&)this->GetViewParameters()).
      SetBackgroundColour(G4Colour(0.,0.,0.));  // Black
  }
  updateQWidget();
}

/**
   SLOT Activate by a click on the transparency menu
@param check : 1 , 0
*/
void G4OpenGLQtViewer::toggleTransparency(bool check) {
  
  if (check) {
    transparency_enabled = false;
  } else {
    transparency_enabled = true;
  }
  SetNeedKernelVisit (true);
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

  updateQWidget();
#ifdef GEANT4_QT_DEBUG
  printf("G4OpenGLQtViewer::toggleRepresentation %d\n",check);
#endif
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

  updateQWidget();

#ifdef GEANT4_QT_DEBUG
  printf("G4OpenGLQtViewer::toggleRepresentation %d\n",check);
#endif
}

/**
   SLOT Activate by a click on the auxiliaire edges menu
@param check : 1 , 0
*/
void G4OpenGLQtViewer::toggleAux(bool check) {
  if (check) {
    fVP.SetAuxEdgeVisible(false);
  } else {
    fVP.SetAuxEdgeVisible(true);
  }
  SetNeedKernelVisit (true);
  updateQWidget();

#ifdef GEANT4_QT_DEBUG
  printf("G4OpenGLQtViewer::toggleRepresentation %d\n",check);
#endif
}

/**
   SLOT Activate by a click on the full screen menu
@param check : 1 , 0
*/
void G4OpenGLQtViewer::toggleFullScreen(bool check) {
  GLWindow->setWindowState(GLWindow->windowState() ^ Qt::WindowFullScreen);

#ifdef GEANT4_QT_DEBUG
  printf("G4OpenGLQtViewer::toggleRepresentation %d\n",check);
#endif
}

/**
   SLOT Activate by a click on the mouse action menu
   @param check : 1 , 0
*/
void G4OpenGLQtViewer::toggleMouseAction(bool check) {
  if (check) { // rotate scene
    fMouseAction = true;
  } else { // move scene
    fMouseAction = false;
  }

#ifdef GEANT4_QT_DEBUG
  printf("G4OpenGLQtViewer::toggleRepresentation %d\n",check);
#endif
}


void G4OpenGLQtViewer::actionControlPanels() {
#ifdef GEANT4_QT_DEBUG
  printf("G4OpenGLQtViewer::actionControlPanels \n");
#endif
}

void G4OpenGLQtViewer::actionExitG4() {
#ifdef GEANT4_QT_DEBUG
  printf("G4OpenGLQtViewer::actionExitG4() \n");
#endif
}

void G4OpenGLQtViewer::actionCreateEPS() {
  QString filters;
#if QT_VERSION < 0x040000
  QStrList listFormat=QImageIO::outputFormats();
  char *tmp=listFormat.first();
  while (tmp!=0) {
    filters += QString(tmp) + ";;";
    tmp=listFormat.next();
  }
#else
  QList<QByteArray> formats =  QImageWriter::supportedImageFormats ();
  for (int i = 0; i < formats.size(); ++i) {
    filters +=formats.at(i) + ";;";
  }
#endif
  filters += "eps;;";
  filters += "ps;;";
  filters += "pdf";
  QString* selectedFilter = new QString();
#if QT_VERSION < 0x040000
  QString nomFich =  QFileDialog::getSaveFileName ( ".",
                                                    filters,
                                                    GLWindow,
                                                    "Save file dialog",
                                                    tr("Save as ..."),
                                                    selectedFilter ); 
#else
  QString nomFich =  QFileDialog::getSaveFileName ( GLWindow,
                                                    tr("Save as ..."),
                                                    ".",
                                                    filters,
                                                    selectedFilter ); 
#endif
  // bmp jpg jpeg png ppm xbm xpm
  if (nomFich == "") {
    return;
  }
#if QT_VERSION < 0x040000
  nomFich += "."+selectedFilter->lower();
#ifdef GEANT4_QT_DEBUG
  printf("G4OpenGLQtViewer::name %s\n",nomFich.ascii());
#endif
#else
  nomFich += "."+selectedFilter->toLower();
#ifdef GEANT4_QT_DEBUG
  printf("G4OpenGLQtViewer::name %s\n",nomFich.toStdString().c_str());
#endif
#endif
  G4OpenGLQtExportDialog* exportDialog= new G4OpenGLQtExportDialog(GLWindow,nomFich,fWindow->height(),fWindow->width());
  if(  exportDialog->exec()) {

    QImage image;
    //    if ((exportDialog->getWidth() !=fWindow->width()) ||
    //        (exportDialog->getHeight() !=fWindow->height())) {
      
      //      rescaleImage(exportDialog->getWidth(),exportDialog->getHeight());// re-scale image
#ifdef GEANT4_QT_DEBUG
      printf("rescaling\n");
#endif
      QGLWidget* glResized = fWindow;
      fWindow->renderPixmap (exportDialog->getWidth()*2,exportDialog->getHeight()*2 ).save("/Users/laurentgarnier/Desktop/zzz.jpg","jpg");
      QPixmap * pixmap = new QPixmap(fWindow->renderPixmap (exportDialog->getWidth(),exportDialog->getHeight() )) ;
      //      image = pixmap.toImage();
      //      glResized->resize(exportDialog->getWidth()*2,exportDialog->getHeight()*2);
#ifdef GEANT4_QT_DEBUG
      printf("rescaling after\n");
#endif
      //      image = glResized->grabFrameBuffer();
      
      //    } else {
      // image = fWindow->grabFrameBuffer();
      //  }    
    // jpeg format
    if (nomFich.endsWith(".jpg") || 
        nomFich.endsWith(".jpeg")) {
      // grabFrameBuffer() :: Returns an image of the frame buffer. If withAlpha is true the alpha channel is included.
      image.save(nomFich,0,exportDialog->getSliderValue());
#ifdef GEANT4_QT_DEBUG
      printf("saving jpeg quality : %d\n",exportDialog->getSliderValue());
#endif
    } else if (nomFich.endsWith(".eps")) {
      generateEPS(nomFich,exportDialog->getNbColor(),image);
    } else if (nomFich.endsWith(".ps") ||nomFich.endsWith(".pdf")) {
      generatePS_PDF(nomFich,exportDialog->getNbColor(),image);
    } else if (nomFich.endsWith(".tif") ||
               nomFich.endsWith(".tiff") ||
               nomFich.endsWith(".jpg") ||
               nomFich.endsWith(".png") ||
               nomFich.endsWith(".bmp") ||
               nomFich.endsWith(".xpm")) {
      image.save(nomFich,0,exportDialog->getSliderValue());
#ifdef GEANT4_QT_DEBUG
      printf("saving ELSE\n");
#endif
    } else {
      G4cerr << "This version of G4UI Could not generate the selected format" << G4endl;
    }
    
  } else { // cancel selected
    return;
  }
  
#ifdef GEANT4_QT_DEBUG
  printf("G4OpenGLQtViewer::actionCreateEPS() \n");
#endif
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




/**
   Save the current mouse press point
   @param p mouse click point
*/
void G4OpenGLQtViewer::G4MousePressEvent(QPoint p)
{
  lastPos = p;
}

/**
   @param pos_x mouse x position
   @param pos_y mouse y position
   @param mButtons mouse button active
*/

#if QT_VERSION < 0x040000
void G4OpenGLQtViewer::G4MouseMoveEvent(int pos_x, int pos_y,Qt::ButtonState mButtons)
#else
void G4OpenGLQtViewer::G4MouseMoveEvent(int pos_x, int pos_y,Qt::MouseButtons mButtons)
#endif
{
  int dx = pos_x - lastPos.x();
  int dy = pos_y - lastPos.y();
  
  if (fMouseAction) {  // rotate
    if (mButtons & Qt::LeftButton) {
      //phi spin stuff here
      
      G4Vector3D vp = fVP.GetViewpointDirection ().unit ();
      G4Vector3D up = fVP.GetUpVector ().unit ();
      
      G4Vector3D yprime = (up.cross(vp)).unit();
      G4Vector3D zprime = (vp.cross(yprime)).unit();
      
      G4double delta_alpha;
      G4double delta_theta;
      
      if (fVP.GetLightsMoveWithCamera()) {
        delta_alpha = dy;
        delta_theta = -dx;
      } else {
        delta_alpha = -dy;
        delta_theta = dx;
      }    

      delta_alpha *= deg;
      delta_theta *= deg;

      G4Vector3D new_vp = std::cos(delta_alpha) * vp + std::sin(delta_alpha) * zprime;
      
      G4Vector3D new_up;
      if (fVP.GetLightsMoveWithCamera()) {
        new_up = (new_vp.cross(yprime)).unit();
        fVP.SetUpVector(new_up);
      } else {
        new_up = up;
      }
      ////////////////
      // Rotates by fixed azimuthal angle delta_theta.

      G4double cosalpha = new_up.dot (new_vp.unit());
      G4double sinalpha = std::sqrt (1. - std::pow (cosalpha, 2));
      yprime = (new_up.cross (new_vp.unit())).unit ();
      G4Vector3D xprime = yprime.cross (new_up);
      // Projection of vp on plane perpendicular to up...
      G4Vector3D a1 = sinalpha * xprime;
      // Required new projection...
      G4Vector3D a2 =
        sinalpha * (std::cos (delta_theta) * xprime + std::sin (delta_theta) * yprime);
      // Required Increment vector...
      G4Vector3D delta = a2 - a1;
      // So new viewpoint is...
      G4Vector3D viewPoint = new_vp.unit() + delta;

      fVP.SetViewAndLights (viewPoint);
      updateQWidget();
      
    } else if (mButtons & Qt::RightButton) {
      // NEVER DONE BECAUSE OF MOUSE MENU
#ifdef GEANT4_QT_DEBUG
      //       printf("G4OpenGLQtViewer::mouseMoveEvent Right \n");
#endif
      //       setXRotation(xRot + dy/2);
      //       setZRotation(zRot + dx/2);
      //       updateQWidget();
    }
  } else {  // move

    float dx = pos_x - lastPos.x();
    float dy = pos_y - lastPos.y();
    
    G4Point3D stp
      = GetSceneHandler()->GetScene()->GetStandardTargetPoint();
    
    G4Point3D tp = stp + fVP.GetCurrentTargetPoint ();
    
    const G4Vector3D& upVector = fVP.GetUpVector ();
    const G4Vector3D& vpVector = fVP.GetViewpointDirection ();
    
    G4Vector3D unitRight = (upVector.cross (vpVector)).unit();
    G4Vector3D unitUp    = (vpVector.cross (unitRight)).unit();
    
    tp += -dx * unitRight + dy * unitUp;
    fVP.SetCurrentTargetPoint (tp - stp);
    
    updateQWidget();
  }
  lastPos = QPoint(pos_x, pos_y);
}

void G4OpenGLQtViewer::rescaleImage(
 int aWidth
,int aHeight
){
#ifdef GEANT4_QT_DEBUG
  printf("should rescale \n");
#endif
}

/**
   Generate Postscript form image
   @param aFilename : name of file
   @param aInColor : numbers of colors : 1->BW 2->RGB 3->RGB+Alpha
   @param aImage : Image to print
*/
bool G4OpenGLQtViewer::generateEPS (
 QString aFilename
,int aInColor
,QImage aImage
)
{
  // FIXME
#ifdef GEANT4_QT_DEBUG
  printf("saving EPS\n");
#endif

  FILE* fp;

  if ((!aImage.isGrayscale ()) &&(aInColor ==1 )) {
#if QT_VERSION < 0x040000
    aImage.convertDepth(1,Qt::MonoOnly);
#else
    aImage.convertToFormat ( aImage.format(), Qt::MonoOnly);
#endif
  }
  const uchar * pixels = aImage.bits ();
    
  if (pixels == NULL)
    return false;
  
#if QT_VERSION < 0x040000
  fp = fopen (aFilename.ascii(), "w");
#else
  fp = fopen (aFilename.toStdString().c_str(), "w");
#endif
  if (fp == NULL) {
    return false;
  }
  
  fprintf (fp, "%%!PS-Adobe-2.0 EPSF-1.2\n");
#if QT_VERSION < 0x040000
  fprintf (fp, "%%%%Title: %s\n", aFilename.ascii());
#else
  fprintf (fp, "%%%%Title: %s\n", aFilename.toStdString().c_str());
#endif
  fprintf (fp, "%%%%Creator: OpenGL pixmap render output\n");
  fprintf (fp, "%%%%BoundingBox: 0 0 %d %d\n", aImage.width(), aImage.height());
  fprintf (fp, "%%%%EndComments\n");
  fprintf (fp, "gsave\n");
  fprintf (fp, "/bwproc {\n");
  fprintf (fp, "    rgbproc\n");
  fprintf (fp, "    dup length 3 idiv string 0 3 0 \n");
  fprintf (fp, "    5 -1 roll {\n");
  fprintf (fp, "    add 2 1 roll 1 sub dup 0 eq\n");
  fprintf (fp, "    { pop 3 idiv 3 -1 roll dup 4 -1 roll dup\n");
  fprintf (fp, "       3 1 roll 5 -1 roll } put 1 add 3 0 \n");
  fprintf (fp, "    { 2 1 roll } ifelse\n");
  fprintf (fp, "    }forall\n");
  fprintf (fp, "    pop pop pop\n");
  fprintf (fp, "} def\n");
  fprintf (fp, "systemdict /colorimage known not {\n");
  fprintf (fp, "   /colorimage {\n");
  fprintf (fp, "       pop\n");
  fprintf (fp, "       pop\n");
  fprintf (fp, "       /rgbproc exch def\n");
  fprintf (fp, "       { bwproc } image\n");
  fprintf (fp, "   }  def\n");
  fprintf (fp, "} if\n");
  fprintf (fp, "/picstr %d string def\n", aImage.width() * aInColor);
  fprintf (fp, "%d %d scale\n", aImage.width(), aImage.height());
  fprintf (fp, "%d %d %d\n", aImage.width(), aImage.height(), 8);
  fprintf (fp, "[%d 0 0 %d 0 0]\n", aImage.width(), aImage.height());
  fprintf (fp, "{currentfile picstr readhexstring pop}\n");
  fprintf (fp, "false %d\n", aInColor);
  fprintf (fp, "colorimage\n");
  

  int width = aImage.width();
  int height = aImage.height();
  int depth = aImage.depth();
  int size = width*height;
  
  if (depth == 1)
    size = (width+7)/8*height;
  else if (aInColor == 1)
    size = size*3;
  
  int i = 0;
  if (depth == 1) {
    //  To be implemented
    //    QImage::Endian bitOrder = aImage.bitOrder();
    /*    for(int y=0; y < height; y++) {
      const uchar * s = aImage.scanLine(y);
      for(int x=0; x < width; x++) {
        // need to copy bit for bit...
        bool b = (bitOrder == QImage::LittleEndian) ?
          (*(s + (x >> 3)) >> (x & 7)) & 1 :
          (*(s + (x >> 3)) << (x & 7)) & 0x80 ;
        if (b)
          pixel[i >> 3] ^= (0x80 >> (i & 7));
        i++;
      }
      // we need to align to 8 bit here
      i = (i+7) & 0xffffff8;
    }
    */
  } else if (depth == 8) {
#ifdef GEANT4_QT_DEBUG
    printf("has 8 bit\n");
#endif
    for(int y=height-1; y >=0 ; y--) {
      const uchar * s = aImage.scanLine(y);
      for(int x=0; x <width; x++) {
        QRgb rgb = aImage.color(s[x]);
        if (aInColor == 1) {
          fprintf (fp, " %02hx ",(unsigned char)qGray(rgb));
          i++;
        } else {
          fprintf (fp, " %02hx %02hx %02hx",
                   (unsigned char) qRed(rgb),
                   (unsigned char) qGreen(rgb),
                   (unsigned char) qBlue(rgb));
          i += 3;
        }
      }
      fprintf (fp, "\n");
    }
  } else {
#if QT_VERSION < 0x040000
  G4cerr << "GenerateEPS:: No alpha channel image with Qt3. This is only supported with Qt4" << G4endl;
#else
    bool alpha = aImage.hasAlphaChannel();
#ifdef GEANT4_QT_DEBUG
    printf("has else %d alpha %d\n",depth,alpha);
#endif
    for(int y=height-1; y >=0 ; y--) {
      QRgb * s = (QRgb*)(aImage.scanLine(y));
      for(int x=0; x <width; x++) {
        QRgb rgb = (*s++);
        if (alpha && qAlpha(rgb) < 0x40) // 25% alpha, convert to white -
          rgb = qRgb(0xff, 0xff, 0xff);
        if (aInColor == 1) {
          fprintf (fp, " %02hx ",(unsigned char)qGray(rgb));
          i++;
        } else {
          fprintf (fp, " %02hx %02hx %02hx",
                   (unsigned char) qRed(rgb),
                   (unsigned char) qGreen(rgb),
                   (unsigned char) qBlue(rgb));
          i += 3;
        }
      }
      fprintf (fp, "\n");
    } 
#endif

  }

  fprintf (fp, "grestore\n");
  fprintf (fp, "showpage\n");
  fclose (fp);

  return true;
}
/**
   Generate Postscript or PDF form image
   @param aFilename : name of file
   @param aInColor : numbers of colors : 1->BW 2->RGB
   @param aImage : Image to print
*/
bool G4OpenGLQtViewer::generatePS_PDF (
 QString aFilename
,int aInColor
,QImage aImage
)
{
#if QT_VERSION < 0x040000
#ifdef Q_WS_MAC || Q_WS_X11
  QPrinter printer;
  //  printer.setPageSize(pageSize);
  if (aInColor == 1) {
    printer.setColorMode(QPrinter::GrayScale);
  } else {
    printer.setColorMode(QPrinter::Color);
  }

  /* FIXME : I don't know which format it will save...
     if (aFilename.endsWith(".ps")) {
     printer.setOutputFormat(QPrinter::PostScriptFormat);
     } else {
     printer.setOutputFormat(QPrinter::PdfFormat);
     }
  */
  printer.setOutputFileName(aFilename);
  //  printer.setFullPage ( true);
  QPainter paint(&printer);
  paint.drawImage (0,0,aImage );
  paint.end();
#else
  G4cerr << "This fonction is only supported on Mac OsX or X11 with Qt3. Full platform supported with Qt4" << G4endl;
#endif
#else
  QPrinter printer;
  //  printer.setPageSize(pageSize);
  if (aInColor == 1) {
    printer.setColorMode(QPrinter::GrayScale);
  } else {
    printer.setColorMode(QPrinter::Color);
  }

  if (aFilename.endsWith(".ps")) {
#if QT_VERSION > 0x040200
    printer.setOutputFormat(QPrinter::PostScriptFormat);
#endif
  } else {
#if QT_VERSION > 0x040100
    printer.setOutputFormat(QPrinter::PdfFormat);
#endif
  }
#if QT_VERSION > 0x040100
  printer.setOutputFileName(aFilename);
#endif
  //  printer.setFullPage ( true);
  QPainter paint(&printer);
  paint.drawImage (0,0,aImage );
  paint.end();
#endif
  return true;
}

#endif

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
