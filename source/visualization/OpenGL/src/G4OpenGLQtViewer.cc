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
// $Id: G4OpenGLQtViewer.cc 104288 2017-05-23 13:23:23Z gcosmo $
//
// 
// G4OpenGLQtViewer : Class to provide Qt specific
//                     functionality for OpenGL in GEANT4
//
// 27/06/2003 : G.Barrand : implementation (at last !).
// 30/06/2014 : M.Kelsey :  Change QPixmap objects to pointers

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
#include "G4UnitsTable.hh"
#include "G4OpenGLStoredQtViewer.hh"
#include "G4Threading.hh"

#include <CLHEP/Units/SystemOfUnits.h>

#include <typeinfo>

#include <qlayout.h>
#include <qlabel.h>
#include <qdialog.h>
#include <qpushbutton.h>
#include <qprocess.h>
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
#include <qmainwindow.h>
#include <qtablewidget.h>
#include <qheaderview.h>
#include <qscrollarea.h>
#include <qsplitter.h>
#include <qcheckbox.h>
#include <qcursor.h>
#include <qthread.h>

#ifdef G4MULTITHREADED
#include "G4Threading.hh"
#endif

#ifdef G4MULTITHREADED
namespace
{
  G4Mutex mWaitForVisSubThreadQtOpenGLContextMoved = G4MUTEX_INITIALIZER;
  G4Mutex mWaitForVisSubThreadQtOpenGLContextInitialized = G4MUTEX_INITIALIZER;
  G4Condition c1_VisSubThreadQtOpenGLContextInitialized = G4CONDITION_INITIALIZER;
  G4Condition c2_VisSubThreadQtOpenGLContextMoved = G4CONDITION_INITIALIZER;
}
#endif

//////////////////////////////////////////////////////////////////////////////
void G4OpenGLQtViewer::CreateMainWindow (
 QGLWidget* glWidget
 ,const QString& name
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{

  if(fGLWidget) return; //Done.

  fGLWidget = glWidget ;
  //  fGLWidget->makeCurrent();

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
        // INIT size
        fWinSize_x = fVP.GetWindowSizeHintX();
        fWinSize_y = fVP.GetWindowSizeHintY();

        isTabbedView = fUiQt->AddTabWidget((QWidget*)fGLWidget,name);
        QObject::connect(fUiQt->GetViewerTabWidget(),
                         SIGNAL(currentChanged(int)),
                         this,
                         SLOT(currentTabActivated(int)));


      }
      createSceneTreeWidget();
      // activate them
    }
  }

  if (!isTabbedView) { // we have to do a dialog

    QWidget *glDialogWidget = getParentWidget();
    if (glDialogWidget == NULL) {
      return;
    }
    glWidget->setParent(glDialogWidget);
    QHBoxLayout *mainLayout = new QHBoxLayout();
    
    mainLayout->setMargin(0);
    mainLayout->setSpacing(0);   
    mainLayout->addWidget(fGLWidget);
    if (fGLWidget->inherits("QMainWindow")) {
      fGLWidget->setWindowTitle( name);
    }
    glDialogWidget->setLayout(mainLayout);

    
    //useful for MACOSX, we have to compt the menuBar height
    int offset = QApplication::desktop()->height() 
      - QApplication::desktop()->availableGeometry().height();
    
    G4int YPos= fVP.GetWindowAbsoluteLocationHintY(QApplication::desktop()->height());
    if (fVP.GetWindowAbsoluteLocationHintY(QApplication::desktop()->height())< offset) {
      YPos = offset;
    }
    glDialogWidget->resize(getWinWidth(), getWinHeight());
    glDialogWidget->move(fVP.GetWindowAbsoluteLocationHintX(QApplication::desktop()->width()),YPos);
    glDialogWidget->show();
  }

  if(!fGLWidget) return;
  
  if (!fContextMenu) 
    createPopupMenu();

}


//////////////////////////////////////////////////////////////////////////////
G4OpenGLQtViewer::G4OpenGLQtViewer (
 G4OpenGLSceneHandler& scene
)
  :G4VViewer (scene, -1)
  ,G4OpenGLViewer (scene)
  ,fGLWidget(NULL)
  ,fRecordFrameNumber(0)
  ,fMouseOnSceneTree(false)
  ,fContextMenu(0)
  ,fLastPickPoint(-1,-1)
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
  ,fUISceneTreeWidget(NULL)
  ,fUIViewerPropertiesWidget(NULL)
  ,fUIPickInfosWidget(NULL)
  ,fNoKeyPress(true)
  ,fAltKeyPress(false)
  ,fControlKeyPress(false)
  ,fShiftKeyPress(false)
  ,fBatchMode(false)
  ,fCheckSceneTreeComponentSignalLock(false)
  ,fViewerPropertiesTableWidgetIsInit(false)
  ,fSceneTreeComponentTreeWidget(NULL)
  ,fSceneTreeWidget(NULL)
  ,fPVRootNodeCreate(false)
  ,fFilterOutput(NULL)
  ,fNbRotation(0)
  ,fTimeRotation(0)
  ,fTouchableVolumes("Touchables")
  ,fShortcutsDialog(NULL)
  ,fViewerPropertiesTableWidget(NULL)
  ,fPickInfosWidget(NULL)
  ,fPickInfosScrollArea(NULL)
  ,fTreeWidgetInfosIgnoredCommands(0)
  ,fSceneTreeDepthSlider(NULL)
  ,fSceneTreeDepth(1)
  ,fModelShortNameItem(NULL)
  ,fMaxPOindexInserted(-1)
  ,fUiQt(NULL)
  ,fSignalMapperMouse(NULL)
  ,fSignalMapperSurface(NULL)
  ,fSignalMapperPicking(NULL)
  ,fTreeIconOpen(NULL)
  ,fTreeIconClosed(NULL)
  ,fLastExportSliderValue(80)
  ,fLastHighlightColor(G4Color(0,0,0,0))
  ,fLastHighlightName(0)
  ,fIsDeleting(false)
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
  fSignalMapperMouse = new QSignalMapper(this);
  fSignalMapperSurface = new QSignalMapper(this);

  // Set default path and format
  fFileSavePath = QDir::currentPath();

  // add available export format
  QList<QByteArray> formats =  QImageWriter::supportedImageFormats ();
  for (int i = 0; i < formats.size(); ++i) {
    addExportImageFormat(formats.at(i).data());
  }

  const char * const icon1[]={
    /* columns rows colors chars-per-pixel */
          "20 20 34 1",
     "  c None",
     ". c #7C7C7C7C7C7C",
     "X c #7D7D7D7D7D7D",
     "o c #828282828282",
     "O c #838383838383",
     "+ c #848484848484",
     "@ c #858585858585",
     "# c #878787878787",
     "$ c #888888888888",
     "% c #8B8B8B8B8B8B",
     "& c #8C8C8C8C8C8C",
     "* c #8F8F8F8F8F8F",
     "= c #909090909090",
     "- c #919191919191",
     "; c #999999999999",
     ": c #9D9D9D9D9D9D",
     "> c #A2A2A2A2A2A2",
     ", c #A3A3A3A3A3A3",
     "< c #A5A5A5A5A5A5",
     "1 c #A6A6A6A6A6A6",
     "2 c #B3B3B3B3B3B3",
     "3 c #B6B6B6B6B6B6",
     "4 c #C2C2C2C2C2C2",
     "5 c #C6C6C6C6C6C6",
     "6 c #CACACACACACA",
     "7 c #CFCFCFCFCFCF",
     "8 c #D0D0D0D0D0D0",
     "9 c #D4D4D4D4D4D4",
     "0 c #D7D7D7D7D7D7",
     "q c #DEDEDEDEDEDE",
     "w c #E0E0E0E0E0E0",
     "e c #E7E7E7E7E7E7",
     "r c #F4F4F4F4F4F4",
     "t c #F7F7F7F7F7F7",
     "                               ",
     "                               ",
     "                               ",
     "                               ",
     "                               ",
     "                               ",
     "    =========>                 ",
     "    7&X+++Oo<e                 ",
     "     2o+@@+-8                  ",
     "     w;.#@+3                   ",
     "      4$o@:q                   ",
     "      r1X%5                    ",
     "       9*,t                    ",
     "        60                     ",
     "                               ",
     "                               ",
     "                               ",
     "                               ",
     "                               ",
     "                               "
  };
  const char * const icon2[]={
    "20 20 68 1",
    "  c None",
    ". c #5F5F10102323",
    "X c #40405F5F1010",
    "o c #696963632E2E",
    "O c #101019194C4C",
    "+ c #101023237070",
    "@ c #70702D2D6363",
    "# c #73732D2D6464",
    "$ c #79792E2E6767",
    "% c #19194C4C5353",
    "& c #2D2D63636161",
    "* c #2E2E61617070",
    "= c #6F6F6E6E4343",
    "- c #707065655F5F",
    "; c #727279795454",
    ": c #535341417070",
    "> c #797954547979",
    ", c #434361617474",
    "< c #414170707070",
    "1 c #686869696363",
    "2 c #6C6C69696363",
    "3 c #656567676F6F",
    "4 c #69696F6F6E6E",
    "5 c #747465656767",
    "6 c #757562626C6C",
    "7 c #70706C6C6969",
    "8 c #616174746565",
    "9 c #656573736969",
    "0 c #616174746969",
    "q c #707075756262",
    "w c #797970706565",
    "e c #636361617474",
    "r c #67676F6F7272",
    "t c #727261617070",
    "y c #616170707070",
    "u c #6F6F72727979",
    "i c #67676E6ED1D1",
    "p c #808080808080",
    "a c #828282828282",
    "s c #838383838383",
    "d c #848484848484",
    "f c #858585858585",
    "g c #868686868686",
    "h c #888888888888",
    "j c #8A8A8A8A8A8A",
    "k c #8D8D8D8D8D8D",
    "l c #8F8F8F8F8F8F",
    "z c #909090909090",
    "x c #949494949494",
    "c c #9C9C9C9C9C9C",
    "v c #9F9F9F9F9F9F",
    "b c #A2A2A2A2A2A2",
    "n c #AEAEAEAEAEAE",
    "m c #B7B7B7B7B7B7",
    "M c #C7C7C7C7C7C7",
    "N c #C9C9C9C9C9C9",
    "B c #D1D1D1D1D1D1",
    "V c #D4D4D4D4D4D4",
    "C c #D9D9D9D9D9D9",
    "Z c #E0E0E0E0E0E0",
    "A c #E2E2E2E2E2E2",
    "S c #EEEEEEEEEEEE",
    "D c #F0F0F0F0F0F0",
    "F c #F5F5F5F5F5F5",
    "G c #F6F6F6F6F6F6",
    "H c #F9F9F9F9F9F9",
    "J c #FCFCFCFCFCFC",
    "K c #FDFDFDFDFDFD",
    "                    ",
    "                    ",
    "                    ",
    "                    ",
    "                    ",
    "     bC             ",
    "     zjnD           ",
    "     ldjjMK         ",
    "     zdhdjcA        ",
    "     zddhdddVK      ",
    "     zghdalBH       ",
    "     zghamSK        ",
    "     lubZH          ",
    "     xMF            ",
    "     G              ",
    "                    ",
    "                    ",
    "                    ",
    "                    ",
    "                    ",
    
  };
  
  const char * const search[]  = {
    /* columns rows colors chars-per-pixel */
    "19 19 8 1",
    "  c #5C5C5C",
    ". c #7D7D7D",
    "X c #9B9B9B",
    "o c #C3C3C3",
    "O c None",
    "+ c #000000",
    "@ c #000000",
    "# c None",
    /* pixels */
    "OOOOOOOOOOOOOOOOOOO",
    "OOOOOOOOOOOOOOOOOOO",
    "OOOOOOOo.  .oOOOOOO",
    "OOOOOOX      XOOOOO",
    "OOOOOo  XOOX  oOOOO",
    "OOOOO. XOOOOX .OOOO",
    "OOOOO  OOOOOO  OOOO",
    "OOOOO  OOOOOO  OOOO",
    "OOOOO. XOOOOo .OOOO",
    "OOOOOo  oOOo  oOOOO",
    "OOOOOOX       XOOOO",
    "OOOOOOOo.  .   XOOO",
    "OOOOOOOOOOOOO.  XOO",
    "OOOOOOOOOOOOOO. XOO",
    "OOOOOOOOOOOOOOOoOOO",
    "OOOOOOOOOOOOOOOOOOO",
    "OOOOOOOOOOOOOOOOOOO",
    "OOOOOOOOOOOOOOOOOOO",
    "OOOOOOOOOOOOOOOOOOO"
  };
  
  fSearchIcon = new QPixmap(search);
  fTreeIconOpen = new QPixmap(icon1);
  fTreeIconClosed = new QPixmap(icon2);

}

//////////////////////////////////////////////////////////////////////////////
G4OpenGLQtViewer::~G4OpenGLQtViewer (
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  fIsDeleting = true;
  
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
  }

  // Delete the open/close icons
  delete fTreeIconOpen;
  delete fTreeIconClosed;

  G4cout <<removeTempFolder().toStdString().c_str() <<G4endl;
}


//
//   Create a popup menu for the widget. This menu is activated by right-mouse click
//
void G4OpenGLQtViewer::createPopupMenu()    {

  fContextMenu = new QMenu("All");

  QMenu *mMouseAction = fContextMenu->addMenu("&Mouse actions");

  fMouseRotateAction = mMouseAction->addAction("Rotate", fSignalMapperMouse, SLOT(map()));
  fMouseMoveAction = mMouseAction->addAction("Move", fSignalMapperMouse, SLOT(map()));
  fMousePickAction = mMouseAction->addAction("Pick", fSignalMapperMouse, SLOT(map()));
  fMouseZoomOutAction = mMouseAction->addAction("Zoom out", fSignalMapperMouse, SLOT(map()));
  fMouseZoomInAction = mMouseAction->addAction("Zoom in", fSignalMapperMouse, SLOT(map()));
  QAction *shortcutsAction = mMouseAction->addAction("Show shortcuts");

  fMouseRotateAction->setCheckable(true);
  fMouseMoveAction->setCheckable(true);
  fMousePickAction->setCheckable(true);
  fMouseZoomOutAction->setCheckable(true);
  fMouseZoomInAction->setCheckable(true);
  shortcutsAction->setCheckable(false);

  connect(fSignalMapperMouse, SIGNAL(mapped(int)),this, SLOT(toggleMouseAction(int)));
  fSignalMapperMouse->setMapping(fMouseRotateAction,1);
  fSignalMapperMouse->setMapping(fMouseMoveAction,2);
  fSignalMapperMouse->setMapping(fMousePickAction,3);
  fSignalMapperMouse->setMapping(fMouseZoomOutAction,4);
  fSignalMapperMouse->setMapping(fMouseZoomInAction,5);

  QObject::connect(shortcutsAction, 
                    SIGNAL(triggered(bool)),
                    this, 
                    SLOT(showShortcuts()));

  // === Style Menu ===
  QMenu *mStyle = fContextMenu->addMenu("&Style");

  QMenu *mProjection = mStyle->addMenu("&Projection");

  fProjectionOrtho = mProjection->addAction("Orthographic", fSignalMapperSurface, SLOT(map()));
  fProjectionPerspective = mProjection->addAction("Persepective", fSignalMapperSurface, SLOT(map()));

 // INIT mProjection
  if (fVP.GetFieldHalfAngle() == 0) {
    createRadioAction(fProjectionOrtho, fProjectionPerspective,SLOT(toggleProjection(bool)),1);
  } else {
    createRadioAction(fProjectionOrtho, fProjectionPerspective,SLOT(toggleProjection(bool)),2);
  }

  // === Drawing Menu ===
  QMenu *mDrawing = mStyle->addMenu("&Drawing");

  fDrawingWireframe = mDrawing->addAction("Wireframe", fSignalMapperSurface, SLOT(map()));

  fDrawingLineRemoval = mDrawing->addAction("Hidden line removal", fSignalMapperSurface, SLOT(map()));

  fDrawingSurfaceRemoval = mDrawing->addAction("Hidden Surface removal", fSignalMapperSurface, SLOT(map()));

  fDrawingLineSurfaceRemoval = mDrawing->addAction("Hidden line and surface removal", fSignalMapperSurface, SLOT(map()));

  fDrawingWireframe->setCheckable(true);
  fDrawingLineRemoval->setCheckable(true);
  fDrawingSurfaceRemoval->setCheckable(true);
  fDrawingLineSurfaceRemoval->setCheckable(true);

  connect(fSignalMapperSurface, SIGNAL(mapped(int)),this, SLOT(toggleSurfaceAction(int)));
  fSignalMapperSurface->setMapping(fDrawingWireframe,1);
  fSignalMapperSurface->setMapping(fDrawingLineRemoval,2);
  fSignalMapperSurface->setMapping(fDrawingSurfaceRemoval,3);
  fSignalMapperSurface->setMapping(fDrawingLineSurfaceRemoval,4);


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
  QAction *movieParameters = mActions->addAction("Save as movie...");
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
  if (!fGLWidget) {
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
      text += "CTRL + Click and move mouse to zoom in/out \n";
      text += "SHIFT + Click and move mouse to change camera point of view \n";
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
    togglePicking();
  } else  if (aAction == 4) {
    fUiQt->SetIconZoomOutSelected();
  } else  if (aAction == 5) {
    fUiQt->SetIconZoomInSelected();
  }

  updateQWidget();
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


void G4OpenGLQtViewer::togglePicking() {
  // FIXME : Not the good way to do, we should handle the multiple cases of Icon/ContextMenu and CheckBox in a better way
  if (fUiQt) {
    if (!fVP.IsPicking()) {
      fUiQt->SetIconPickSelected();
    } else {
      fUiQt->SetIconRotateSelected();
    }
  }
  
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI != NULL)  {
    if (!fVP.IsPicking()) {
      UI->ApplyCommand(std::string("/vis/viewer/set/picking true"));
    } else {
      UI->ApplyCommand(std::string("/vis/viewer/set/picking false"));
    }
  }

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
  if (check != fGLWidget->isFullScreen()) { //toggle
    fGLWidget->setWindowState(fGLWidget->windowState() ^ Qt::WindowFullScreen);
  }
}


void G4OpenGLQtViewer::savePPMToTemp() {
  if (fMovieTempFolderPath == "") {
    return;
  }
  QGLWidget* qGLW = dynamic_cast<QGLWidget*> (fGLWidget) ;
  if (! qGLW) {
    return;
  }
  QString fileName ="Test"+QString::number(fRecordFrameNumber)+".ppm";
  QString filePath =fMovieTempFolderPath+fileName;

  QImage image;
  image = qGLW->grabFrameBuffer();
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
  for (unsigned int i = 0; i < fExportImageFormatVector.size(); ++i) {
    filters += QString("*.") + fExportImageFormatVector.at(i).c_str() + ";;";
  }

  QString* selectedFormat = new QString(fDefaultExportImageFormat.c_str());
  QString qFilename;
  qFilename =  QFileDialog::getSaveFileName ( fGLWidget,
                                         tr("Save as ..."),
                                         fFileSavePath,
                                         filters,
                                         selectedFormat );

  
  std::string name = qFilename.toStdString().c_str();

  // bmp jpg jpeg png ppm xbm xpm
  if (name.empty()) {
    return;
  }

  fFileSavePath = QFileInfo(qFilename).path();
  
  std::string format = selectedFormat->toLower().toStdString().c_str();
  
  // set the format to current
  fExportImageFormat = format.substr(format.find_last_of(".") + 1);

  std::string filename = name;
  std::string extension = "";
  if (name.find_last_of(".") != std::string::npos) {
    filename = name.substr(0,name.find_last_of(".") + 1);
    extension = name.substr(name.find_last_of(".") + 1);
  } else {
    extension = fExportImageFormat;
  }
  
  filename+= "."+ extension;
  
  if (!setExportFilename(filename.c_str(),0)) {
    return;
  }
  
  G4OpenGLQtExportDialog* exportDialog= new G4OpenGLQtExportDialog(fGLWidget,format.c_str(),fGLWidget->height(),fGLWidget->width());
  if(  exportDialog->exec()) {

    if ((exportDialog->getWidth() !=fGLWidget->width()) ||
        (exportDialog->getHeight() !=fGLWidget->height())) {
      setExportSize(exportDialog->getWidth(),exportDialog->getHeight());
      
    }
    if (fExportImageFormat == "eps") {
      fVectoredPs = exportDialog->getVectorEPS();
    } else if (fExportImageFormat == "ps") {
      fVectoredPs = true;
    }
    fLastExportSliderValue = exportDialog->getSliderValue();

    if (exportImage(filename)) {
      // set the default format to current
      fDefaultExportImageFormat = format;
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
  const QColor color = QColor(QColorDialog::getRgba (QColor(Qt::black).rgba(),&a,fGLWidget));
#else
  const QColor color =
    QColorDialog::getColor(Qt::black,
                           fGLWidget,
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
  const QColor color = QColor(QColorDialog::getRgba (QColor(Qt::yellow).rgba(),&a,fGLWidget));
#else
  const QColor& color =
    QColorDialog::getColor(Qt::yellow,
                           fGLWidget,
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
  const QColor color = QColor(QColorDialog::getRgba (QColor(Qt::white).rgba(),&a,fGLWidget));
#else
  const QColor& color =
    QColorDialog::getColor(Qt::white,
                           fGLWidget,
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
    fMovieParametersDialog= new G4OpenGLQtMovieDialog(this,fGLWidget);
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



void G4OpenGLQtViewer::FinishView()
{
 /* From Apple doc: 
  CGLFlushDrawable : Copies the back buffer of a double-buffered context to the front buffer.
  If the backing store attribute is set to false, the buffers can be exchanged rather than copied
 */
  glFlush ();

  // L. Garnier 10/2009 : Not necessary and cause problems on mac OS X 10.6
  //  fGLWidget->swapBuffers ();
}

/**
   Save the current mouse press point
   @param p mouse click point
*/
void G4OpenGLQtViewer::G4MousePressEvent(QMouseEvent *evnt)
{
  if (evnt->button() == Qt::RightButton) {
    return;
  }
  if ((evnt->button() & Qt::LeftButton) && (! (evnt->modifiers() & Qt::ControlModifier ))){
    fGLWidget->setMouseTracking(true);
    fAutoMove = false; // stop automove
    fLastPos1 = evnt->pos();
    fLastPos2 = fLastPos1;
    fLastPos3 = fLastPos2;
    fLastEventTime->start();
    if (fUiQt != NULL) {

      if (fUiQt->IsIconZoomInSelected()) {  // zoomIn
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

      } else if (fUiQt->IsIconRotateSelected() ) {
 
        if (fShiftKeyPress) { // move
          fGLWidget->setCursor(QCursor(Qt::SizeAllCursor));

        } else {  // rotate
          fGLWidget->setCursor(QCursor(Qt::ClosedHandCursor));
        }
      } else if (fUiQt->IsIconMoveSelected()) {
        fGLWidget->setCursor(QCursor(Qt::SizeAllCursor));
      } else if (fUiQt->IsIconPickSelected()) {
        fGLWidget->setCursor(QCursor(Qt::PointingHandCursor));
      }
    }
  }
}

/**
 */
void G4OpenGLQtViewer::G4MouseReleaseEvent(QMouseEvent *evnt)
{
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);

  // factorX == factorY
  double factorX =  ((double)viewport[2]/fGLWidget->width());
  double factorY =  ((double)viewport[3]/fGLWidget->height());
  fSpinningDelay = fLastEventTime->elapsed();
  QPoint delta = (fLastPos3-fLastPos1)*factorX;
  
  // reset cursor state
  fGLWidget->setCursor(QCursor(Qt::ArrowCursor));
  
  if (fVP.IsPicking()){  // pick
    if ((delta.x() != 0) || (delta.y() != 0)) {
      return;
    }
    updatePickInfosWidget(evnt->pos().x()*factorX,evnt->pos().y()*factorY);
    
  } else if (fSpinningDelay < fLaunchSpinDelay ) {
    if ((delta.x() == 0) && (delta.y() == 0)) {
      return;
    }

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
        // prevent from closing widget when rotating (cause a crash)
        if (fIsDeleting) {
          return;
        }
        
        if (rotate) {  // rotate
          if (fNoKeyPress) {
            rotateQtScene(((float)delta.x())/correctionFactor,((float)delta.y())/correctionFactor);
          } else if (fAltKeyPress) {
            rotateQtSceneToggle(((float)delta.x())/correctionFactor,((float)delta.y())/correctionFactor);
          }
          
        } else if (move) {  // move
          moveScene(-((float)delta.x())/correctionFactor,-((float)delta.y())/correctionFactor,0,true);
        }
      }
      ((QApplication*)G4Qt::getInstance ())->processEvents();
    }
  }
  fGLWidget->setMouseTracking(false);

}


void G4OpenGLQtViewer::G4MouseDoubleClickEvent()
{
  fGLWidget->setMouseTracking(true);
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

  bool move = false;
  if (fUiQt != NULL) {
    if (fUiQt->IsIconMoveSelected()) {  // move
      move = true;
    }
  }
  if (!move) {  // rotate, pick, zoom...
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
    fGLWidget->setCursor(QCursor(Qt::SizeAllCursor));

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
    fGLWidget->setCursor(QCursor(Qt::ClosedHandCursor));

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


void G4OpenGLQtViewer::G4keyReleaseEvent (QKeyEvent *)
{
  fGLWidget->setCursor(QCursor(Qt::ArrowCursor));
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
    return "ppmtompeg is needed to encode in video format. It is available here: http://netpbm.sourceforge.net ";
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

/** 
  Export image. Try to get the format according to the file extention.
  If not present, the last one choosen by /vis/ogl/set/exportFormat
  If not, will take the default format : eps
  Best format actually available is pdf (vectored and allow transparency)
  If name is not set, it will take the default name value given by /vis/ogl/set/printFilename
 */
bool G4OpenGLQtViewer::exportImage(std::string name, int width, int height) {

  QGLWidget* qGLW = dynamic_cast<QGLWidget*> (fGLWidget) ;
  if (! qGLW) {
    return false;
  }
  // If there is already an extention
  bool increaseFileNumber = true;
  // if
  if (name.size() != name.substr(name.find_last_of(".") + 1).size()) {
    increaseFileNumber = false;
  }
  if (! setExportFilename(name,increaseFileNumber)) {
    return false;
  }
  if ((width !=-1) && (height != -1)) {
    setExportSize(width, height);
  }
  // first, try to do it with generic function
  if (G4OpenGLViewer::exportImage(name, width, height)) {
    return true;

  // Then try Qt saving functions
  } else {
    QImage image;
    image = qGLW->grabFrameBuffer();

    bool res = image.save(QString(getRealPrintFilename().c_str()),0,fLastExportSliderValue);
    
    if (!res) {
      G4cerr << "Error saving file... " << getRealPrintFilename().c_str() << G4endl;
      return false;
    } else {
      G4cout << "File " << getRealPrintFilename().c_str() << " size: " << fGLWidget->width() << "x" << fGLWidget->height() << " has been saved " << G4endl;
      fExportFilenameIndex++;
    }
  }
  return true;
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
  fprintf (fp,"PATTERN		I\n");
  fprintf (fp,"OUTPUT		%s\n",getSaveFileName().toStdString().c_str());
  fprintf (fp,"\n");
  fprintf (fp,"# You must specify the type of the input files.  The choices are:\n");
  fprintf (fp,"#    YUV, PPM, JMOVIE, Y, JPEG, PNM\n");
  fprintf (fp,"#	(must be upper case)\n");
  fprintf (fp,"#\n");
  fprintf (fp,"BASE_FILE_FORMAT	PPM\n");
  fprintf (fp,"\n");
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
  fprintf (fp,"GOP_SIZE	1\n");
  fprintf (fp,"\n");
  fprintf (fp,"# number of slices in a frame\n");
  fprintf (fp,"#\n");
  fprintf (fp,"# 1 is a good number.  another possibility is the number of macroblock rows\n");
  fprintf (fp,"# (which is the height divided by 16)\n");
  fprintf (fp,"#\n");
  fprintf (fp,"SLICES_PER_FRAME	1\n");
  fprintf (fp,"PIXEL HALF");
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
  QDialog* dialog = NULL;
  // create window
  if (((QApplication*)interactorManager->GetMainInteractor())) {
    // look for the main window
    QWidgetList wl = QApplication::allWidgets();
    QWidget *widget = NULL;
    for (int i=0; i < wl.size(); i++) {
      widget = wl.at(i);
      if ((found== false) && (widget->inherits("QMainWindow"))) {
        dialog = new QDialog(widget,Qt::WindowTitleHint | Qt::WindowSystemMenuHint | Qt::WindowMinMaxButtonsHint);
        found = true;
      }
    }
    
    if (found==false) {
      dialog = new QDialog();
    }
  } else {
    dialog= new QDialog();
  }
  if (found) {
    return dialog;
  } else {
    return NULL;
  }
}


void G4OpenGLQtViewer::createSceneTreeWidget() {
  fUISceneTreeWidget = fUiQt->GetSceneTreeWidget();
  
  if (!fUISceneTreeWidget) {
    return;
  }
  
  // do not remove previous widgets, hide them!
  QLayoutItem * wItem;
  bool found = false;
  if (fUISceneTreeWidget->layout()->count() ) {
    for(int idx = 0; idx < fUISceneTreeWidget->layout()->count(); idx++){
      wItem = fUISceneTreeWidget->layout()->itemAt(idx);
      if (fSceneTreeWidget) {
        if(dynamic_cast<QWidget *>(wItem->widget())) {
          if (wItem->widget()->windowTitle() == fSceneTreeWidget->windowTitle()) {
            wItem->widget()->show();
            found = true;
          } else {
            wItem->widget()->hide();
          }
        }
      } else {
        wItem->widget()->hide();
      }
    }
  }

  if (!found) {
    // initialize scene tree / viewer properties / picking
    fSceneTreeWidget = new QWidget();
    QVBoxLayout* layoutSceneTree = new QVBoxLayout();
    fSceneTreeWidget->setStyleSheet ("padding: 0px ");
    
    fSceneTreeWidget->setLayout(layoutSceneTree);
    fSceneTreeWidget->layout()->setContentsMargins(5,5,5,5);
    fSceneTreeWidget->setWindowTitle(QString(GetName().data()));

    if (fUISceneTreeWidget != NULL) {
      fUISceneTreeWidget->layout()->addWidget(fSceneTreeWidget);
    }
    
    // not available for Immediate mode
    if (dynamic_cast<G4OpenGLStoredQtViewer*> (this)) {
      createSceneTreeComponent();
    }
  }
}


void G4OpenGLQtViewer::createSceneTreeComponent(){
  
  QLayout* vLayout = fSceneTreeWidget->layout();
  
  // Search line
  QWidget* coutButtonWidget = new QWidget();
  QHBoxLayout* layoutCoutTBButtons = new QHBoxLayout();

  fFilterOutput = new QLineEdit();
  fFilterOutput->setToolTip("Filter output by...");
  fFilterOutput->setStyleSheet ("padding: 0px ");

  QPixmap* searchIcon = fUiQt->getSearchIcon();
#if QT_VERSION > 0x050100
  fFilterOutput->addAction(*searchIcon,QLineEdit::TrailingPosition);
  fFilterOutput->setStyleSheet ("border-radius:7px;");
#else
  QPushButton *coutTBFilterButton = new QPushButton();
  coutTBFilterButton->setIcon(*searchIcon);
  coutTBFilterButton->setStyleSheet ("padding-left: 0px; border:0px;");
  fFilterOutput->setStyleSheet ("padding-right: 0px;");
#endif
  layoutCoutTBButtons->addWidget(fFilterOutput);
  
#if QT_VERSION <= 0x050100
  layoutCoutTBButtons->addWidget(coutTBFilterButton);
#endif

  coutButtonWidget->setLayout(layoutCoutTBButtons);
  vLayout->addWidget(coutButtonWidget);
  
  // reduce margins
  vLayout->setContentsMargins(0,0,0,0);


  fSceneTreeComponentTreeWidget = new QTreeWidget();
  fSceneTreeComponentTreeWidget->setSelectionMode(QAbstractItemView::SingleSelection);
  fSceneTreeComponentTreeWidget->setHeaderLabel ("Scene tree : "+QString(GetName().data()));
  fSceneTreeComponentTreeWidget->setColumnHidden (1,true);  // copy number
  fSceneTreeComponentTreeWidget->setColumnHidden (2,true);  // PO index
  fSceneTreeComponentTreeWidget->setColumnHidden (3,true);  // Informations
  //   data(0) : POindex
  //   data(1) : copy number
  //   data(2) : g4color

  vLayout->addWidget(fSceneTreeComponentTreeWidget);

  connect(fSceneTreeComponentTreeWidget,SIGNAL(itemChanged(QTreeWidgetItem*, int)),SLOT(sceneTreeComponentItemChanged(QTreeWidgetItem*, int)));
  connect(fSceneTreeComponentTreeWidget,SIGNAL(itemSelectionChanged ()),SLOT(sceneTreeComponentSelected()));
  connect(fSceneTreeComponentTreeWidget,SIGNAL(itemDoubleClicked ( QTreeWidgetItem*, int)),SLOT(changeColorAndTransparency( QTreeWidgetItem*, int)));

  
  // Depth slider
  QWidget *helpWidget = new QWidget();
  QHBoxLayout *helpLayout = new QHBoxLayout();

  QWidget* depthWidget = new QWidget();
  QWidget *showBox = new QWidget(depthWidget);
  QHBoxLayout *showBoxLayout = new QHBoxLayout();

  // reduce margins
  showBoxLayout->setContentsMargins(5,5,5,5);

  QLabel *zero = new QLabel();
  zero->setText("Show all");          
  QLabel *one = new QLabel();
  one->setText("Hide all");
  fSceneTreeDepthSlider = new QSlider ( Qt::Horizontal);
  fSceneTreeDepthSlider->setMaximum (1000);
  fSceneTreeDepthSlider->setMinimum (0);
  fSceneTreeDepthSlider->setTickPosition(QSlider::TicksAbove);
  // set a minimum size
  fSceneTreeDepthSlider->setMinimumWidth (40);

  showBoxLayout->addWidget(zero);
  showBoxLayout->addWidget(fSceneTreeDepthSlider);
  showBoxLayout->addWidget(one);

  showBox->setLayout(showBoxLayout);

  helpLayout->addWidget(showBox);
  helpWidget->setLayout(helpLayout);
  helpLayout->setContentsMargins(0,0,0,0);

  vLayout->addWidget(helpWidget);
  
  connect( fSceneTreeDepthSlider, SIGNAL( valueChanged(int) ), this, SLOT( changeDepthInSceneTree(int) ) );
  connect( fFilterOutput, SIGNAL( textEdited ( const QString &) ), this, SLOT(changeSearchSelection()));
  fTreeItemModels.clear();
  
  fPVRootNodeCreate = false;
  
  fMaxPOindexInserted = -1;
  

}


void G4OpenGLQtViewer::createViewerPropertiesWidget() {
  
  // Get the pointer to the Viewer Properties widget
  fUIViewerPropertiesWidget = fUiQt->GetViewerPropertiesWidget();

  if (!fUIViewerPropertiesWidget) {
    return;
  }

  // remove previous widgets
  QLayoutItem * wItem;
  if (fUIViewerPropertiesWidget->layout()->count()) {
    while ((wItem = fUIViewerPropertiesWidget->layout()->takeAt(0)) != 0) {
      delete wItem->widget();
      delete wItem;
    }
  }

  // add properties
  QGroupBox *groupBox = new QGroupBox();
  groupBox->setTitle(GetName().data());
  QVBoxLayout *vbox = new QVBoxLayout;
  
  // add properties content
  fViewerPropertiesTableWidget = new QTableWidget();

  QSizePolicy vPolicy = fViewerPropertiesTableWidget->sizePolicy();
  vPolicy.setVerticalStretch(4);
  
  vbox->addWidget(fViewerPropertiesTableWidget);
  groupBox->setLayout(vbox);
  fUIViewerPropertiesWidget->layout()->addWidget(groupBox);

  connect(fViewerPropertiesTableWidget, SIGNAL(itemChanged(QTableWidgetItem*)),this, SLOT(tableWidgetViewerSetItemChanged(QTableWidgetItem *)));
  
  updateViewerPropertiesTableWidget();

  QDialog* dial = static_cast<QDialog*> (fUIViewerPropertiesWidget->parent());
  if (dial) {
    // change name
    dial->setWindowTitle(QString("Viewer properties - ")+GetName());
  }
}


void G4OpenGLQtViewer::createPickInfosWidget(){
  
  // Get the pointer to the Pick infos widget
  fUIPickInfosWidget = fUiQt->GetPickInfosWidget();
  
  if (!fUIPickInfosWidget) {
    return;
  }
  
  // remove previous widgets
  QLayoutItem * wItem;
  if (fUIPickInfosWidget->layout()->count()) {
    while ((wItem = fUIPickInfosWidget->layout()->takeAt(0)) != 0) {
      delete wItem->widget();
      delete wItem;
    }
  }
  
  QGroupBox *groupBox = new QGroupBox("");
  QVBoxLayout *vbox = new QVBoxLayout;

  // add picking infos
  QWidget *pickingInfoWidget = new QWidget();
  QHBoxLayout *pickingInfoLayout = new QHBoxLayout();
  
  pickingInfoWidget->setStyleSheet ("padding-left: 0px; border:0px;");
  pickingInfoWidget->setLayout(pickingInfoLayout);
  
  vbox->addWidget(pickingInfoWidget);
  // add picking content
  
  fPickInfosScrollArea = new QScrollArea();
  fPickInfosScrollArea->setWidgetResizable(true);

  
  fPickInfosWidget = new QWidget();
  fPickInfosWidget->setStyleSheet ("padding: 0px ");
  
  QVBoxLayout* vLayout = new QVBoxLayout();
  fPickInfosWidget->setLayout (vLayout);
  fPickInfosScrollArea->setWidget(fPickInfosWidget);
  
  QSizePolicy vPolicy = fPickInfosWidget->sizePolicy();
  vPolicy.setVerticalStretch(4);
  vbox->addWidget(fPickInfosScrollArea);
  pickingInfoLayout->setContentsMargins(0,0,0,0);
  vLayout->setContentsMargins(0,0,0,0);
  vbox->setContentsMargins(1,1,1,1);

  groupBox->setLayout(vbox);
  fUIPickInfosWidget->layout()->addWidget(groupBox);
  
  updatePickInfosWidget(fLastPickPoint.x(),fLastPickPoint.y());
}


// set the component to check/unchecked, also go into its child
// and set the same status to all his childs
void G4OpenGLQtViewer::setCheckComponent(QTreeWidgetItem* item,bool check)
{
  if (item) {

    const PVPath& fullPath = fTreeItemModels[item->data(0,Qt::UserRole).toInt()];
    // If a physical volume
    if (fullPath.size() > 0) {
      SetTouchable(fullPath);
      TouchableSetVisibility(fullPath, check);
      fMouseOnSceneTree = true;
    }
  }

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
  QGLWidget* qGLW = dynamic_cast<QGLWidget*> (fGLWidget) ;
  if (! qGLW) {
    return;
  }
  if (isGl2psWriting()) {

    G4OpenGLViewer::DrawText(g4text);

  } else {

    if (!fGLWidget) return;

#ifdef G4MULTITHREADED
    if (G4Threading::G4GetThreadId() != G4Threading::MASTER_ID) return;
#endif

    G4VSceneHandler::MarkerSizeType sizeType;
    G4double size = fSceneHandler.GetMarkerSize(g4text,sizeType);

    QFont font = QFont();
    font.setPointSizeF(size);

    const G4Colour& c = fSceneHandler.GetTextColour(g4text);
    glColor4d(c.GetRed(),c.GetGreen(),c.GetBlue(),c.GetAlpha());
    
    G4Point3D position = g4text.GetPosition();

    const G4String& textString = g4text.GetText();
    const char* textCString = textString.c_str();
  
    glRasterPos3d(position.x(),position.y(),position.z());

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

    qGLW->renderText
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
    createSceneTreeComponent();
  }

  // if no UI
  if (fSceneTreeComponentTreeWidget == NULL) {
    return;
  }

  fSceneTreeComponentTreeWidget->blockSignals(true);

  // Create the "volume" node if not
  //  if (fSceneTreeComponentTreeWidget->topLevelItemCount () == 0) {
  if (!fPVRootNodeCreate) {
    const G4Colour& color = fSceneHandler.GetColour();
    
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
    newItem = new QTreeWidgetItem();
    fSceneTreeComponentTreeWidget->addTopLevelItem(newItem);
  } else {
    newItem = new QTreeWidgetItem(parentTreeNode);
    fSceneTreeComponentTreeWidget->addTopLevelItem(parentTreeNode);
  }


  newItem->setText(0,name);
  newItem->setData(1,Qt::UserRole,copyNb);
  newItem->setText(2,QString::number(POIndex));
  newItem->setData(0, Qt::UserRole, POIndex);
  newItem->setText(3,logicalName);
  newItem->setFlags(newItem->flags()|Qt::ItemIsUserCheckable);
  newItem->setCheckState(0,state);
  newItem->setExpanded(true);
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
  } else {
    // Set a tootip
    newItem->setToolTip (0,QString("double-click to change the color"));
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

  const G4Colour& color = fSceneHandler.GetColour();

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
    color = fSceneHandler.GetColour();
  }
  if (modelShortName == "") {
    return ;
  }
  // try to init it
  if (fSceneTreeComponentTreeWidget == NULL) {
    createSceneTreeComponent();
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
  // Prevent to get out the std::map
  if (fLastSceneTreeWidgetAskForIterator != fLastSceneTreeWidgetAskForIteratorEnd) {
    fLastSceneTreeWidgetAskForIterator++;
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


std::string G4OpenGLQtViewer::parseSceneTreeAndSaveState(){
  std::string commandLine = "";
  for (int b=0;b<fSceneTreeComponentTreeWidget->topLevelItemCount();b++) {
    commandLine += parseSceneTreeElementAndSaveState(fSceneTreeComponentTreeWidget->topLevelItem(b),1)+"\n";
  }
  if (commandLine != "") {
    commandLine = std::string("# Disable auto refresh and quieten vis messages whilst scene and\n") +
    "# trajectories are established:\n" +
    "/vis/viewer/set/autoRefresh false\n" +
    "/vis/verbose errors" +
    commandLine +
    "# Re-establish auto refreshing and verbosity:\n" +
    "/vis/viewer/set/autoRefresh true\n" +
    "/vis/verbose confirmations\n";
  }
  return commandLine;
}


std::string G4OpenGLQtViewer::parseSceneTreeElementAndSaveState(QTreeWidgetItem* item, unsigned int level){
  // parse current item
  std::string str( level, ' ' );
  std::string commandLine = "\n#"+ str + "PV Name: " + item->text(0).toStdString();

  if (item->text(3) != "") {
    commandLine += " LV Name: "+item->text(3).toStdString()+"\n";
    // save check state
    commandLine += "/vis/geometry/set/visibility " + item->text(3).toStdString() + " ! "; // let default value for depth
    if (item->checkState(0) == Qt::Checked) {
      commandLine += "1";
    }
    if (item->checkState(0) == Qt::Unchecked) {
      commandLine += "0";
    }
    commandLine +="\n";
    
    // save color
    const QColor& c = item->data(2,Qt::UserRole).value<QColor>();
    std::stringstream red;
    red << ((double)c.red())/255;
    std::stringstream green;
    green << (double)c.green()/255;
    std::stringstream blue;
    blue << ((double)c.blue())/255;
    std::stringstream alpha;
    alpha << ((double)c.alpha())/255;
    
    commandLine += "/vis/geometry/set/colour " + item->text(3).toStdString() + " ! " + red.str() + " " + green.str() + " " + blue.str() + " " + alpha.str()+"\n";

  } else {
    commandLine += "\n";
  }

  // parse childs
  for (int b=0;b< item->childCount();b++) {
    commandLine += parseSceneTreeElementAndSaveState(item->child(b),level+1);
  }
  
  return commandLine;
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

    fCheckSceneTreeComponentSignalLock = false;
  }
}


void G4OpenGLQtViewer::sceneTreeComponentSelected() {
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

    changeColorAndTransparency(item->data(0,Qt::UserRole).toInt(),
                               G4Colour (((G4double)color.red())/255,
                                         ((G4double)color.green())/255,
                                         ((G4double)color.blue())/255,
                                         ((G4double)color.alpha())/255));

    // set scene tree parameters
    changeQColorForTreeWidgetItem(item,color);
  }
}


void G4OpenGLQtViewer::changeColorAndTransparency(GLuint index, G4Color color) {

  // change vis attributes to set new colour
  G4int iPO = index;
  if (iPO >= 0 && fTreeItemModels.find(iPO) != fTreeItemModels.end()) {
    const PVPath& fullPath = fTreeItemModels[iPO];
    // If a physical volume
    if (fullPath.size()) {
      SetTouchable(fullPath);
      TouchableSetColour(fullPath, color);
      fMouseOnSceneTree = true;
    }
  }
}


G4Colour G4OpenGLQtViewer::getColorForPoIndex(int poIndex) {
  // FIXME 09/2014 : Could be optimize by searching in a tab instead of item->data
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


const std::vector<G4ModelingParameters::VisAttributesModifier>*
G4OpenGLQtViewer::GetPrivateVisAttributesModifiers() const
{
  static std::vector<G4ModelingParameters::VisAttributesModifier>
  privateVisAttributesModifiers;

  privateVisAttributesModifiers.clear();

// I don't think we need this.  (JA Sep 2016).
//  // For each modified touchable...
//  std::map<int,PVPath>::const_iterator i;
//  for (i = fTreeItemModels.begin();
//       i != fTreeItemModels.end();
//       ++i) {
//
//    // How do I know if it's been modified or not?
//
//    int iPO = i->first;
//    const PVPath& fullPath = i->second;
//
//    // If a physical volume
//    if (fullPath.size()) {
//
//      //    const G4bool& visibilityChanged = ???
//      //    const G4bool& visibility = ???
//      //    const G4bool& colourChanged = ???
//      //    const QColor& colour = ???
//      //    G4Colour g4colour(((G4double)colour.red())/255,
//      //                      ((G4double)colour.green())/255,
//      //                      ((G4double)colour.blue())/255,
//      //                      ((G4double)colour.alpha())/255);
//      // Next 4 lines are for testing, to be replaced by the above...
//      G4bool visibilityChanged = true;
//      G4bool visibility = true;
//      G4bool colourChanged = true;
//      G4Colour g4colour(G4Colour::Red());
//
//      // Instantiate a working copy of a G4VisAttributes object...
//      G4VisAttributes workingVisAtts;
//      // ...and use it to create vis attribute modifiers...
//      if (visibilityChanged) {
//        workingVisAtts.SetVisibility(visibility);
//        privateVisAttributesModifiers.push_back
//        (G4ModelingParameters::VisAttributesModifier
//         (workingVisAtts,
//          G4ModelingParameters::VASVisibility,
//          fullPath));
//      }
//      if (colourChanged) {
//        workingVisAtts.SetColour(g4colour);
//        privateVisAttributesModifiers.push_back
//        (G4ModelingParameters::VisAttributesModifier
//         (workingVisAtts,
//          G4ModelingParameters::VASColour,
//          fullPath));
//      }
//    }
//  }

  return &privateVisAttributesModifiers;
}


void G4OpenGLQtViewer::changeSearchSelection()
{
  const QString& searchText = fFilterOutput->text();
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
      fOldLastSceneTreeWidgetAskForIterator = fOldPositivePoIndexSceneTreeWidgetQuickMap.begin();
      fOldLastSceneTreeWidgetAskForIteratorEnd = fOldPositivePoIndexSceneTreeWidgetQuickMap.end();
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
     fLastSceneTreeWidgetAskForIterator = fPositivePoIndexSceneTreeWidgetQuickMap.end();
     fLastSceneTreeWidgetAskForIteratorEnd = fPositivePoIndexSceneTreeWidgetQuickMap.end();
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

  if (fLastSceneTreeWidgetAskForIterator != fLastSceneTreeWidgetAskForIteratorEnd) {
    if (POindex == fLastSceneTreeWidgetAskForIterator->first) {
      if (fLastSceneTreeWidgetAskForIterator->second != NULL) {
        return fLastSceneTreeWidgetAskForIterator->second;
      }
    }
  }
  
  // if not, use the "find" algorithm
  fLastSceneTreeWidgetAskForIterator = fPositivePoIndexSceneTreeWidgetQuickMap.find(POindex);
  fLastSceneTreeWidgetAskForIteratorEnd = fPositivePoIndexSceneTreeWidgetQuickMap.end();
  
  if (fLastSceneTreeWidgetAskForIterator != fPositivePoIndexSceneTreeWidgetQuickMap.end()) {
    return fLastSceneTreeWidgetAskForIterator->second;
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
  // Prevent to get out the std::map
  if (fOldLastSceneTreeWidgetAskForIterator != fOldLastSceneTreeWidgetAskForIteratorEnd) {
    fOldLastSceneTreeWidgetAskForIterator++;
  }
  
  if (fOldLastSceneTreeWidgetAskForIterator != fOldPositivePoIndexSceneTreeWidgetQuickMap.end()) {
    if (POindex == fOldLastSceneTreeWidgetAskForIterator->first) {
      if (fOldLastSceneTreeWidgetAskForIterator->second != NULL) {
        return fOldLastSceneTreeWidgetAskForIterator->second;
      }
    }
  }
  
  // if not, use the "find" algorithm
  fOldLastSceneTreeWidgetAskForIterator = fOldPositivePoIndexSceneTreeWidgetQuickMap.find(POindex);
  fOldLastSceneTreeWidgetAskForIteratorEnd = fOldPositivePoIndexSceneTreeWidgetQuickMap.end();
  
  if (fOldLastSceneTreeWidgetAskForIterator != fOldPositivePoIndexSceneTreeWidgetQuickMap.end()) {
    return fOldLastSceneTreeWidgetAskForIterator->second;
  }
  return NULL;
}



/**
   Should replace actual tree by the one in this class
   and update tree
*/
void G4OpenGLQtViewer::displaySceneTreeComponent() {
  // no UI
  if (fUISceneTreeWidget == NULL) {
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

/**
 Update the scene tree widget
 */
void G4OpenGLQtViewer::updateSceneTreeWidget() {
  // Ensure case where closing a UI tab close the widget
  if (!fSceneTreeWidget) {
    createSceneTreeWidget();
  }
}


  /**
 Update the viewer properties component widget
 Clear it only if the number of command is less than the previous table widget row count
 */
void G4OpenGLQtViewer::updateViewerPropertiesTableWidget() {
  
  if (!isCurrentWidget()) {
    return;
  }
  
  // Ensure case where closing a UI tab close the widget
  if (!fViewerPropertiesTableWidget) {
    createViewerPropertiesWidget();
  }
  int treeWidgetInfosIgnoredCommands = 0;
  G4UImanager* UI = G4UImanager::GetUIpointer();
  G4UIcommandTree * commandTreeTop = UI->GetTree();
  G4UIcommandTree* path = commandTreeTop->FindCommandTree("/vis/viewer/set/");

  if (!path) {
    return;
  }

  // clear old table
  if ((path->GetCommandEntry()-fTreeWidgetInfosIgnoredCommands) != fViewerPropertiesTableWidget->rowCount()) {
    fViewerPropertiesTableWidget->clear();
  }
  
  fViewerPropertiesTableWidget->blockSignals(true);
  // TODO : Could be optimized by comparing current command to old commands. That should not change so much

  fViewerPropertiesTableWidget->setColumnCount (2);
  fViewerPropertiesTableWidget->setRowCount (path->GetCommandEntry()-fTreeWidgetInfosIgnoredCommands);
  fViewerPropertiesTableWidget->setHorizontalHeaderLabels(QStringList() << tr("Property")
                                                        << tr("Value"));
  fViewerPropertiesTableWidget->verticalHeader()->setVisible(false);
  fViewerPropertiesTableWidget->setAlternatingRowColors (true);
  
  // For the moment, we do only command that have a "set" command in UI
  
  for (int a=0;a<path->GetCommandEntry();a++) {
    G4UIcommand* commandTmp = path->GetCommand(a+1);
  
    // get current parameters
    QString params = "";
    
    if(commandTmp->GetCommandName() == "autoRefresh") {
      if (fVP.IsAutoRefresh()) {
        params = "True";
      } else {
        params = "False";
      }
    } else if(commandTmp->GetCommandName() == "auxiliaryEdge") {
      if (fVP.IsAuxEdgeVisible()) {
        params = "True";
      } else {
        params = "False";
      }
    } else if(commandTmp->GetCommandName() == "background") {
      params = QString().number(fVP.GetBackgroundColour().GetRed()) + "  "+
      QString().number(fVP.GetBackgroundColour().GetGreen()) + "  "+
      QString().number(fVP.GetBackgroundColour().GetBlue()) + "  "+
      QString().number(fVP.GetBackgroundColour().GetAlpha());
      
    } else if(commandTmp->GetCommandName() == "culling") {
      params = QString().number(fVP. IsCulling ());
    } else if(commandTmp->GetCommandName() == "cutawayMode") {
      if (fVP.GetCutawayMode() == G4ViewParameters::cutawayUnion) {
        params = "union";
      } else {
        params = "intersection";
      }

    } else if(commandTmp->GetCommandName() == "defaultColour") {
      params = QString().number(fVP.GetDefaultVisAttributes()->GetColor().GetRed()) + "  "+
      QString().number(fVP.GetDefaultVisAttributes()->GetColor().GetGreen()) + "  "+
      QString().number(fVP.GetDefaultVisAttributes()->GetColor().GetBlue()) + "  "+
      QString().number(fVP.GetDefaultVisAttributes()->GetColor().GetAlpha());
      
    } else if(commandTmp->GetCommandName() == "defaultTextColour") {
      params = QString().number(fVP.GetDefaultTextVisAttributes()->GetColor().GetRed()) + "  "+
      QString().number(fVP.GetDefaultTextVisAttributes()->GetColor().GetGreen()) + "  "+
      QString().number(fVP.GetDefaultTextVisAttributes()->GetColor().GetBlue()) + "  "+
      QString().number(fVP.GetDefaultTextVisAttributes()->GetColor().GetAlpha());
      
    } else if(commandTmp->GetCommandName() == "edge") {
      G4ViewParameters::DrawingStyle existingStyle = fVP.GetDrawingStyle();
      params = "False";
      if (existingStyle == G4ViewParameters::hsr) {
        params = "True";
      }

    } else if(commandTmp->GetCommandName() == "explodeFactor") {
      params = QString().number(fVP.GetExplodeFactor()) +  "  " + QString(G4String(G4BestUnit(fVP.GetExplodeFactor(),"Length")).data());

    } else if(commandTmp->GetCommandName() == "globalLineWidthScale") {
      params = QString().number(fVP.GetGlobalLineWidthScale());
      
    } else if(commandTmp->GetCommandName() == "globalMarkerScale") {
      params = QString().number(fVP.GetGlobalMarkerScale());

    } else if(commandTmp->GetCommandName() == "hiddenEdge") {
      G4ViewParameters::DrawingStyle style = fVP.GetDrawingStyle();
      if ((style == G4ViewParameters::hlr) ||
          (style == G4ViewParameters::hlhsr)) {
        params = "True";
      } else {
        params = "False";
      }

    } else if(commandTmp->GetCommandName() == "hiddenMarker") {
      if (fVP.IsMarkerNotHidden()) {
        params = "False";
      } else {
        params = "True";
      }

    } else if(commandTmp->GetCommandName() == "lightsMove") {
      if (fVP.GetLightsMoveWithCamera()) {
        params = "camera";
      } else {
        params = "object";
      }
    } else if(commandTmp->GetCommandName() == "lightsThetaPhi") {
      G4Vector3D direction = fVP.GetLightpointDirection();
      // degree
      params = QString().number(direction.theta()/CLHEP::degree)+ "  "+ QString().number(direction.phi()/CLHEP::degree)+"  deg";
      if (commandTmp->GetParameterEntries() == 3) {
        if (commandTmp->GetParameter(2)->GetDefaultValue() != "deg") {
          params = QString().number(direction.theta())+ "  "+ QString().number(direction.phi())+" "+commandTmp->GetParameter(2)->GetDefaultValue().data();
        }
      }
    } else if(commandTmp->GetCommandName() == "lightsVector") {
     params = QString().number(fVP.GetLightpointDirection().x()) + "  "+
      QString().number(fVP.GetLightpointDirection().y()) + "  "+
      QString().number(fVP.GetLightpointDirection().z());

    } else if(commandTmp->GetCommandName() == "lineSegmentsPerCircle") {
      params = QString().number(fVP.GetNoOfSides());
      
    } else if(commandTmp->GetCommandName() == "picking") {
      if (fVP.IsPicking()) {
        params = "True";
      } else {
        params = "False";
      }

    } else if(commandTmp->GetCommandName() == "projection") {
      if (fVP.GetFieldHalfAngle() == 0.) {
        params = "orthogonal";
      } else {
        params = QString("perspective ") + QString().number(fVP.GetFieldHalfAngle()/CLHEP::degree) + " deg";
      }

    } else if(commandTmp->GetCommandName() == "rotationStyle") {
      if (fVP.GetRotationStyle() == G4ViewParameters::constrainUpDirection) {
        params = "constrainUpDirection";
      } else {
        params = "freeRotation";
      }

    } else if(commandTmp->GetCommandName() == "sectionPlane") {
      if (fVP.IsSection()) {
        params =  QString("on ") +
        G4String(G4BestUnit(fVP.GetSectionPlane().point(),"Length")).data()+
        QString().number(fVP.GetSectionPlane().normal().x())
        + " " + QString().number(fVP.GetSectionPlane().normal().y())
        + " " + QString().number(fVP.GetSectionPlane().normal().z());
      } else {
        params = "off";
      }

    } else if(commandTmp->GetCommandName() == "style") {
      if (fVP.GetDrawingStyle() == G4ViewParameters::wireframe || fVP.GetDrawingStyle() == G4ViewParameters::hlr) {
         params = "wireframe";
      } else {
        params = "surface";
      }

      
    } else if(commandTmp->GetCommandName() == "targetPoint") {
      G4Point3D point = fVP.GetCurrentTargetPoint();
      if (fSceneHandler.GetScene()) {
        G4String b = G4BestUnit(fSceneHandler.GetScene()->GetStandardTargetPoint() + fVP.GetCurrentTargetPoint(),"Length");
        params = b.data();
      }
    } else if(commandTmp->GetCommandName() == "upThetaPhi") {
      G4Vector3D up = fVP.GetUpVector();
      // degree
      params = QString().number(up.theta()/CLHEP::degree)+ "  "+ QString().number(up.phi()/CLHEP::degree)+"  deg";
      if (commandTmp->GetParameterEntries() == 3) {
        if (commandTmp->GetParameter(2)->GetDefaultValue() != "deg") {
          params = QString().number(up.theta())+ "  "+ QString().number(up.phi())+" "+commandTmp->GetParameter(2)->GetDefaultValue().data();
        }
      }
    } else if(commandTmp->GetCommandName() == "upVector") {
      G4Vector3D up = fVP.GetUpVector();
      params = QString().number(up.x())+ "  "+ QString().number(up.y())+"  "+QString().number(up.z())+ "  ";
      
    } else if(commandTmp->GetCommandName() == "viewpointThetaPhi") {
      G4Vector3D direction = fVP.GetViewpointDirection();
      // degree
      params = QString().number(direction.theta()/CLHEP::degree)+ "  "+ QString().number(direction.phi()/CLHEP::degree)+"  deg";
      if (commandTmp->GetParameterEntries() == 3) {
        if (commandTmp->GetParameter(2)->GetDefaultValue() != "deg") {
          params = QString().number(direction.theta())+ "  "+ QString().number(direction.phi())+" "+commandTmp->GetParameter(2)->GetDefaultValue().data();
        }
      }
    } else if(commandTmp->GetCommandName() == "viewpointVector") {
      G4Vector3D direction = fVP.GetViewpointDirection();
      params = QString().number(direction.x())+ "  "+ QString().number(direction.y())+" "+QString().number(direction.z());
    } else {
      // No help
    }
    
    /* DO NOT DISPLAY COMMANDS WITHOUT ANY PARAMETERS SET
    if (params == "") {
      // TODO : display default parameters // should not be editable ?
      
      for( G4int i_thParameter=0; i_thParameter<commandTmp->GetParameterEntries(); i_thParameter++ ) {
        commandParam = commandTmp->GetParameter(i_thParameter);
        
        if (QString(QChar(commandParam->GetParameterType())) == "b") {
          if (commandParam->GetDefaultValue().data()) {
            params += "True";
          } else {
            params += "False";
          }
        } else {
          params += QString((char*)(commandParam->GetDefaultValue()).data());
        }
        if (i_thParameter<commandTmp->GetParameterEntries()-1) {
          params += " ";
        }
      }
    }
    */
    
    if (params != "") {

      QTableWidgetItem *nameItem;
      QTableWidgetItem *paramItem;
      
      // already present ?
      QList<QTableWidgetItem *>	list = fViewerPropertiesTableWidget->findItems (commandTmp->GetCommandName().data(),Qt::MatchExactly);
      if (list.size() == 1) {
        nameItem = list.first();
        paramItem = fViewerPropertiesTableWidget->item(nameItem->row(),1);
        
      } else {
        nameItem = new QTableWidgetItem();
        paramItem = new QTableWidgetItem();
        fViewerPropertiesTableWidget->setItem(a-treeWidgetInfosIgnoredCommands, 0, nameItem);
        fViewerPropertiesTableWidget->setItem(a-treeWidgetInfosIgnoredCommands, 1, paramItem);
        
        // Set Guidance
        QString guidance;
        G4int n_guidanceEntry = commandTmp->GetGuidanceEntries();
        for( G4int i_thGuidance=0; i_thGuidance < n_guidanceEntry; i_thGuidance++ ) {
          guidance += QString((char*)(commandTmp->GetGuidanceLine(i_thGuidance)).data()) + "\n";
        }
        
        nameItem->setToolTip(guidance);
        paramItem->setToolTip(GetCommandParameterList(commandTmp));
        
        fViewerPropertiesTableWidget->setRowHeight(a-treeWidgetInfosIgnoredCommands,15);
      }
      
      // set current name and parameters
      nameItem->setText(commandTmp->GetCommandName().data());
      paramItem->setText(params);
      
      nameItem->setFlags(Qt::NoItemFlags);
      nameItem->setForeground(QBrush());
      
    } else {
      treeWidgetInfosIgnoredCommands++;
    }
  }
  // remove empty content row
  for (int i=0; i<treeWidgetInfosIgnoredCommands; i++) {
    fViewerPropertiesTableWidget->removeRow (fViewerPropertiesTableWidget->rowCount() - 1);
  }

  // The resize should done only at creation
  if (!fViewerPropertiesTableWidgetIsInit) {
    fViewerPropertiesTableWidgetIsInit = true;
    
    fViewerPropertiesTableWidget->resizeColumnsToContents();
    
    int x = fViewerPropertiesTableWidget->horizontalHeader()->length();
    int y = fViewerPropertiesTableWidget->verticalHeader()->length()+ fViewerPropertiesTableWidget->horizontalHeader()->sizeHint().height() + 2;
    
    //  fViewerPropertiesTableWidget->verticalHeader()->setResizeMode(QHeaderView::ResizeToContents);
    //  fViewerPropertiesTableWidget->horizontalHeader()->setResizeMode(QHeaderView::ResizeToContents);
    
    // resize to fit content
    QDialog* dial = static_cast<QDialog*> (fUIViewerPropertiesWidget->parent());
    if (dial) {
      dial->resize(x+56,y+46); // more or less (margins) ...
    }
  }
  fViewerPropertiesTableWidget->blockSignals(false);

  fTreeWidgetInfosIgnoredCommands = treeWidgetInfosIgnoredCommands;
}


/**
 Update the pick infos component widget
 */
void G4OpenGLQtViewer::updatePickInfosWidget(int aX, int aY) {
  fLastPickPoint = QPoint(aX,aY);
  
  if (!isCurrentWidget()) {
    return;
  }
  // Ensure case where closing a UI tab close the widget
  if (!fPickInfosWidget) {
    createPickInfosWidget();
  }

  const std::vector < G4OpenGLViewerPickMap* > & pickMap = GetPickDetails(aX,aY);
  
  // remove all previous widgets
  if (fPickInfosWidget) {
    QLayoutItem * wItem;
    if (fPickInfosWidget->layout()->count()) {
      while ((wItem = fPickInfosWidget->layout()->takeAt(0)) != 0) {
        delete wItem->widget();
        delete wItem;
      }
    }
  } else {
    // Ensure case where closing a UI tab close the widget
    if (!fPickInfosWidget) {
      createPickInfosWidget();
    }
  }

  // Create a new signalMapper
  delete fSignalMapperPicking;
  fSignalMapperPicking = new QSignalMapper(this);
  
  // parse all pick results
  for (unsigned int a=0; a< pickMap.size(); a++) {
    // Add a box inside the pick viewer box
    std::ostringstream label;
    if (pickMap[a]->getAttributes().size() > 0) {
      std::string txt = pickMap[a]->getAttributes()[0].data();

      // Look for Volumes
      std::size_t pos = txt.find("Physical Volume Path (PVPath):") + 31;
      if (pos != 31-1) {
        std::size_t pos2 = txt.find("\n",pos);
        label << "Volume: " + txt.substr(pos,pos2-pos);
      } else  {
        // Look for tracks
        if (pickMap[a]->getAttributes().size() > 1) {
          std::string txt1 = pickMap[a]->getAttributes()[1].data();
          
          pos = txt.find("Run ID (RunID):") + 16;
          if (pos != 16-1) {
            std::size_t pos2 = txt.find("\n",pos);
            std::size_t posEID = txt.find("Event ID (EventID):",pos2) + 20;
            std::size_t posEID2 = txt.find("\n",posEID);
            std::size_t posPN = txt1.find("Particle Name (PN):") + 20;
            std::size_t posPN2 = txt1.find("\n",posPN);
            std::size_t posCh = txt1.find("Charge (Ch):",posPN2) + 13;
            std::size_t posCh2 = txt1.find("\n",posCh);
            label << "RunID:" << txt.substr(pos,pos2-pos)
            << " EventID:" << txt.substr(posEID,posEID2-posEID)
            << " PN:" << txt1.substr(posPN,posPN2-posPN)
            << " Ch:" << txt1.substr(posCh,posCh2-posCh);
          }
        } else {
          label << "Hit number:" << a << ", PickName: " << pickMap[a]->getPickName();
        }
          
      }
    } else {
      label << "Hit number:" << a << ", PickName: " << pickMap[a]->getPickName();
    }
    QPushButton* pickCoutButton = new QPushButton(label.str().c_str());
    pickCoutButton->setStyleSheet ("text-align: left; padding: 1px; border: 0px;");
    pickCoutButton->setIcon(*fTreeIconClosed);
    fPickInfosWidget->layout()->addWidget(pickCoutButton);
    
    QStringList newStr;
    
    // Add to stringList
    newStr = QStringList(QString(pickMap[a]->print().data()).trimmed());
    
    QTextEdit* ed = new QTextEdit();
    ed->setFontFamily("Courier");
    ed->setFontPointSize(12);
    ed->setReadOnly(true);
    fPickInfosWidget->layout()->addWidget(ed);
    ed->setVisible((false));
    ed->append(newStr.join(""));
    
    connect(pickCoutButton, SIGNAL(clicked()), fSignalMapperPicking, SLOT(map()));
    fSignalMapperPicking->setMapping(pickCoutButton,fPickInfosWidget->layout()->count()-1);
  }
  
  connect(fSignalMapperPicking, SIGNAL(mapped(int)),this, SLOT(toggleSceneTreeComponentPickingCout(int)));
  
  // add a label to push everything up!
  QLabel * pushUp = new QLabel("");
  QSizePolicy vPolicy = QSizePolicy(QSizePolicy::Minimum,QSizePolicy::Minimum);
  vPolicy.setVerticalStretch(10);
  pushUp->setSizePolicy(vPolicy);
  fPickInfosWidget->layout()->addWidget(pushUp);
  
/*  // highlight the first one :
  
  // first un-highlight the last selected
  changeColorAndTransparency(fLastHighlightName,fLastHighlightColor);
  
  if (pickMap.size() > 0 ) {
    // get the new one
    fLastHighlightName = pickMap[0]->getPickName();
    fLastHighlightColor = getColorForPoIndex(fLastHighlightName);
    // set the new one
    changeColorAndTransparency(fLastHighlightName,G4Color(1,1,1,1));
    
    updateQWidget();
  }
 */
  QDialog* dial = static_cast<QDialog*> (fUIPickInfosWidget->parent());
  if (dial) {
    // change name
    std::ostringstream oss;
    if (pickMap.size() == 0) {
      oss << "No object";
    } else if (pickMap.size() == 1) {
      oss << "1 object";
    } else {
      oss << QString::number(pickMap.size()).toStdString() << " objects";
    }
    oss << " selected - " << GetName();
    dial->setWindowTitle(oss.str().c_str());
  }
  // set picking cout visible
  fPickInfosScrollArea->setVisible(true);
}


void G4OpenGLQtViewer::toggleSceneTreeComponentPickingCout(int pickItem) {

  QWidget* w;
  // close other items, it could take too much space
  
  for (int a=0; a<fPickInfosWidget->layout()->count(); a++) {
    w = fPickInfosWidget->layout()->itemAt(a)->widget();
    QTextEdit* ed = dynamic_cast<QTextEdit*>(w);
    QPushButton* button;
    if (ed) {
      if (a == pickItem) {
        w->setVisible(!w->isVisible());
      } else {
        w->setVisible(false);
      }
      if (a >= 1) {
        button = dynamic_cast<QPushButton*>(fPickInfosWidget->layout()->itemAt(a-1)->widget());
        if (button) {
          if (button->isVisible()) {
            button->setIcon(*fTreeIconOpen);
          } else {
            button->setIcon(*fTreeIconClosed);
          }
        }
      }
    }
  }
}


void G4OpenGLQtViewer::currentTabActivated(int currentTab) {
  if (fUiQt->GetViewerTabWidget()->tabText(currentTab) == GetName())  {
    createViewerPropertiesWidget();
    createPickInfosWidget();
    createSceneTreeWidget();
  }
}


void G4OpenGLQtViewer::tableWidgetViewerSetItemChanged(QTableWidgetItem * item) {
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI != NULL)  {
    QTableWidgetItem* previous = fViewerPropertiesTableWidget->item(fViewerPropertiesTableWidget->row(item),0);
    if (previous) {
      fViewerPropertiesTableWidget->blockSignals(true);
      UI->ApplyCommand((std::string("/vis/viewer/set/")
                        + previous->text().toStdString()
                        + " "
                        + item->text().toStdString()).c_str());
      fViewerPropertiesTableWidget->blockSignals(false);
    }
  }
}

bool G4OpenGLQtViewer::isCurrentWidget(){
  // Prevent from repainting a hidden tab (the current tab name has to be the one of th GL viewer)
  if ( GetName() != fUiQt->GetViewerTabWidget()->tabText(fUiQt->GetViewerTabWidget()->currentIndex()).toStdString().c_str()) {
    return false;
  }
  return true;
}

/**   Build the parameter list parameters in a QString<br>
 Reimplement partialy the G4UIparameter.cc
 @param aCommand : command to list parameters
 @see G4UIparameter::List()
 @see G4UIcommand::List()
 @return the command list parameters, or "" if nothing
 */
QString G4OpenGLQtViewer::GetCommandParameterList (
                                  const G4UIcommand *aCommand
                                  )
{
  G4int n_parameterEntry = aCommand->GetParameterEntries();
  QString txt;
  
  if( n_parameterEntry > 0 ) {
    G4UIparameter *param;
    
    // Re-implementation of G4UIparameter.cc
    
    for( G4int i_thParameter=0; i_thParameter<n_parameterEntry; i_thParameter++ ) {
      param = aCommand->GetParameter(i_thParameter);
      txt += "\nParameter : " + QString((char*)(param->GetParameterName()).data()) + "\n";
      if( ! param->GetParameterGuidance().isNull() )
        txt += QString((char*)(param->GetParameterGuidance()).data())+ "\n" ;
      txt += " Parameter type  : " + QString(QChar(param->GetParameterType())) + "\n";
      if(param->IsOmittable()){
        txt += " Omittable       : True\n";
      } else {
        txt += " Omittable       : False\n";
      }
      if( param->GetCurrentAsDefault() ) {
        txt += " Default value   : taken from the current value\n";
      } else if( ! param->GetDefaultValue().isNull() ) {
        txt += " Default value   : " + QString((char*)(param->GetDefaultValue()).data())+ "\n";
      }
      if( ! param->GetParameterRange().isNull() ) {
        txt += " Parameter range : " + QString((char*)(param->GetParameterRange()).data())+ "\n";
      }
      if( ! param->GetParameterCandidates().isNull() ) {
        txt += " Candidates      : " + QString((char*)(param->GetParameterCandidates()).data())+ "\n";
      }
    }
  }
  return txt;
}

#ifdef G4MULTITHREADED

void G4OpenGLQtViewer::DoneWithMasterThread()
{
  // Called by Main Thread !
  
  // Useful to avoid two vis thread at the same time
  G4MUTEXLOCK(&mWaitForVisSubThreadQtOpenGLContextInitialized);
}

void G4OpenGLQtViewer::SwitchToVisSubThread()
{
  // Called by VisSub Thread !

  QGLWidget* qGLW = dynamic_cast<QGLWidget*> (fGLWidget) ;
  if (! qGLW) {
    return;
  }
  
  // Set the current QThread to its static variable
  SetQGLContextVisSubThread(QThread::currentThread());
  
  // - Wait for the vis thread to set its QThread
  G4CONDITIONBROADCAST(&c1_VisSubThreadQtOpenGLContextInitialized);
  
  // Unlock the vis thread if it is Qt Viewer
  G4CONDITIONWAIT(&c2_VisSubThreadQtOpenGLContextMoved, &mWaitForVisSubThreadQtOpenGLContextMoved);

  // make context current
  qGLW->makeCurrent();
}

void G4OpenGLQtViewer::DoneWithVisSubThread()
{
  // Called by vis sub thread
  QGLWidget* qGLW = dynamic_cast<QGLWidget*> (fGLWidget) ;
  if (! qGLW) {
    return;
  }

  // finish with this vis sub thread context
  qGLW->doneCurrent();

#if QT_VERSION > 0x050000
  // and move it back to the main thread
  qGLW->context()->moveToThread(fQGLContextMainThread);
#endif
}

void G4OpenGLQtViewer::SwitchToMasterThread()
{
  // Called by VisSub Thread !

  QGLWidget* qGLW = dynamic_cast<QGLWidget*> (fGLWidget) ;
  if (! qGLW) {
    return;
  }

  // Useful to avoid two vis thread at the same time
  G4MUTEXUNLOCK(&mWaitForVisSubThreadQtOpenGLContextInitialized);
  qGLW->makeCurrent();
}


void G4OpenGLQtViewer::MovingToVisSubThread(){
  // Called by Main Thread !

  QGLWidget* qGLW = dynamic_cast<QGLWidget*> (fGLWidget) ;
  if (! qGLW) {
    return;
  }

  // - Wait for the vis sub thread to set its QThread
  G4CONDITIONWAIT( &c1_VisSubThreadQtOpenGLContextInitialized, &mWaitForVisSubThreadQtOpenGLContextInitialized);


  // Set current QThread for the way back
  SetQGLContextMainThread(QThread::currentThread());
  
  // finish with this main thread context
  qGLW->doneCurrent();
#if QT_VERSION > 0x050000
  qGLW->context()->moveToThread(fQGLContextVisSubThread);
#endif
  
  G4CONDITIONBROADCAST(&c2_VisSubThreadQtOpenGLContextMoved);
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
#endif
