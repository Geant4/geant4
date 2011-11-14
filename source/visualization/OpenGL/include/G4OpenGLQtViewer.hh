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
// $Id: G4OpenGLQtViewer.hh,v 1.25 2010-10-08 10:07:31 lgarnier Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// G4OpenGLQtViewer : Class to provide WindowsNT specific
//                       functionality for OpenGL in GEANT4

#ifdef G4VIS_BUILD_OPENGLQT_DRIVER

#ifndef G4OPENGLQTVIEWER_HH
#define G4OPENGLQTVIEWER_HH

#include "globals.hh"

#include "G4OpenGLViewer.hh"

#include <qobject.h>
#include <qpoint.h>

class G4OpenGLSceneHandler;
class G4UImanager;

class QGLWidget;
class QDialog;
class QContextMenuEvent;
class QMenu;
class QImage;
class QAction;
class QMouseEvent;
class QKeyEvent;
class QWheelEvent;
class QProcess;
class QTime;
class QVBoxLayout;
class QTreeWidgetItem;
class QTreeWidget;
class G4OpenGLSceneHandler;
class G4OpenGLQtMovieDialog;

class G4OpenGLQtViewer: public QObject, virtual public G4OpenGLViewer {

  Q_OBJECT

public:
  G4OpenGLQtViewer (G4OpenGLSceneHandler& scene);
  virtual ~G4OpenGLQtViewer ();
  virtual void updateQWidget()=0;
  QString setEncoderPath(QString path);
  QString getEncoderPath();
  QString setTempFolderPath(QString path);
  QString getTempFolderPath();
  QString setSaveFileName(QString path);
  QString getSaveFileName();
  bool isRecording();
  bool isStopped();
  bool isPaused();
  bool isEncoding();
  bool isWaiting();
  bool isFailed();
  void setWaiting();
  bool isBadEncoder();
  bool isBadOutput();
  bool isBadTmp();
  bool isSuccess();
  void setBadTmp();
  void setBadOutput();
  void setBadEncoder();
  bool isReadyToEncode();
  void resetRecording();
  void encodeVideo();
  void stopVideo();
  void saveVideo();
  bool generateMpegEncoderParameters();
  void displayRecordingStatus();
  void DrawText(const char * ,double x,double y,double z, double size);
  void ResetView ();
  void addTreeElement(const G4String model,std::vector < std::pair<std::string,std::pair <unsigned int, unsigned int> > >);
  bool isTouchableVisible(unsigned int POindex);

public:
  void G4MousePressEvent(QMouseEvent *event);
  void G4wheelEvent (QWheelEvent * event); 
  void G4keyPressEvent (QKeyEvent * event); 
  void G4MouseDoubleClickEvent();
  void G4MouseReleaseEvent();
  void G4MouseMoveEvent(QMouseEvent *event);

protected:
  void CreateGLQtContext ();
  void CreateFontLists ();
  virtual void CreateMainWindow (QGLWidget*,QString);
  void G4manageContextMenuEvent(QContextMenuEvent *e);
  void rotateQtScene(float, float);
  void rotateQtSceneToggle(float, float);
  void moveScene(float, float, float,bool);
  void FinishView();
  void updateKeyModifierState(Qt::KeyboardModifiers);
  void displayViewComponentTree();


protected:
  QGLWidget* fWindow;
  QWidget* fGLWindow;
  bool hasPendingEvents();
  void savePPMToTemp();
  int fRecordFrameNumber;

  bool fHasToRepaint;
  bool fReadyToPaint;
  bool fIsRepainting;

private:
  enum mouseActions {STYLE1,STYLE2,STYLE3,STYLE4}; 
  enum RECORDING_STEP {WAIT,START,PAUSE,CONTINUE,STOP,READY_TO_ENCODE,ENCODING,FAILED,SUCCESS,BAD_ENCODER,BAD_OUTPUT,BAD_TMP,SAVE}; 

  void createPopupMenu();
  void createRadioAction(QAction *,QAction *, const std::string&,unsigned int a=1);
  void rescaleImage(int, int);
  bool printPDF(const std::string,int,QImage);  
  void showMovieParametersDialog();
  void initMovieParameters();
  QString createTempFolder();
  QString removeTempFolder();
  void setRecordingStatus(RECORDING_STEP);
  void setRecordingInfos(QString);
  QString getProcessErrorMsg();
  QWidget* getParentWidget();
  bool parseAndInsertInTree(QTreeWidgetItem *,std::vector < std::pair<std::string,std::pair <unsigned int, unsigned int> > > treeVect,QString parentRoot);
  void setCheckComponent(QTreeWidgetItem* item,bool check);
  void initViewComponent();
  bool parseAndCheckVisibility(QTreeWidgetItem * treeNode,unsigned int POindex);

  QMenu *fContextMenu;

  mouseActions fMouseAction; // 1: rotate 2:move 3:pick 4:shortcuts 
  QPoint fLastPos1;
  QPoint fLastPos2;
  QPoint fLastPos3;

  /** delta of depth move. This delta is put in % of the scene view */
  G4double fDeltaDepth;
  /** delta of zoom move. This delta is put in % of the scene view */
  G4double fDeltaZoom;
  /** To ensure key event are keep one by one */
  bool fHoldKeyEvent;
  /** To ensure move event are keep one by one */
  bool fHoldMoveEvent;
  /** To ensure rotate event are keep one by one */
  bool fHoldRotateEvent;
  bool fAutoMove;
  QString fEncoderPath;
  QString fTempFolderPath;
  QString fMovieTempFolderPath;
  QString fSaveFileName;
  QString fParameterFileName;
  QAction *fRotateAction;
  QAction *fMoveAction;
  QAction *fPickAction;
  QAction *fFullScreenOn;
  QAction *fFullScreenOff;
  QAction *fDrawingWireframe;
  QAction *fDrawingLineRemoval;
  QAction *fDrawingSurfaceRemoval;
  QAction *fDrawingLineSurfaceRemoval;
  G4OpenGLQtMovieDialog* fMovieParametersDialog;
  RECORDING_STEP fRecordingStep;
  QProcess *fProcess;
  QTime *fLastEventTime;
  int fSpinningDelay;
  int fNbMaxFramesPerSec;
  float fNbMaxAnglePerSec;
  int fLaunchSpinDelay;
  QWidget* fUIViewComponentsTBWidget;
  QVBoxLayout* fLayoutViewComponentsTBWidget;
  bool fNoKeyPress;
  bool fAltKeyPress;
  bool fControlKeyPress;
  bool fShiftKeyPress;
  bool fBatchMode;
  bool fCheckViewComponentLock;
  QTreeWidget * fViewerComponentTreeWidget;
  int fNbRotation ;
  int fTimeRotation;

public Q_SLOTS :
  void startPauseVideo();

private Q_SLOTS :
  void actionMouseRotate();
  void actionMouseMove();
  void actionMousePick();
  void actionDrawingWireframe();
  void actionDrawingLineRemoval();
  void actionDrawingSurfaceRemoval();
  void actionDrawingLineSurfaceRemoval();
  void actionSaveImage();
  void actionChangeBackgroundColor();
  void actionChangeTextColor();
  void actionChangeDefaultColor();
  void actionMovieParameters();

  void showShortcuts();
  void toggleDrawingAction(int);
  void toggleMouseAction(mouseActions);
  void toggleRepresentation(bool);
  void toggleProjection(bool);
  void toggleTransparency(bool);
  void toggleAntialiasing(bool);
  void toggleHaloing(bool);
  void toggleAux(bool);
  void toggleFullScreen(bool);
  void processEncodeFinished();
  void processLookForFinished();
  void processEncodeStdout();
  void viewComponentItemChanged(QTreeWidgetItem* item, int id);
  // Only use for Qt>4.0
  //  void dialogClosed();
};

#endif

#endif
