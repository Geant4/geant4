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
// $Id: G4OpenGLQtViewer.hh,v 1.9 2008-03-11 16:05:56 lgarnier Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// G4OpenGLQtViewer : Class to provide WindowsNT specific
//                       functionality for OpenGL in GEANT4

#ifdef G4VIS_BUILD_OPENGLQT_DRIVER

#ifndef G4OPENGLQTVIEWER_HH
#define G4OPENGLQTVIEWER_HH

#include "globals.hh"

#include "G4VViewer.hh"
#include "G4OpenGLSceneHandler.hh"

#include <qobject.h>
#include <qpoint.h>

class QGLWidget;
class QDialog;
class QContextMenuEvent;
#if QT_VERSION < 0x040000
class QPopupMenu;
#else
class QMenu;
#endif
class QImage;
class QAction;
class QMouseEvent;
class QKeyEvent;
class QProcess;

class G4OpenGLSceneHandler;
class G4OpenGLQtMovieDialog;

class G4OpenGLQtViewer: public QObject, virtual public G4OpenGLViewer {

  Q_OBJECT

public:
  G4OpenGLQtViewer (G4OpenGLSceneHandler& scene);
  virtual ~G4OpenGLQtViewer ();
  void SetView ();
  void ShowView ();
  virtual void updateQWidget()=0;
  void setupViewport(int, int);
  QString setEncoderPath(QString path);
  QString getEncoderPath();
  QString setTempFolderPath(QString path);
  QString getTempFolderPath();
  QString setSaveFileName(QString path);
  QString getSaveFileName();
  bool isRecording();
  bool isStopped();
  bool isReadyToEncode();
  void resetRecording();
  void encodeVideo();
  bool generateMpegEncoderParameters();
  void displayRecordingStatus();

protected:
  void CreateGLQtContext ();
  virtual void CreateMainWindow (QGLWidget*);
  void manageContextMenuEvent(QContextMenuEvent *e);
#if QT_VERSION < 0x040000
  void G4MousePressEvent(QPoint, Qt::ButtonState);
#else
  void G4MousePressEvent(QPoint, Qt::MouseButtons);
#endif
  void G4MouseDoubleClickEvent(QPoint p);
#if QT_VERSION < 0x040000
  void G4MouseMoveEvent(int, int, Qt::ButtonState,bool mAutoMove = false);
#else
  void G4MouseMoveEvent(int, int, Qt::MouseButtons,bool mAutoMove = false);
#endif
  void G4keyPressEvent (QKeyEvent * event); 
  void rotateScene(G4double, G4double,bool mAutoRotate=false);
  void moveScene(G4double, G4double, G4double,bool,bool mAutoMove=false);


protected:
  G4int WinSize_x;
  G4int WinSize_y;
  QGLWidget* fWindow;
  QDialog* GLWindow;
  bool hasPendingEvents();
  void savePPMToTemp();
  int fRecordFrameNumber;

private:
  enum mouseActions {STYLE1,STYLE2,STYLE3,STYLE4}; 
  enum RECORDING_STEP {WAIT,START,PAUSE,CONTINUE,STOP,READY_TO_ENCODE,ENCODING,FAILED,SUCCESS}; 

  void createPopupMenu();
  void createRadioAction(QAction *,QAction *, const std::string&,unsigned int a=1);
  void rescaleImage(int, int);
  bool generateEPS(QString,int,QImage);  
  bool generateVectorEPS (QString,int,int,QImage);
  bool generatePS_PDF(QString,int,QImage);  
  void showMovieParametersDialog();
  void initMovieParameters();
  QString createTempFolder();
  QString removeTempFolder();
  void startPauseVideo();
  void stopVideo();
  void setRecordingStatus(RECORDING_STEP);
  void setRecordingInfos(QString);
  QString getProcessErrorMsg();


#if QT_VERSION < 0x040000
  QPopupMenu *fContextMenu;
#else
  QMenu *fContextMenu;
#endif

  mouseActions fMouseAction; // 1: rotate 2:move 3:pick 4:shortcuts 
  QPoint fLastPos;
  /** delta X of move event */
  G4double fDeltaPosX;
  /** delta Y of move event */
  G4double fDeltaPosY;
  /** delta of scene rotation. This delta is put in degree */
  G4double fDeltaRotation;
  /** delta of scene translation. This delta is put in % of the scene view */
  G4double fDeltaSceneTranslation;
  /** delta of depth move. This delta is put in % of the scene view */
  G4double fDeltaDepth;
  /** delta of zoom move. This delta is put in % of the scene view */
  G4double fDeltaZoom;
  /** delta of auto move/rotation. This delta is put in % of the move/rotation param */
  G4double fDeltaMove;
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

private slots :
  void actionMouseRotate();
  void actionMouseMove();
  void actionMousePick();
  void actionDrawingWireframe();
  void actionDrawingLineRemoval();
  void actionDrawingSurfaceRemoval();
  void actionDrawingLineSurfaceRemoval();
  void actionSaveImage();
  void actionMovieParameters();

  void showShortcuts();
  void toggleDrawingAction(int);
  void toggleMouseAction(mouseActions);
  void toggleRepresentation(bool);
  void toggleProjection(bool);
  void toggleBackground(bool);
  void toggleTransparency(bool);
  void toggleAntialiasing(bool);
  void toggleHaloing(bool);
  void toggleAux(bool);
  void toggleFullScreen(bool);
  void processEncodeFinished();
  void processLookForFinished();
  void processEncodeStdout();
  // Only use for Qt>4.0
  //  void dialogClosed();
};

#endif

#endif
