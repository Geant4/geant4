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
// $Id: G4OpenGLWtViewer.hh 75567 2013-11-04 11:35:11Z gcosmo $
//
// 
// G4OpenGLWtViewer : Class to provide WindowsNT specific
//                       functionality for OpenGL in GEANT4

#ifdef G4VIS_BUILD_OPENGLWT_DRIVER

#ifndef G4OPENGLWTVIEWER_HH
#define G4OPENGLWTVIEWER_HH

#include "globals.hh"

#include "G4OpenGLViewer.hh"

#include <Wt/WObject>
#include <Wt/WPoint>
#include <Wt/WTime>
#include <Wt/WContainerWidget>
#include <Wt/WMatrix4x4>

class G4OpenGLSceneHandler;
class G4UImanager;

class WDialog;
#ifdef _A_FINIR_FIXME
class WContextMenuEvent;
#endif
class WMenu;
class WImage;
#ifdef _A_FINIR_FIXME
class WWheelEvent;
#endif
class WProcess;
class G4UIWt;

class G4OpenGLSceneHandler;
class G4OpenGLWtMovieDialog;

class G4OpenGLWtViewer:  virtual public G4OpenGLViewer {

public:
  G4OpenGLWtViewer (G4OpenGLSceneHandler& scene);
  virtual ~G4OpenGLWtViewer ();
  virtual void updateWWidget()=0;

  Wt::WMatrix4x4 mMatrix;

#ifdef _A_FINIR_FIXME
  Wt::WString setEncoderPath(Wt::WString path);
  Wt::WString getEncoderPath();
  Wt::WString setTempFolderPath(Wt::WString path);
  Wt::WString getTempFolderPath();
  Wt::WString setSaveFileName(Wt::WString path);
  Wt::WString getSaveFileName();
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
#endif
  void drawText(const char * ,int x,int y,int z, int size);
protected:
  void CreateGLWtContext ();
  virtual void CreateMainWindow (Wt::WGLWidget*,Wt::WString);
#ifdef _A_FINIR_FIXME
  void G4manageContextMenuEvent(Wt::WContextMenuEvent *e);
#endif
  void G4MousePressEvent(Wt::WMouseEvent *event);
  void G4MouseReleaseEvent();
  void G4MouseDoubleClickEvent();
  void G4MouseMoveEvent(Wt::WMouseEvent *event);
//  void G4wheelEvent (Wt::WWheelEvent * event);
  void G4keyPressEvent (Wt::WKeyEvent * event); 
  void rotateWtScene(float, float);
  void rotateWtSceneToggle(float, float);
  void moveScene(float, float, float,bool);
#ifdef _A_FINIR_FIXME
  void updateKeyModifierState(Wt::KeyboardModifiers);
#endif
  
  inline Wt::WGLWidget* getGLWindow() {
    return fWindow;
  }
  

protected:
  Wt::WGLWidget* fWindow;
  Wt::WWidget* fGLWindow;
  bool hasPendingEvents();
  void resizeGL(int width, int height);
  
#ifdef _A_FINIR_FIXME
  void savePPMToTemp();
#endif
  int fRecordFrameNumber;

  bool fHasToRepaint;
  bool fReadyToPaint;
  bool fIsRepainting;

private:
  enum mouseActions {STYLE1,STYLE2,STYLE3,STYLE4}; 
  enum RECORDING_STEP {WAIT,START,PAUSE,CONTINUE,STOP,READY_TO_ENCODE,ENCODING,FAILED,SUCCESS,BAD_ENCODER,BAD_OUTPUT,BAD_TMP,SAVE}; 

#ifdef _A_FINIR_FIXME
  void createPopupMenu();
  void createRadioAction(Wt::WAction *,Wt::WAction *, const std::string&,unsigned int a=1);
#endif
  void rescaleImage(int, int);
#ifdef _A_FINIR_FIXME
  bool printPDF(const std::string,int,WImage);  
  void showMovieParametersDialog();
  void initMovieParameters();
  Wt::WString createTempFolder();
  Wt::WString removeTempFolder();
  void setRecordingStatus(RECORDING_STEP);
  void setRecordingInfos(Wt::WString);
  Wt::WString getProcessErrorMsg();

  WMenu *fContextMenu;
#endif

  mouseActions fMouseAction; // 1: rotate 2:move 3:pick 4:shortcuts 
  Wt::WPoint fLastPos1;
  Wt::WPoint fLastPos2;
  Wt::WPoint fLastPos3;
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
  Wt::WString fEncoderPath;
  Wt::WString fTempFolderPath;
  Wt::WString fMovieTempFolderPath;
  Wt::WString fSaveFileName;
  Wt::WString fParameterFileName;
#ifdef _A_FINIR_FIXME
  WAction *fRotateAction;
  WAction *fMoveAction;
  WAction *fPickAction;
  WAction *fFullScreenOn;
  WAction *fFullScreenOff;
  WAction *fDrawingWireframe;
  WAction *fDrawingLineRemoval;
  WAction *fDrawingSurfaceRemoval;
  WAction *fDrawingLineSurfaceRemoval;
#endif
  G4OpenGLWtMovieDialog* fMovieParametersDialog;
  RECORDING_STEP fRecordingStep;
  WProcess *fProcess;
  Wt::WTime *fLastEventTime;
  int fSpinningDelay;
  int fNbMaxFramesPerSec;
  float fNbMaxAnglePerSec;
  int fLaunchSpinDelay;
  Wt::WTabWidget* fUISceneTreeComponentsTBWidget;

  G4double fXRot;
  G4double fYRot;
  bool fNoKeyPress;
  bool fAltKeyPress;
  bool fControlKeyPress;
  bool fShiftKeyPress;
  bool fBatchMode;
  G4UIWt* fUiWt;

#ifdef _A_FINIR_FIXME
 void rotateTheta(int);
 void rotatePhi(int);
 void moveX(int);
 void moveY(int);
 void moveZ(int);
#endif

public :
  void startPauseVideo();

private :
#ifdef _A_FINIR_FIXME
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
#endif

  void showShortcuts();
#ifdef _A_FINIR_FIXME
  void toggleDrawingAction(int);
  void toggleMouseAction(mouseActions);
  void toggleRepresentation(bool);
  void toggleProjection(bool);
  void toggleTransparency(bool);
  void toggleAntialiasing(bool);
  void toggleHaloing(bool);
  void toggleAux(bool);
  void toggleFullScreen(bool);
#endif
  void processEncodeFinished();
  void processLookForFinished();
  void processEncodeStdout();
  // Only use for Wt>4.0
  //  void dialogClosed();


};

#endif

#endif
