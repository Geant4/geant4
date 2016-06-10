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
// $Id: G4OpenGLQtMovieDialog.hh 66373 2012-12-18 09:41:34Z gcosmo $
// GEANT4 tag $Name: 
//
// 

#ifdef G4VIS_BUILD_OPENGLQT_DRIVER

#ifndef G4OPENGLQTMOVIEDIALOG_HH
#define G4OPENGLQTMOVIEDIALOG_HH

#include <qdialog.h>

class QPushButton;
class QLabel;
class QLineEdit;
class G4OpenGLQtViewer;

class QGroupBox;

/** The G4OpenGLQtMovieDialog class provide a Dialog displaying differents options
    for each file format
*/
class G4OpenGLQtMovieDialog : public QDialog
{
  Q_OBJECT

 public:
  /** Construct a G4OpenGLQtMovieDialog
      @param parent : parent widget
  */
  G4OpenGLQtMovieDialog(G4OpenGLQtViewer*,QWidget*);

  /** Destroys G4OpenGLQtMovieDialog */
  ~G4OpenGLQtMovieDialog();
  void setRecordingStatus(QString);
  void setRecordingInfos(QString);

private :
  QLineEdit* fEncoderPath;
  QLineEdit* fTempFolderPath;
  QLineEdit* fSaveFileName;
  G4OpenGLQtViewer *fParentViewer;
  QLabel *fEncoderStatus;
  QLabel *fTempFolderStatus;
  QLabel *fSaveFileStatus;
  QLabel *fRecordingStatus;
  QLabel *fRecordingInfos;
  QPushButton *fButtonStopFinishClose;
  QPushButton *fButtonSave;
  QPushButton *fButtonStartPause;

public Q_SLOTS :
  void stopFinishClose();
  void save();
  bool checkEncoderSwParameters();
  bool checkSaveFileNameParameters();
  bool checkTempFolderParameters();

private Q_SLOTS :
  void selectEncoderPathAction();
  void selectTempPathAction();
  void selectSaveFileNameAction();
  void resetRecording();
  void enabledApplyButton();
};

#endif

#endif
