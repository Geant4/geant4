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
// $Id: G4OpenGLQtMovieDialog.cc,v 1.8 2008-10-24 13:49:19 lgarnier Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifdef G4VIS_BUILD_OPENGLQT_DRIVER

#include "G4OpenGLQtMovieDialog.hh"
#include "G4OpenGLQtViewer.hh"

#include <qpushbutton.h>
#include <qpalette.h>
#include <qlabel.h>
#include <qgroupbox.h>
#include <qlayout.h>
#include <qlineedit.h>
#include <qfiledialog.h>
#include <qprocess.h>


// +---------------------------------------+
// +        Path for encoder               +
// +  _______                              +
// + | select| ____________________        +
// +  -------                              +
// +        Temp path                      +
// +  _______                              +
// + | select| ____________________        +
// +  -------                              +
// +                                       +
// + max number of frames  ________        +
// + ....                                  +
// +                                       +
// +     Label : X frames Saves/Encoding   +
// +         Cancel        Encode          +
// +---------------------------------------+

G4OpenGLQtMovieDialog::G4OpenGLQtMovieDialog(
 G4OpenGLQtViewer* parentViewer,
 QWidget* parentWidget
)
  : QDialog( parentWidget ),
    fParentViewer(parentViewer)
{
#if QT_VERSION > 0x030200
  setModal(false);
#endif
#if QT_VERSION < 0x040000
  setCaption( tr( " Movie parameters" ));
#else
  setWindowTitle( tr( " Movie parameters" ));
#endif


  // global layout
  QVBoxLayout* globalVLayout = new QVBoxLayout(this);
  globalVLayout->setMargin(10);
  globalVLayout->setSpacing(10);
  
  // Encoder group box
  QGroupBox *encoderGroupBox = new QGroupBox(tr("Encoder path"),this);		
  QVBoxLayout *encoderVGroupBoxLayout = new QVBoxLayout(encoderGroupBox);

  // Encoder Path 
  QWidget *encoderHBox = new QWidget(encoderGroupBox);
  QHBoxLayout *encoderHBoxLayout = new QHBoxLayout(encoderHBox);
  fEncoderPath = new QLineEdit("",encoderHBox);

  QPushButton *encoderButton = new QPushButton(tr("..."),encoderHBox);
  encoderButton->setMaximumWidth (30);

  fEncoderStatus = new QLabel(encoderGroupBox);

#if QT_VERSION > 0x040000
  fEncoderStatus->setWordWrap(true);
  encoderVGroupBoxLayout->setMargin(15);
#else
  fEncoderStatus->setAlignment ( Qt::AlignAuto |Qt::WordBreak );
  encoderGroupBox->setInsideMargin(15);
#endif

  fEncoderStatus->setText("");

#if QT_VERSION < 0x040000
  encoderHBoxLayout->add(fEncoderPath);
  encoderHBoxLayout->add(encoderButton);
  encoderVGroupBoxLayout->add(encoderHBox);
  encoderVGroupBoxLayout->add(fEncoderStatus);

  globalVLayout->add(encoderGroupBox);
#else
  encoderHBoxLayout->addWidget(fEncoderPath);
  encoderHBoxLayout->addWidget(encoderButton);
  encoderVGroupBoxLayout->addWidget(encoderHBox);
  encoderVGroupBoxLayout->addWidget(fEncoderStatus);

  encoderGroupBox->setLayout(encoderVGroupBoxLayout);
  globalVLayout->addWidget(encoderGroupBox);
#endif

  connect( encoderButton, SIGNAL( clicked( ) ), this, SLOT(selectEncoderPathAction() ) );


  // temp folder group box
  QGroupBox *tempFolderGroupBox = new QGroupBox(tr("Temporary folder path"),this);
  QVBoxLayout *tempFolderVGroupBoxLayout = new QVBoxLayout(tempFolderGroupBox);

  // temp folder Path 
  QWidget *tempFolderHBox = new QWidget(tempFolderGroupBox);
  QHBoxLayout *tempFolderHBoxLayout = new QHBoxLayout(tempFolderHBox);

  fTempFolderPath = new QLineEdit("",tempFolderHBox);

  QPushButton *tempButton = new QPushButton(tr("..."),tempFolderHBox);
  tempButton->setMaximumWidth (30);

  fTempFolderStatus = new QLabel(tempFolderGroupBox);
#if QT_VERSION > 0x040000
  fTempFolderStatus->setWordWrap(true);
  tempFolderVGroupBoxLayout->setMargin(15);
#else
  fTempFolderStatus->setAlignment ( Qt::AlignAuto |Qt::WordBreak );
  tempFolderGroupBox->setInsideMargin(15);
#endif
  fTempFolderStatus->setText("");

#if QT_VERSION < 0x040000
  tempFolderHBoxLayout->add(fTempFolderPath);
  tempFolderHBoxLayout->add(tempButton);
  tempFolderVGroupBoxLayout->add(tempFolderHBox);
  tempFolderVGroupBoxLayout->add(fTempFolderStatus);

  globalVLayout->add(tempFolderGroupBox);
#else
  tempFolderHBoxLayout->addWidget(fTempFolderPath);
  tempFolderHBoxLayout->addWidget(tempButton);
  tempFolderVGroupBoxLayout->addWidget(tempFolderHBox);
  tempFolderVGroupBoxLayout->addWidget(fTempFolderStatus);

  tempFolderGroupBox->setLayout(tempFolderVGroupBoxLayout);
  globalVLayout->addWidget(tempFolderGroupBox);
#endif

  connect( tempButton, SIGNAL( clicked( ) ), this, SLOT(selectTempPathAction() ) );




  // save file group box
  QGroupBox *saveFileGroupBox = new QGroupBox(tr("Save as"),this);
  QVBoxLayout *saveFileVGroupBoxLayout = new QVBoxLayout(saveFileGroupBox);

  // save file 
  QWidget *saveFileHBox = new QWidget(saveFileGroupBox);
  QHBoxLayout *saveFileHBoxLayout = new QHBoxLayout(saveFileHBox);

  fSaveFileName = new QLineEdit("G4Movie.mpeg",saveFileHBox);

  QPushButton *saveButton = new QPushButton(tr("..."),saveFileHBox);
  saveButton->setMaximumWidth (30);

  fSaveFileStatus = new QLabel(saveFileGroupBox);
#if QT_VERSION > 0x040000
  fSaveFileStatus->setWordWrap(true);
  saveFileVGroupBoxLayout->setMargin(15);
#else
  fSaveFileStatus->setAlignment ( Qt::AlignAuto |Qt::WordBreak );
  saveFileGroupBox->setInsideMargin(15);
#endif
  fSaveFileStatus->setText("");

#if QT_VERSION < 0x040000
  saveFileHBoxLayout->add(fSaveFileName);
  saveFileHBoxLayout->add(saveButton);
  saveFileVGroupBoxLayout->add(saveFileHBox);
  saveFileVGroupBoxLayout->add(fSaveFileStatus);

  globalVLayout->add(saveFileGroupBox);
#else
  saveFileHBoxLayout->addWidget(fSaveFileName);
  saveFileHBoxLayout->addWidget(saveButton);
  saveFileVGroupBoxLayout->addWidget(saveFileHBox);
  saveFileVGroupBoxLayout->addWidget(fSaveFileStatus);

  saveFileGroupBox->setLayout(saveFileVGroupBoxLayout);
  globalVLayout->addWidget(saveFileGroupBox);
#endif

  connect( saveButton, SIGNAL( clicked( ) ), this, SLOT(selectSaveFileNameAction() ) );



  // label

  QLabel *infoLabel = new QLabel("  Press SPACE to Start/Pause video recording \n  Press RETURN to Stop video recording",this);

  // global status
  QGroupBox *statusGroupBox = new QGroupBox(tr("Status"),this);
  QVBoxLayout *statusVGroupBoxLayout = new QVBoxLayout(statusGroupBox);

  fRecordingStatus = new QLabel(statusGroupBox);
#if QT_VERSION > 0x040000
  statusVGroupBoxLayout->setMargin(15);
  fRecordingStatus->setWordWrap(true);
#else
  statusGroupBox->setInsideMargin(15);
  fRecordingStatus->setAlignment ( Qt::AlignAuto |Qt::WordBreak );
#endif
  QPalette palette( fRecordingStatus->palette() );
#if QT_VERSION > 0x040000
  palette.setColor( QPalette::Text, Qt::green);
#else
  palette.setColor( QColorGroup::Text, Qt::green);
#endif
  fRecordingStatus->setPalette(palette);

  fRecordingInfos = new QLabel(statusGroupBox);
#if QT_VERSION > 0x040000
  fRecordingInfos->setWordWrap(true);
#else
  fRecordingInfos->setAlignment ( Qt::AlignAuto |Qt::WordBreak );
#endif
  setRecordingInfos("");

#if QT_VERSION < 0x040000
  statusVGroupBoxLayout->add(fRecordingStatus);
  statusVGroupBoxLayout->add(fRecordingInfos);

  globalVLayout->add(infoLabel);
  globalVLayout->add(statusGroupBox);
#else
  statusVGroupBoxLayout->addWidget(fRecordingStatus);
  statusVGroupBoxLayout->addWidget(fRecordingInfos);

  statusGroupBox->setLayout(statusVGroupBoxLayout);
  globalVLayout->addWidget(infoLabel);
  globalVLayout->addWidget(statusGroupBox);
#endif

  // buttons
  QWidget *buttonBox = new QWidget(this);

  QHBoxLayout *buttonBoxLayout = new QHBoxLayout(buttonBox);

  QPushButton *buttonReset = new QPushButton( tr( "&Reset" ),buttonBox );
  buttonReset->setAutoDefault( TRUE );
  buttonBoxLayout->addWidget(buttonReset);

  fButtonStartPause = new QPushButton( tr( "  &Start " ),buttonBox );
  fButtonStartPause->setEnabled(true);
  fButtonStartPause->setAutoDefault( TRUE );
  buttonBoxLayout->addWidget(fButtonStartPause);

  fButtonStopFinishClose = new QPushButton( tr( "&Stop" ),buttonBox );
  fButtonStopFinishClose->setEnabled(false);
  fButtonStopFinishClose->setAutoDefault( TRUE );
  buttonBoxLayout->addWidget(fButtonStopFinishClose);

  fButtonSave = new QPushButton( tr( "&Save" ),buttonBox );
  fButtonSave->setEnabled(false);
  fButtonSave->setAutoDefault( TRUE );
  buttonBoxLayout->addWidget(fButtonSave);

  QPushButton *buttonCancel = new QPushButton( tr( "&Cancel" ),buttonBox );
  buttonCancel->setAutoDefault( TRUE );
  buttonBoxLayout->addWidget(buttonCancel);

#if QT_VERSION >= 0x040000
  buttonBox->setLayout(buttonBoxLayout);
#endif
  globalVLayout->addWidget(buttonBox);



#if QT_VERSION >= 0x040000
  setLayout(globalVLayout);
#endif

  // signals and slots connections
  connect( fButtonStartPause, SIGNAL( clicked() ), fParentViewer, SLOT(    startPauseVideo() ) );
  connect( buttonReset, SIGNAL( clicked() ), this, SLOT( resetRecording() ) );
  connect( buttonCancel, SIGNAL( clicked() ), this, SLOT( reject() ) );
  connect( fButtonStopFinishClose, SIGNAL( clicked() ), this, SLOT( stopFinishClose() ) );
  connect( fButtonSave, SIGNAL( clicked() ), this, SLOT( save() ) );

  // fill
  setRecordingStatus("");
  fEncoderPath->setText(fParentViewer->getEncoderPath());
  fTempFolderPath->setText(fParentViewer->getTempFolderPath());

  // connect line edit signals
  connect (fEncoderPath,SIGNAL(textChanged ( const QString&)),this,SLOT(checkEncoderSwParameters()));
  connect (fTempFolderPath,SIGNAL(textChanged ( const QString&)),this,SLOT(checkTempFolderParameters()));
  connect (fSaveFileName,SIGNAL(textChanged ( const QString&)),this,SLOT(checkSaveFileNameParameters()));

#if QT_VERSION >= 0x040000
  connect (fEncoderPath,SIGNAL(editingFinished ()),this,SLOT(checkEncoderSwParameters()));
  connect (fTempFolderPath,SIGNAL(editingFinished ()),this,SLOT(checkTempFolderParameters()));
  connect (fSaveFileName,SIGNAL(editingFinished ()),this,SLOT(checkSaveFileNameParameters()));
#else
  connect (fEncoderPath,SIGNAL(lostFocus ()),this,SLOT(checkEncoderSwParameters()));
  connect (fTempFolderPath,SIGNAL(lostFocus ()),this,SLOT(checkTempFolderParameters()));
  connect (fSaveFileName,SIGNAL(lostFocus ()),this,SLOT(checkSaveFileNameParameters()));
#endif

}



G4OpenGLQtMovieDialog::~G4OpenGLQtMovieDialog()
{
}

void G4OpenGLQtMovieDialog::selectEncoderPathAction()
{
#if QT_VERSION < 0x040000
  QString nomFich =  QFileDialog::getOpenFileName ( ".",
                                                    NULL,
                                                    this,
                                                    "Select your encoder",
                                                    tr("Select your encoder ..."));
#else
  QString nomFich =  QFileDialog::getOpenFileName ( this,
                                                    "Select your encoder",
                                                    tr("Select your encoder ...")); 

#endif

  if (nomFich == "") {
    return;
  }
  fEncoderPath->setText(nomFich);
  checkEncoderSwParameters();
 }


void G4OpenGLQtMovieDialog::selectTempPathAction()
{
#if QT_VERSION < 0x040000
  QString nomFich =  QFileDialog::getExistingDirectory (NULL,
                                                    this,
                                                    "Select temporary folder",
                                                    tr("Select temporary folder ..."));
#else
  QString nomFich =  QFileDialog::getExistingDirectory ( this,
                                                    "Select temporary folder",
                                                    tr("Select temporary folder ...")); 

#endif
  if (nomFich == "") {
    return;
  }
  fTempFolderPath->setText(nomFich);
  checkTempFolderParameters();
 }


void G4OpenGLQtMovieDialog::selectSaveFileNameAction()
{
#if QT_VERSION < 0x040000
  QString nomFich =  QFileDialog::getSaveFileName ( ".",
                                                    NULL,
                                                    this,
                                                    "Select saved file",
                                                    tr("Select saved file ..."));
#else
  QString nomFich =  QFileDialog::getSaveFileName ( this,
                                                    "Select saved file",
                                                    tr("Select saved file ...")); 

#endif
  if (nomFich == "") {
    return;
  }
  fSaveFileName->setText(nomFich);
  checkSaveFileNameParameters();
 }


void G4OpenGLQtMovieDialog::stopFinishClose() {
  fParentViewer->stopVideo();
}

void G4OpenGLQtMovieDialog::save() {
  if (((fParentViewer->isPaused()) || fParentViewer->isRecording() || fParentViewer->isStopped())) {
    fParentViewer->saveVideo();
  }
}

	/**
 * If one of parameter is incorrect, put it in red and don't valid it
 * If valid, save it
 */
bool G4OpenGLQtMovieDialog::checkEncoderSwParameters() {

  bool status = true;
  QPalette palette( fEncoderPath->palette() );

  QString temp = fParentViewer->setEncoderPath(fEncoderPath->text());
  setRecordingInfos("");
  fEncoderStatus->setText(temp);
  if (temp != "") {
#if QT_VERSION > 0x040000
    palette.setColor( QPalette::Base, Qt::red);
#else
    palette.setColor( QColorGroup::Base, Qt::red);
#endif
    if (fParentViewer->isReadyToEncode()) {
      setRecordingInfos("No valid encode defined, screen capture had been saved in the temp folder in ppm format.\nPlease define a encoder and clic on Apply button");
	}
    status = false;
  } else {
#if QT_VERSION > 0x040000
    palette.setColor( QPalette::Base, Qt::white);
#else
    palette.setColor( QColorGroup::Base, Qt::white);
#endif
    fEncoderPath->setText(fParentViewer->getEncoderPath());
  }
  fEncoderPath->setPalette(palette);
  return status;
}


/**
 * If one of parameter is incorrect, put it in red and don't valid it
 * If valid, save it
 */
bool G4OpenGLQtMovieDialog::checkTempFolderParameters() {

  bool status = true;
  QPalette palette( fTempFolderPath->palette() );

  QString temp = fParentViewer->setTempFolderPath(fTempFolderPath->text());
  fTempFolderStatus->setText(temp);
  if (temp != "") {
#if QT_VERSION > 0x040000
    palette.setColor( QPalette::Base, Qt::red);
#else
    palette.setColor( QColorGroup::Base, Qt::red);
#endif
    status = false;
  } else {
#if QT_VERSION > 0x040000
    palette.setColor( QPalette::Base, Qt::white);
#else
    palette.setColor( QColorGroup::Base, Qt::white);
#endif
    fTempFolderPath->setText(fParentViewer->getTempFolderPath());
  }
  fTempFolderPath->setPalette(palette);
  return status;
}


/**
 * If one of parameter is incorrect, put it in red and don't valid it
 * If valid, save it
 */
bool G4OpenGLQtMovieDialog::checkSaveFileNameParameters() {

  bool status = true;
  QPalette palette( fSaveFileName->palette() );

  QString temp = fParentViewer->setSaveFileName(fSaveFileName->text());
  fSaveFileStatus->setText(temp);
  if (temp != "") { 
#if QT_VERSION > 0x040000
    palette.setColor( QPalette::Base, Qt::red);
#else
    palette.setColor( QColorGroup::Base, Qt::red);
#endif
    status = false;
  } else {
#if QT_VERSION > 0x040000
    palette.setColor( QPalette::Base, Qt::white);
#else
    palette.setColor( QColorGroup::Base, Qt::white);
#endif
    fSaveFileName->setText(fParentViewer->getSaveFileName());
  }
  fSaveFileName->setPalette(palette);
  return status;
}


void G4OpenGLQtMovieDialog::resetRecording() {
  fParentViewer->resetRecording();
}


void G4OpenGLQtMovieDialog::setRecordingStatus(QString txt) {
  fRecordingStatus->setText(txt);
  if (fParentViewer->isWaiting()) {
    fButtonStartPause->setText("  &Start ");
    fButtonStartPause->setEnabled(true);
    fButtonStopFinishClose->setEnabled(false);
    fButtonSave->setEnabled(false);

  } else if (fParentViewer->isPaused()) {

    fButtonStartPause->setText("  &Continue ");
    fButtonStartPause->setEnabled(true);
    fButtonStopFinishClose->setEnabled(true);
    fButtonSave->setEnabled(false);

  } else if (fParentViewer->isRecording()) {

    fButtonStartPause->setText("  &Pause ");
    fButtonStartPause->setEnabled(true);
    fButtonStopFinishClose->setEnabled(true);
    fButtonSave->setEnabled(false);

  } else if (fParentViewer->isBadOutput()) {

    fButtonStartPause->setText("  &Start ");
    fButtonStartPause->setEnabled(true);
    fButtonStopFinishClose->setEnabled(false);
    fButtonSave->setEnabled(false);

  } else if (fParentViewer->isBadTmp()) {

    fButtonStartPause->setText("  &Start ");
    fButtonStartPause->setEnabled(false);
    fButtonStopFinishClose->setEnabled(false);
    fButtonSave->setEnabled(false);

  } else if (fParentViewer->isBadEncoder()) {

    fButtonStartPause->setText("  &Start ");
    fButtonStartPause->setEnabled(true);
    fButtonStopFinishClose->setEnabled(false);
    fButtonSave->setEnabled(false);

  } else if (fParentViewer->isSuccess()) {

    fButtonStartPause->setText("  &Start ");
    fButtonStartPause->setEnabled(false);
    fButtonStopFinishClose->setEnabled(false);
    fButtonSave->setEnabled(false);

  } else if (fParentViewer->isFailed()) {

    fButtonStartPause->setText("  &Start ");
    fButtonStartPause->setEnabled(false);
    fButtonStopFinishClose->setEnabled(false);
    fButtonSave->setEnabled(false);

  } else if (fParentViewer->isStopped()) {

    fButtonStartPause->setText("  &Start ");
    fButtonStartPause->setEnabled(false);
    fButtonStopFinishClose->setEnabled(false);
    fButtonSave->setEnabled(true);
  }
}


void G4OpenGLQtMovieDialog::setRecordingInfos(QString txt) {
  fRecordingInfos->setText(txt);
}


void G4OpenGLQtMovieDialog::enabledApplyButton() {
  fButtonStartPause->setEnabled(true);
}

#endif
