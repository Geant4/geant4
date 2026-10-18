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
//
//
// Alex Vergara Gil  17th October 2026
// Ported to Vulkan from the Geant4 OpenGL visualization driver.

#ifndef G4VULKANQTMOVIEDIALOG_HH
#define G4VULKANQTMOVIEDIALOG_HH

#include <qdialog.h>

class QPushButton;
class QLabel;
class QLineEdit;
class G4VulkanQtViewer;

class QGroupBox;

/** The G4VulkanQtMovieDialog class provide a Dialog displaying differents options
    for each file format
*/
class G4VulkanQtMovieDialog : public QDialog
{
    Q_OBJECT

  public:

    /** Construct a G4VulkanQtMovieDialog
        @param parent : parent widget
    */
    G4VulkanQtMovieDialog(G4VulkanQtViewer*, QWidget*);

    /** Destroys G4VulkanQtMovieDialog */
    ~G4VulkanQtMovieDialog();
    void setRecordingStatus(QString);
    void setRecordingInfos(QString);

  private:

    QLineEdit* fEncoderPath;
    QLineEdit* fTempFolderPath;
    QLineEdit* fSaveFileName;
    G4VulkanQtViewer* fParentViewer;
    QLabel* fEncoderStatus;
    QLabel* fTempFolderStatus;
    QLabel* fSaveFileStatus;
    QLabel* fRecordingStatus;
    QLabel* fRecordingInfos;
    QPushButton* fButtonStopFinishClose;
    QPushButton* fButtonSave;
    QPushButton* fButtonStartPause;

  public Q_SLOTS:
    void stopFinishClose();
    void save();
    bool checkEncoderSwParameters();
    bool checkSaveFileNameParameters();
    bool checkTempFolderParameters();

  private Q_SLOTS:
    void selectEncoderPathAction();
    void selectTempPathAction();
    void selectSaveFileNameAction();
    void resetRecording();
    void enabledApplyButton();
};

#endif
