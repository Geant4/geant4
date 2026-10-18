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

#ifndef G4VULKANQTEXPORTDIALOG_HH
#define G4VULKANQTEXPORTDIALOG_HH

#include <qdialog.h>

class QButtonGroup;
class QPushButton;
class QRadioButton;
class QCheckBox;
class QSlider;
class QComboBox;
class QLabel;
class QLineEdit;

class QGroupBox;

/** The G4VulkanQtExportDialog class provide a Dialog displaying differents options
    for each file format
*/
class G4VulkanQtExportDialog : public QDialog
{
    Q_OBJECT

  public:

    /** Construct a G4VulkanQtExportDialog
        @param parentw : parent widget
        @param format  : format of save file in lower case
        @param height  : height of the original file
        @param width   : width of the original file
    */
    G4VulkanQtExportDialog(QWidget* parentw, QString format, int height = 0, int width = 0);

    /** Destroys G4VulkanQtExportDialog */
    ~G4VulkanQtExportDialog();

    /** @return the value of the slider if format has a slider widget, instead return -1 */
    int getSliderValue();

    /** return the new width for file if format has a width widget, instead return
        the original value */
    int getWidth();

    /** return the new height for file if format has a height widget, instead return
        the original value  */
    int getHeight();

    /** return if vector EPS is checked, if button does'nt exist, return 0 */
    bool getVectorEPS();

  public Q_SLOTS:

    /** Called by a clic on modify/original size button.This will
   invert buttons and hide/unhide size
     */
    void changeSizeBox();

    /** Called by a clic on vectorEPS check box.If vectorEPS checkBox is checked,
        it will enable change size buttons. Else it will disable them.
    */
    void changeVectorEPS();

    /** Called by changing value in height lineEdit. If ratio is keep, will also change the width
     */
    void textWidthChanged(const QString&);

    /** Called by changing value in width lineEdit. If ratio is keep, will also change the height
     */
    void textHeightChanged(const QString&);

  private:

    QString f_name, f_type, f_dir;
    //  QLabel* qualityLabel;
    //  bool expAll;
    QPushButton* buttonOk;
    QPushButton* buttonCancel;

    QGroupBox* sizeGroupBox;

    QCheckBox /** transparencyEPS,*boxTransparency,*/* vectorEPSCheckBox;
    QCheckBox* ratioCheckBox;
    QSlider* qualitySlider;
    //  QLabel *formatLabel;
    QRadioButton *colorButton, *BWButton;
    QRadioButton *original, *modify;
    QLineEdit *height, *width;
    QWidget *heightWidget, *widthWidget;
    int originalWidth;
    int originalHeight;
    bool isChangingSize;
};

#endif
