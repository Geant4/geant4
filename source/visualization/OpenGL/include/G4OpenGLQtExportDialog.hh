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
// $Id: G4OpenGLQtExportDialog.hh,v 1.7 2009-02-04 16:48:40 lgarnier Exp $
// GEANT4 tag $Name: 
//
// 

#ifdef G4VIS_BUILD_OPENGLQT_DRIVER

#ifndef G4OPENGLQTEXPORTDIALOG_HH
#define G4OPENGLQTEXPORTDIALOG_HH

#include <qdialog.h>

class QButtonGroup;
class QPushButton;
class QRadioButton;
class QCheckBox;
class QSlider;
class QComboBox;
class QLabel;
class QLineEdit;

#if QT_VERSION >= 0x040000
class QGroupBox;
#endif

/** The G4OpenGLQtExportDialog class provide a Dialog displaying differents options
    for each file format
*/
class G4OpenGLQtExportDialog : public QDialog
{
  Q_OBJECT

 public:
  /** Construct a G4OpenGLQtExportDialog
      @param parent : parent widget
      @param format : format of save file in lower case
      @param height : height of the original file
      @param width : width of the original file
  */
  G4OpenGLQtExportDialog(QWidget* parent, QString format, int height =0, int width=0);

  /** Destroys G4OpenGLQtExportDialog */
  ~G4OpenGLQtExportDialog();

  /** @return the value of the slider if format has a slider widget, instead return -1 */
  int getSliderValue();

  /** return the new width for file if format has a width widget, instead return 
      the original value */
  int getWidth();

  /** return the new height for file if format has a height widget, instead return
      the original value  */
  int getHeight();

  /** return the value of the transparency if format has a transparency, instead return -1 */
  int getTransparency();

  /** return the numbers of colors in the picture if format has a BW/Color widget, instead return -1 */
  int getNbColor();

  /** return if vector EPS is checked, if button does'nt exist, return 0 */
  bool getVectorEPS();

	    public slots:

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
  void textWidthChanged(const QString &); 

  /** Called by changing value in width lineEdit. If ratio is keep, will also change the height
   */
  void textHeightChanged(const QString &); 

 private:
  QString f_name, f_type, f_dir;
  QLabel* qualityLabel;
  bool expAll;
  QPushButton* buttonOk;
  QPushButton* buttonCancel;

#if QT_VERSION >= 0x040000
  QGroupBox * sizeGroupBox;
#endif

  QCheckBox* transparencyEPS,*boxTransparency,*vectorEPSCheckBox;
  QCheckBox* ratioCheckBox;
  QSlider * qualitySlider;
  QLabel *formatLabel;
  QRadioButton* colorButton,*BWButton;
  QRadioButton* original,* modify;
  QLineEdit* height,*width;
  QWidget* heightWidget,* widthWidget;
  int originalWidth;
  int originalHeight;
  bool isChangingSize;
};

#endif

#endif
