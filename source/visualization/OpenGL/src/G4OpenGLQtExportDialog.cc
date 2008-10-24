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
// $Id: G4OpenGLQtExportDialog.cc,v 1.9 2008-10-24 14:17:10 lgarnier Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifdef G4VIS_BUILD_OPENGLQT_DRIVER

#include "G4OpenGLQtExportDialog.hh"

#include <qvariant.h>
#include <qpushbutton.h>
#include <qcheckbox.h>
#include <qlabel.h>
#include <qcombobox.h>
#include <qslider.h>
#include <qlayout.h>
#include <qgroupbox.h>
#include <qradiobutton.h>
#include <qimage.h>
#include <qlineedit.h>
#include <qbuttongroup.h>

G4OpenGLQtExportDialog::G4OpenGLQtExportDialog(
 QWidget* parent
,QString format
 ,int aHeight
 ,int aWidth
)
  : QDialog( parent ),
    isChangingSize(false)
{
#if QT_VERSION < 0x040000
  setCaption( tr( " Export options" ));
#else
  setWindowTitle( tr( " Export options" ));
#endif
  originalWidth = aWidth;
  originalHeight = aHeight;

  // Initializations
  qualitySlider = NULL;
  width = NULL;
  height = NULL;
  colorButton = NULL;
  BWButton = NULL;

  // global layout
  QVBoxLayout* globalVLayout = new QVBoxLayout(this);
  globalVLayout->setMargin(10);
  globalVLayout->setSpacing(10);
  


  // FIXME : L. Garnier 4/12/07
  // Not implented. Should deal with alpha channel

//   if((format == "tif") ||
//      (format == "tiff") ||
//      (format == "jpg") ||
//      (format == "jpeg") ||
//      (format == "png") ||
//      (format == "xpm")) {

//     QGroupBox *transparencyGroupBox = new QGroupBox(tr("Transparency"),this);
//     QVBoxLayout *transparencyGroupBoxLayout = new QVBoxLayout(transparencyGroupBox);

//     boxTransparency = new QCheckBox("Save transparency",transparencyGroupBox);
//     boxTransparency->setChecked( false );

//     transparencyGroupBoxLayout->addWidget(boxTransparency);    
// #if QT_VERSION >= 0x040000
//     transparencyGroupBox->setLayout(transparencyGroupBoxLayout);
// #endif
//     globalVLayout->addWidget(transparencyGroupBox);

//   }

  // FIXME : L. Garnier 4/12/07
  // This is not working for PS and PDF images, it does nothing.
  // Image is staying in color mode
  //  if ((format == "ps") || (format == "pdf") || (format == "eps")) {


  // size box

  QWidget * sizeWidget = new QWidget(this); // widget containing group button
  QVBoxLayout * sizeWidgetLayout = new QVBoxLayout(sizeWidget);
  sizeWidgetLayout->setMargin (10); 

  // original and modify radiobuttons
#if QT_VERSION < 0x040000
  QButtonGroup * sizeButtonGroupBox = new QButtonGroup ( 2,Qt::Vertical, tr("Size"),this);
  sizeButtonGroupBox->setInsideMargin (15); 

  original = new QRadioButton("Original",sizeButtonGroupBox);
  modify = new QRadioButton("Modify",sizeButtonGroupBox);

  sizeButtonGroupBox->insert(original);
  sizeButtonGroupBox->insert(modify);
  sizeButtonGroupBox->setExclusive(true);
  sizeWidgetLayout->add(sizeButtonGroupBox);

  connect( sizeButtonGroupBox, SIGNAL( clicked(int) ), this, SLOT( changeSizeBox()) );
#else
  
  sizeGroupBox = new QGroupBox(tr("Size"));
  QVBoxLayout *sizeGroupBoxLayout = new QVBoxLayout(sizeGroupBox);
  QButtonGroup * sizeButtonGroupBox = new QButtonGroup();
  sizeGroupBoxLayout->setMargin (15); 

  original = new QRadioButton("Original");
  modify = new QRadioButton("Modify");

  sizeButtonGroupBox->addButton(original);
  sizeButtonGroupBox->addButton(modify);
  sizeButtonGroupBox->setExclusive(true);

  sizeGroupBoxLayout->addWidget(original);    
  sizeGroupBoxLayout->addWidget(modify);    

  sizeGroupBox->setLayout(sizeGroupBoxLayout);
  sizeWidgetLayout->addWidget(sizeGroupBox);
  
  connect( sizeButtonGroupBox, SIGNAL( buttonClicked(QAbstractButton*) ), this, SLOT( changeSizeBox()) );
#endif
  original->setChecked( true );


  // height
  heightWidget = new QWidget(sizeWidget);

  QHBoxLayout *heightLineLayout = new QHBoxLayout(heightWidget);

  QString tmp;
 
  heightLineLayout->addWidget(new QLabel("Height",heightWidget));
  height = new QLineEdit(tmp.setNum(originalHeight),heightWidget);
  height->setMaxLength(5);
#if QT_VERSION < 0x040000
  heightLineLayout->add(height);
#else
  heightLineLayout->addWidget(height);
#endif

#if QT_VERSION >= 0x040000
  heightWidget->setLayout(heightLineLayout);
#endif

#if QT_VERSION < 0x040000
  sizeWidgetLayout->add(heightWidget);
#else
  sizeWidgetLayout->addWidget(heightWidget);
#endif
  connect( height, SIGNAL( textChanged ( const QString& ) ), this, SLOT( textHeightChanged(const QString &) ) );


  // width
  widthWidget = new QWidget(sizeWidget);

  QHBoxLayout *widthLineLayout = new QHBoxLayout(widthWidget);

#if QT_VERSION < 0x040000
  widthLineLayout->add(new QLabel("Width ",widthWidget));
#else
  widthLineLayout->addWidget(new QLabel("Width ",widthWidget));
#endif
  width = new QLineEdit(tmp.setNum(originalWidth),widthWidget);
  width->setMaxLength(5);
#if QT_VERSION < 0x040000
  widthLineLayout->add(width);
#else
  widthLineLayout->addWidget(width);
#endif
#if QT_VERSION >= 0x040000
  widthWidget->setLayout(widthLineLayout);
#endif
#if QT_VERSION < 0x040000
  sizeWidgetLayout->add(widthWidget);
#else
  sizeWidgetLayout->addWidget(widthWidget);
#endif
  connect( width, SIGNAL( textChanged ( const QString& ) ), this, SLOT( textWidthChanged(const QString &) ) );



  // ratio check box

  ratioCheckBox = new QCheckBox( "Keep ratio",sizeWidget);
  ratioCheckBox->setChecked( true );

#if QT_VERSION < 0x040000
  sizeWidgetLayout->add(ratioCheckBox);
#else
  sizeWidgetLayout->addWidget(ratioCheckBox);
#endif

#if QT_VERSION < 0x040000
  ratioCheckBox->setEnabled ( false );
  heightWidget->setEnabled ( false );
  widthWidget->setEnabled ( false );
#else
  ratioCheckBox->hide();
  heightWidget->hide();
  widthWidget->hide();
#endif

#if QT_VERSION >= 0x040000
  sizeWidget->setLayout(sizeWidgetLayout);
#endif
  globalVLayout->addWidget(sizeWidget);

 if (format == "eps") {

   QGroupBox *EPSWidgetGroupBox = new QGroupBox(tr("EPS options"),this); // widget containing group button


#if QT_VERSION < 0x040000

    EPSWidgetGroupBox->setInsideMargin (15); 

    //    QButtonGroup * EPSColorButtonGroupBox = new QButtonGroup( 2,Qt::Vertical, tr("EPS options"),this);
    //    EPSGroupBoxLayout = new QVBoxLayout(EPSColorButtonGroupBox);
    //     colorButton = new QRadioButton("Color",EPSColorButtonGroupBox);
    //     BWButton = new QRadioButton("Grayscale",EPSColorButtonGroupBox);
    //     EPSColorButtonGroupBox->setInsideMargin (15); 
    //     EPSColorButtonGroupBox->insert(colorButton);
    //     EPSColorButtonGroupBox->insert(BWButton);
    //     EPSColorButtonGroupBox->setExclusive(true);
    //     EPSWidgetGroupBox->add(EPSColorButtonGroupBox);

    vectorEPSCheckBox = new QCheckBox( "Vector EPS File",EPSWidgetGroupBox);

#else
    QVBoxLayout * EPSGroupBoxLayout = new QVBoxLayout(EPSWidgetGroupBox);
     EPSGroupBoxLayout->setMargin (15); 

//     colorButton = new QRadioButton("Color",EPSWidgetGroupBox);
//     BWButton = new QRadioButton("Grayscale",EPSWidgetGroupBox);

//     QButtonGroup * EPSColorButtonGroupBox = new QButtonGroup();
//     EPSColorButtonGroupBox->addButton(colorButton);
//     EPSColorButtonGroupBox->addButton(BWButton);
//     EPSColorButtonGroupBox->setExclusive(true);

//     EPSGroupBoxLayout->addWidget(colorButton);    
//     EPSGroupBoxLayout->addWidget(BWButton);    

    vectorEPSCheckBox = new QCheckBox( "Vector EPS File",EPSWidgetGroupBox);
    EPSGroupBoxLayout->addWidget(vectorEPSCheckBox);

    EPSWidgetGroupBox->setLayout(EPSGroupBoxLayout);
#endif
    //    colorButton->setChecked( true );
    vectorEPSCheckBox->setChecked( true );
    
    globalVLayout->addWidget(EPSWidgetGroupBox);
    connect( vectorEPSCheckBox, SIGNAL( clicked() ), this, SLOT( changeVectorEPS()) );

  }

  if ((format == "jpg") || 
      (format == "jpeg")) {
    
    QGroupBox *imageGroupBox = new QGroupBox(tr("Image quality"),this);
    QHBoxLayout *hSliderLayout = new QHBoxLayout(imageGroupBox);
    hSliderLayout->setMargin (15); 

    qualitySlider= new QSlider(Qt::Horizontal,imageGroupBox);
#if QT_VERSION < 0x040000
    qualitySlider->setMinValue(0);
    qualitySlider->setMaxValue(100);
    qualitySlider->setTickmarks(QSlider::Below); 
#else
    qualitySlider->setMinimum(0);
    qualitySlider->setMaximum(100);
    qualitySlider->setTickPosition(QSlider::TicksBelow);
#endif
    qualitySlider->setValue(60);
    hSliderLayout->addWidget(new QLabel("Low ",imageGroupBox));
    hSliderLayout->addWidget(qualitySlider);
    hSliderLayout->addWidget(new QLabel(" Maximum",imageGroupBox));
    
#if QT_VERSION >= 0x040000
    imageGroupBox->setLayout(hSliderLayout);
#endif

    globalVLayout->addWidget(imageGroupBox);
  }


  // button ok/cancel box

  QWidget *buttonBox = new QWidget(this);

  QHBoxLayout *buttonBoxLayout = new QHBoxLayout(buttonBox);

  buttonOk = new QPushButton( tr( "&OK" ),buttonBox );
  buttonOk->setAutoDefault( TRUE );
  buttonOk->setDefault( TRUE );
  buttonBoxLayout->addWidget(buttonOk);

  buttonCancel = new QPushButton( tr( "&Cancel" ),buttonBox );
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
  connect( buttonOk, SIGNAL( clicked() ), this, SLOT( accept() ) );
  connect( buttonCancel, SIGNAL( clicked() ), this, SLOT( reject() ) );
}



int G4OpenGLQtExportDialog::getSliderValue()
{
  if (!qualitySlider) return -1;
  return qualitySlider->value();
}

int G4OpenGLQtExportDialog::getHeight()
{
  if (!height) return originalHeight;
  return height->text().toInt();
}

int G4OpenGLQtExportDialog::getWidth()
{
  if (!width) return originalWidth;
  return width->text().toInt();
}

int G4OpenGLQtExportDialog::getTransparency()
{
  if (!boxTransparency) return -1;
  return boxTransparency->isChecked();
}

int G4OpenGLQtExportDialog::getNbColor()
{
  if (!colorButton) return -1;
  // Black and white
  if (!colorButton->isChecked())
    return 1;
  // rgb color
  return 3;
}

bool G4OpenGLQtExportDialog::getVectorEPS()
{
  if (!vectorEPSCheckBox) return 0;
  return vectorEPSCheckBox->isChecked();
}


void G4OpenGLQtExportDialog::changeVectorEPS()
{
  if (!vectorEPSCheckBox) return;
  if (vectorEPSCheckBox->isChecked()) {
#if QT_VERSION < 0x040000
    original->setEnabled ( true );
    modify->setEnabled ( true );
#else
    sizeGroupBox->show();
    original->show();
    modify->show();
#endif
    changeSizeBox();
  } else {
#if QT_VERSION < 0x040000
    original->setEnabled ( false );
    modify->setEnabled ( false );
    ratioCheckBox->setEnabled ( false );
    heightWidget->setEnabled ( false );
    widthWidget->setEnabled ( false );
#else
    sizeGroupBox->hide();
    original->hide();
    modify->hide();
    ratioCheckBox->hide();
    heightWidget->hide();
    widthWidget->hide();
#endif
  }
}


void G4OpenGLQtExportDialog::changeSizeBox()
{
  if (!original) return;
  if (!heightWidget) return;
  if (!widthWidget) return;
  if (!ratioCheckBox) return;

  if ( original->isChecked()) {
#if QT_VERSION < 0x040000
    ratioCheckBox->setEnabled ( false );
    heightWidget->setEnabled ( false );
    widthWidget->setEnabled ( false );
#else
    ratioCheckBox->hide();
    heightWidget->hide();
    widthWidget->hide();
#endif
  } else {
#if QT_VERSION < 0x040000
    ratioCheckBox->setEnabled ( true );
    heightWidget->setEnabled ( true );
    widthWidget->setEnabled ( true );
#else
    heightWidget->show();
    widthWidget->show();
    ratioCheckBox->show();
#endif
  }
}


void G4OpenGLQtExportDialog::textWidthChanged(
 const QString & s
 )
{
  if (!ratioCheckBox) return;
  if (!width) return;
  if (isChangingSize == true) return; // exclusive slot

  if (ratioCheckBox->isChecked()){
    isChangingSize = true;
    QString tmp;
  height->setText(tmp.setNum((int)(s.toInt()*(double)((double)originalHeight/(double)originalWidth))));
  isChangingSize = false;
  }
}

void G4OpenGLQtExportDialog::  textHeightChanged(
 const QString & s
)
{
  if (!ratioCheckBox) return;
  if (!width) return;
  if (isChangingSize == true) return; // exclusive slot

  if (ratioCheckBox->isChecked()){
  isChangingSize = true;
    QString tmp;
    width->setText(tmp.setNum(s.toInt()*originalWidth/originalHeight));
  isChangingSize = false;
  }
} 

G4OpenGLQtExportDialog::~G4OpenGLQtExportDialog()
{
}


#endif
