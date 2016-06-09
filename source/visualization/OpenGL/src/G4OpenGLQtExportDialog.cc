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
// $Id: G4OpenGLQtExportDialog.cc,v 1.4 2007/11/13 17:48:51 lgarnier Exp $
// GEANT4 tag $Name: geant4-09-01 $
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

G4OpenGLQtExportDialog::G4OpenGLQtExportDialog(
 QWidget* parent
,QString nomFich
 ,int aHeight
 ,int aWidth
)
  : QDialog( parent )
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
#if QT_VERSION < 0x040000
  QVBoxLayout* globalVLayout = new QVBoxLayout(this);
#else
  QVBoxLayout* globalVLayout = new QVBoxLayout();
#endif

  
  if (nomFich.endsWith(".jpg") || 
      nomFich.endsWith(".jepg")) {
    
    QGroupBox *imageGroupBox = new QGroupBox(tr("Image quality"));
#if QT_VERSION < 0x040000
    QVBoxLayout *imageGroupBoxLayout = new QVBoxLayout(imageGroupBox);
#else
    QVBoxLayout *imageGroupBoxLayout = new QVBoxLayout;
#endif
    QWidget *sliderBox = new QWidget;

#if QT_VERSION < 0x040000
    QHBoxLayout *hSlider = new QHBoxLayout(sliderBox);
#else
    QHBoxLayout *hSlider = new QHBoxLayout;
#endif

    //    qualityLabel =  new QLabel( tr( "Image quality" ) );
    //    imageGroupBoxLayout->addWidget(qualityLabel);
    qualitySlider= new QSlider(Qt::Horizontal,0);
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
    hSlider->addWidget(new QLabel("low",0));
    hSlider->addWidget(qualitySlider);
    hSlider->addWidget(new QLabel("Maximum",0));
#if QT_VERSION >= 0x040000
    sliderBox->setLayout(hSlider);
#endif
    imageGroupBoxLayout->addWidget(sliderBox);

#if QT_VERSION >= 0x040000
    imageGroupBox->setLayout(imageGroupBoxLayout);
#endif
    globalVLayout->addWidget(imageGroupBox);
  }
  
  if(nomFich.endsWith(".eps")) {
    QGroupBox *EPSGroupBox = new QGroupBox(tr("EPS options"));

#if QT_VERSION < 0x040000
    QVBoxLayout *EPSGroupBoxLayout = new QVBoxLayout(EPSGroupBox);
#else
    QVBoxLayout *EPSGroupBoxLayout = new QVBoxLayout;
#endif

    //    transparencyEPS = new QCheckBox( "transparencyEPS" );
    //    transparencyEPS->setText( "save background" );
    //    transparencyEPS->setChecked( true );

    colorButton = new QRadioButton("Color",0);
    BWButton = new QRadioButton("Grayscale",0);
    colorButton->setChecked( true );
    BWButton->setChecked( false );


    //    EPSGroupBoxLayout->addWidget(transparencyEPS);    
    EPSGroupBoxLayout->addWidget(colorButton);    
    EPSGroupBoxLayout->addWidget(BWButton);    
#if QT_VERSION >= 0x040000
    EPSGroupBox->setLayout(EPSGroupBoxLayout);
#endif
    globalVLayout->addWidget(EPSGroupBox);

  }

  if(nomFich.endsWith(".tif") ||
     nomFich.endsWith(".tiff") ||
     nomFich.endsWith(".jpg") ||
     nomFich.endsWith(".png") ||
     nomFich.endsWith(".xpm")) {

    QGroupBox *transparencyGroupBox = new QGroupBox(tr("Transparency"));
#if QT_VERSION < 0x040000
    QVBoxLayout *transparencyGroupBoxLayout = new QVBoxLayout(transparencyGroupBox);
#else
    QVBoxLayout *transparencyGroupBoxLayout = new QVBoxLayout;
#endif

    boxTransparency = new QCheckBox("Save transparency",0);
    boxTransparency->setChecked( false );
    //    boxTransparency->setEnabled(false);

    transparencyGroupBoxLayout->addWidget(boxTransparency);    
#if QT_VERSION >= 0x040000
    transparencyGroupBox->setLayout(transparencyGroupBoxLayout);
#endif
    globalVLayout->addWidget(transparencyGroupBox);

  }

  // size box
  QGroupBox *sizeGroupBox = new QGroupBox(tr("Size"));
  QWidget* modifyAndRatioWidget = new QWidget;

#if QT_VERSION < 0x040000
  QHBoxLayout *modifyAndRatioLayout = new QHBoxLayout(modifyAndRatioWidget);
  QVBoxLayout *sizeGroupBoxLayout = new QVBoxLayout(sizeGroupBox);
#else
  QHBoxLayout *modifyAndRatioLayout = new QHBoxLayout;
  QVBoxLayout *sizeGroupBoxLayout = new QVBoxLayout;
#endif

  // original button
  original = new QRadioButton("Original",0);
  original->setChecked( true );
  sizeGroupBoxLayout->addWidget(original);

  // modify and ratio
  modify = new QRadioButton("Modify",0);
  modify->setChecked( false );

  ratioCheckBox = new QCheckBox( "Keep ratio",0 );
  ratioCheckBox->setChecked( true );

  modifyAndRatioLayout->addWidget(modify);
  modifyAndRatioLayout->addWidget(ratioCheckBox);
#if QT_VERSION >= 0x040000
  modifyAndRatioWidget->setLayout(modifyAndRatioLayout);
#endif
  sizeGroupBoxLayout->addWidget(modifyAndRatioWidget);
  if (modify->isChecked()) {
    ratioCheckBox->show();
  } else {
    ratioCheckBox->hide();
  }

  connect( original, SIGNAL( clicked(bool) ), this, SLOT( changeSizeBox(true)) );
  connect( modify, SIGNAL( clicked(bool) ), this, SLOT( changeSizeBox(false) ) );

  // height
  heightWidget = new QWidget;

#if QT_VERSION < 0x040000
  QHBoxLayout *heightLineLayout = new QHBoxLayout(heightWidget);
#else
  QHBoxLayout *heightLineLayout = new QHBoxLayout;
#endif

  QString tmp;
 
  heightLineLayout->addWidget(new QLabel("Height",0));
  height = new QLineEdit(tmp.setNum(originalHeight),0);
  height->setMaxLength(5);
  heightLineLayout->addWidget(height);
#if QT_VERSION >= 0x040000
  heightWidget->setLayout(heightLineLayout);
#endif
  sizeGroupBoxLayout->addWidget(heightWidget);
  connect( height, SIGNAL( textChanged ( const QString& ) ), this, SLOT( textHeightChanged(const QString &) ) );


  // width
  widthWidget = new QWidget;

#if QT_VERSION < 0x040000
  QHBoxLayout *widthLineLayout = new QHBoxLayout(widthWidget);
#else
  QHBoxLayout *widthLineLayout = new QHBoxLayout;
#endif

  widthLineLayout->addWidget(new QLabel("Width ",0));
  width = new QLineEdit(tmp.setNum(originalWidth),0);
  width->setMaxLength(5);
  widthLineLayout->addWidget(width);
#if QT_VERSION >= 0x040000
  widthWidget->setLayout(widthLineLayout);
#endif
  sizeGroupBoxLayout->addWidget(widthWidget);
  connect( width, SIGNAL( textChanged ( const QString& ) ), this, SLOT( textWidthChanged(const QString &) ) );

#if QT_VERSION >= 0x040000
  sizeGroupBox->setLayout(sizeGroupBoxLayout);
#endif
  globalVLayout->addWidget(sizeGroupBox);

  heightWidget->hide();
  widthWidget->hide();

  // button ok/cancel box

  QGroupBox *buttonGroupBox = new QGroupBox();

#if QT_VERSION < 0x040000
  QHBoxLayout *buttonGroupBoxLayout = new QHBoxLayout(buttonGroupBox);
#else
  QHBoxLayout *buttonGroupBoxLayout = new QHBoxLayout;
#endif

  buttonOk = new QPushButton( tr( "&OK" ),0 );
  buttonOk->setAutoDefault( TRUE );
  buttonOk->setDefault( TRUE );
  buttonGroupBoxLayout->addWidget(buttonOk);

  buttonCancel = new QPushButton( tr( "&Cancel" ),0 );
  buttonCancel->setAutoDefault( TRUE );
  buttonGroupBoxLayout->addWidget(buttonCancel);

#if QT_VERSION >= 0x040000
  buttonGroupBox->setLayout(buttonGroupBoxLayout);
#endif
  globalVLayout->addWidget(buttonGroupBox);


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
  if (!height) return -1;
  return height->text().toInt();
}

int G4OpenGLQtExportDialog::getWidth()
{
  if (!width) return -1;
  return width->text().toInt();
}

bool G4OpenGLQtExportDialog::getTransparency()
{
  if (!boxTransparency) return -1;
  return boxTransparency->isChecked();
}

int G4OpenGLQtExportDialog::getNbColor()
{
  // Black and white
  if (!colorButton->isChecked())
    return 1;
  // rgb color
  return 3;
}


void G4OpenGLQtExportDialog::changeSizeBox(bool aClick)
{
  if (aClick) {
    modify->toggle();
  } else {
    original->toggle();
  }
  if ( original->isChecked()) {
    heightWidget->hide();
    widthWidget->hide();
    ratioCheckBox->hide();
  } else {
    heightWidget->show();
    widthWidget->show();
    ratioCheckBox->show();
  }
}

void G4OpenGLQtExportDialog::textWidthChanged(
 const QString & s
 )
{
  if (ratioCheckBox->isChecked()){
    QString tmp;
    width->setText(tmp.setNum(s.toInt()*originalHeight/originalHeight));
  }
}

void G4OpenGLQtExportDialog::  textHeightChanged(
 const QString & s
)
{
  if (ratioCheckBox->isChecked()){
    QString tmp;
    width->setText(tmp.setNum(s.toInt()*originalWidth/originalWidth));
  }
} 

G4OpenGLQtExportDialog::~G4OpenGLQtExportDialog()
{
}


#endif
