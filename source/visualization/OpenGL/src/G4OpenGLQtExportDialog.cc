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
// $Id: G4OpenGLQtExportDialog.cc,v 1.3 2007-11-09 15:03:21 lgarnier Exp $
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

G4OpenGLQtExportDialog::G4OpenGLQtExportDialog(
 QWidget* parent
,QString nomFich
 ,int aHeight
 ,int aWidth
)
  : QDialog( parent )
{
  setWindowTitle( tr( " Export options" ));
  originalWidth = aWidth;
  originalHeight = aHeight;

  // Initializations
  qualitySlider = NULL;
  width = NULL;
  height = NULL;
  colorButton = NULL;
  BWButton = NULL;

  // global layout
  QVBoxLayout* globalVLayout = new QVBoxLayout();

  
  if (nomFich.endsWith(".jpg") || 
      nomFich.endsWith(".jepg")) {
    
    QGroupBox *imageGroupBox = new QGroupBox(tr("Image quality"));
    QVBoxLayout *imageGroupBoxLayout = new QVBoxLayout;

    QWidget *sliderBox = new QWidget;
    QHBoxLayout *hSlider = new QHBoxLayout;
    //    qualityLabel =  new QLabel( tr( "Image quality" ) );
    //    imageGroupBoxLayout->addWidget(qualityLabel);
    qualitySlider= new QSlider(Qt::Horizontal);
    qualitySlider->setMinimum(0);
    qualitySlider->setMaximum(100);
    qualitySlider->setTickPosition(QSlider::TicksBelow);
    qualitySlider->setValue(60);
    hSlider->addWidget(new QLabel("low"));
    hSlider->addWidget(qualitySlider);
    hSlider->addWidget(new QLabel("Maximum"));
    sliderBox->setLayout(hSlider);
    imageGroupBoxLayout->addWidget(sliderBox);

    imageGroupBox->setLayout(imageGroupBoxLayout);
    globalVLayout->addWidget(imageGroupBox);
  }
  
  if(nomFich.endsWith(".eps")) {
    QGroupBox *EPSGroupBox = new QGroupBox(tr("EPS options"));
    QVBoxLayout *EPSGroupBoxLayout = new QVBoxLayout;

    //    transparencyEPS = new QCheckBox( "transparencyEPS" );
    //    transparencyEPS->setText( "save background" );
    //    transparencyEPS->setChecked( true );

    colorButton = new QRadioButton("Color");
    BWButton = new QRadioButton("Grayscale");
    colorButton->setChecked( true );
    BWButton->setChecked( false );


    //    EPSGroupBoxLayout->addWidget(transparencyEPS);    
    EPSGroupBoxLayout->addWidget(colorButton);    
    EPSGroupBoxLayout->addWidget(BWButton);    
    EPSGroupBox->setLayout(EPSGroupBoxLayout);
    globalVLayout->addWidget(EPSGroupBox);

  }

  if(nomFich.endsWith(".tif") ||
     nomFich.endsWith(".tiff") ||
     nomFich.endsWith(".jpg") ||
     nomFich.endsWith(".png") ||
     nomFich.endsWith(".xpm")) {

    QGroupBox *transparencyGroupBox = new QGroupBox(tr("Transparency"));
    QVBoxLayout *transparencyGroupBoxLayout = new QVBoxLayout;

    boxTransparency = new QCheckBox("Save transparency");
    boxTransparency->setChecked( false );
    //    boxTransparency->setEnabled(false);

    transparencyGroupBoxLayout->addWidget(boxTransparency);    
    transparencyGroupBox->setLayout(transparencyGroupBoxLayout);
    globalVLayout->addWidget(transparencyGroupBox);

  }

  // size box
  QGroupBox *sizeGroupBox = new QGroupBox(tr("Size"));
  QVBoxLayout *sizeGroupBoxLayout = new QVBoxLayout;
  
  QHBoxLayout *modifyAndRatioLayout = new QHBoxLayout;
  QWidget* modifyAndRatioWidget = new QWidget;

  // original button
  original = new QRadioButton("Original");
  original->setChecked( true );
  sizeGroupBoxLayout->addWidget(original);

  // modify and ratio
  modify = new QRadioButton("Modify");
  modify->setChecked( false );

  ratioCheckBox = new QCheckBox( "Keep ratio" );
  ratioCheckBox->setChecked( true );

  modifyAndRatioLayout->addWidget(modify);
  modifyAndRatioLayout->addWidget(ratioCheckBox);
  modifyAndRatioWidget->setLayout(modifyAndRatioLayout);
  sizeGroupBoxLayout->addWidget(modifyAndRatioWidget);
  ratioCheckBox->setVisible(modify->isChecked());

  connect( original, SIGNAL( clicked(bool) ), this, SLOT( changeSizeBox(true)) );
  connect( modify, SIGNAL( clicked(bool) ), this, SLOT( changeSizeBox(false) ) );

  // height
  QHBoxLayout *heightLineLayout = new QHBoxLayout;
  heightWidget = new QWidget;
  QString tmp;
 
  heightLineLayout->addWidget(new QLabel("Height"));
  height = new QLineEdit(tmp.setNum(originalHeight));
  height->setMaxLength(5);
  heightLineLayout->addWidget(height);
  heightWidget->setLayout(heightLineLayout);
  sizeGroupBoxLayout->addWidget(heightWidget);
  connect( height, SIGNAL( textChanged ( const QString& ) ), this, SLOT( textHeightChanged(const QString &) ) );


  // width
  QHBoxLayout *widthLineLayout = new QHBoxLayout;
  widthWidget = new QWidget;

  widthLineLayout->addWidget(new QLabel("Width "));
  width = new QLineEdit(tmp.setNum(originalWidth));
  width->setMaxLength(5);
  widthLineLayout->addWidget(width);
  widthWidget->setLayout(widthLineLayout);
  sizeGroupBoxLayout->addWidget(widthWidget);
  connect( width, SIGNAL( textChanged ( const QString& ) ), this, SLOT( textWidthChanged(const QString &) ) );

  sizeGroupBox->setLayout(sizeGroupBoxLayout);
  globalVLayout->addWidget(sizeGroupBox);

  heightWidget->setVisible(false);
  widthWidget->setVisible(false);

  // button ok/cancel box

  QGroupBox *buttonGroupBox = new QGroupBox();
  QHBoxLayout *buttonGroupBoxLayout = new QHBoxLayout;

  buttonOk = new QPushButton( tr( "&OK" ) );
  buttonOk->setAutoDefault( TRUE );
  buttonOk->setDefault( TRUE );
  buttonGroupBoxLayout->addWidget(buttonOk);

  buttonCancel = new QPushButton( tr( "&Cancel" ) );
  buttonCancel->setAutoDefault( TRUE );
  buttonGroupBoxLayout->addWidget(buttonCancel);

  buttonGroupBox->setLayout(buttonGroupBoxLayout);
  globalVLayout->addWidget(buttonGroupBox);


  setLayout(globalVLayout);

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
    heightWidget->setVisible(false);
    widthWidget->setVisible(false);
    ratioCheckBox->setVisible(false);
  } else {
    heightWidget->setVisible(true);
    widthWidget->setVisible(true);
    ratioCheckBox->setVisible(true);
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
