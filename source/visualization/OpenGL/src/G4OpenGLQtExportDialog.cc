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
// $Id: G4OpenGLQtExportDialog.cc 82764 2014-07-08 14:24:04Z gcosmo $
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
 QWidget* parentw
,QString format
 ,int aHeight
 ,int aWidth
)
  : QDialog( parentw ),
    isChangingSize(false)
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
  QVBoxLayout* globalVLayout = new QVBoxLayout(this);
  globalVLayout->setMargin(10);
  globalVLayout->setSpacing(10);
  



  // size box

  QWidget * sizeWidget = new QWidget(this); // widget containing group button
  QVBoxLayout * sizeWidgetLayout = new QVBoxLayout(sizeWidget);
  sizeWidgetLayout->setMargin (10); 

  // original and modify radiobuttons
  
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
  original->setChecked( true );


  // height
  heightWidget = new QWidget(sizeWidget);

  QHBoxLayout *heightLineLayout = new QHBoxLayout(heightWidget);

  QString tmp;
 
  heightLineLayout->addWidget(new QLabel("Height",heightWidget));
  height = new QLineEdit(tmp.setNum(originalHeight),heightWidget);
  height->setMaxLength(5);
  heightLineLayout->addWidget(height);

  heightWidget->setLayout(heightLineLayout);

  sizeWidgetLayout->addWidget(heightWidget);
  connect( height, SIGNAL( textChanged ( const QString& ) ), this, SLOT( textHeightChanged(const QString &) ) );


  // width
  widthWidget = new QWidget(sizeWidget);

  QHBoxLayout *widthLineLayout = new QHBoxLayout(widthWidget);

  widthLineLayout->addWidget(new QLabel("Width ",widthWidget));
  width = new QLineEdit(tmp.setNum(originalWidth),widthWidget);
  width->setMaxLength(5);
  widthLineLayout->addWidget(width);
  widthWidget->setLayout(widthLineLayout);
  sizeWidgetLayout->addWidget(widthWidget);
  connect( width, SIGNAL( textChanged ( const QString& ) ), this, SLOT( textWidthChanged(const QString &) ) );



  // ratio check box

  ratioCheckBox = new QCheckBox( "Keep ratio",sizeWidget);
  ratioCheckBox->setChecked( true );

  sizeWidgetLayout->addWidget(ratioCheckBox);

  ratioCheckBox->hide();
  heightWidget->hide();
  widthWidget->hide();

  sizeWidget->setLayout(sizeWidgetLayout);
  globalVLayout->addWidget(sizeWidget);

 if (format == "eps") {

   QGroupBox *EPSWidgetGroupBox = new QGroupBox(tr("EPS options"),this); // widget containing group button


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
    qualitySlider->setMinimum(0);
    qualitySlider->setMaximum(100);
    qualitySlider->setTickPosition(QSlider::TicksBelow);
    qualitySlider->setValue(60);
    hSliderLayout->addWidget(new QLabel("Low ",imageGroupBox));
    hSliderLayout->addWidget(qualitySlider);
    hSliderLayout->addWidget(new QLabel(" Maximum",imageGroupBox));
    
    imageGroupBox->setLayout(hSliderLayout);

    globalVLayout->addWidget(imageGroupBox);
  }


  // button ok/cancel box

  QWidget *buttonBox = new QWidget(this);

  QHBoxLayout *buttonBoxLayout = new QHBoxLayout(buttonBox);

  buttonOk = new QPushButton( tr( "&OK" ),buttonBox );
  buttonOk->setAutoDefault( true );
  buttonOk->setDefault( true );
  buttonBoxLayout->addWidget(buttonOk);

  buttonCancel = new QPushButton( tr( "&Cancel" ),buttonBox );
  buttonCancel->setAutoDefault( true );
  buttonBoxLayout->addWidget(buttonCancel);

  buttonBox->setLayout(buttonBoxLayout);
  globalVLayout->addWidget(buttonBox);



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
  if (!height) return originalHeight;
  return height->text().toInt();
}

int G4OpenGLQtExportDialog::getWidth()
{
  if (!width) return originalWidth;
  return width->text().toInt();
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
    sizeGroupBox->show();
    original->show();
    modify->show();
    changeSizeBox();
  } else {
    sizeGroupBox->hide();
    original->hide();
    modify->hide();
    ratioCheckBox->hide();
    heightWidget->hide();
    widthWidget->hide();
  }
}


void G4OpenGLQtExportDialog::changeSizeBox()
{
  if (!original) return;
  if (!heightWidget) return;
  if (!widthWidget) return;
  if (!ratioCheckBox) return;

  if ( original->isChecked()) {
    ratioCheckBox->hide();
    heightWidget->hide();
    widthWidget->hide();
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
