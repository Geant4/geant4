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

#include "G4RTQtScanner.hh"

//#include "G4TheRayTracer.hh"
//#include "G4RayTracerQtViewer.hh"
#include "G4ViewParameters.hh"
#include "G4UImanager.hh"
#include "G4UIQt.hh"

#include <QLabel>
#include <QImage>
#include <QPainter>
#include <QColor>

#define G4warn G4cout

G4RTQtScanner::G4RTQtScanner()
: G4VRTScanner()
, theNRow(0), theNColumn(0)
, theIRow(0), theIColumn(0)
, fpImageLabel(nullptr)
, fpImage(nullptr)
{}

G4RTQtScanner::~G4RTQtScanner() {}

void G4RTQtScanner::Initialize(G4int nRow, G4int nColumn) {
  theNRow = nRow;
  theNColumn = nColumn;
  theIRow = 0;
  theIColumn = -1;
}

G4bool G4RTQtScanner::Coords(G4int& iRow, G4int& iColumn)
{
  // Increment column and, if necessary, increment row...
  ++theIColumn;
  if (theIColumn >= theNColumn) {
    theIColumn = 0;
    ++theIRow;
  }

  // Return if finished...
  if (theIRow >= theNRow) {
    // ...and paint the image
    fpImageLabel->setPixmap(QPixmap::fromImage(*fpImage));
    QPainter windowPainter(fpImageLabel);
    windowPainter.drawImage(0, 0, *fpImage);
    return false;
  }

  // Return current row and column...
  iRow = theIRow;
  iColumn = theIColumn;
  return true;
}

G4bool G4RTQtScanner::GetQtWindow(const G4String& name, G4ViewParameters& vp)
{
  auto UI = G4UImanager::GetUIpointer();
  auto uiQt = dynamic_cast<G4UIQt*>(UI->GetG4UIWindow());
  if (!uiQt) {
    G4warn << "G4RTQtScanner::GetQtWindow: RayTracerQt requires G4UIQt"
    << G4endl;
    return false;
  }
  uiQt->AddTabWidget(this,QString(name));
  setBackgroundRole(QPalette::Dark);

  theWindowSizeX = vp.GetWindowSizeHintX();  // As used by the ray tracer
  theWindowSizeY = vp.GetWindowSizeHintY();  // and to make fpImage
  fpImage = new QImage(theWindowSizeX, theWindowSizeY, QImage::Format_RGB32);
  fpImageLabel = new QLabel;
  fpImageLabel->setPixmap(QPixmap::fromImage(*fpImage));
  setWidget(fpImageLabel);
  setAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
  ensureVisible(theWindowSizeX/2, theWindowSizeY/2);

  return true;
}

void G4RTQtScanner::Draw
(unsigned char red, unsigned char green, unsigned char blue)
// Add a coloured point to the image at current position.
{
  QPainter painter(fpImage);
  painter.setPen(QPen(QColor(int(red), int(green), int(blue))));
  painter.drawPoint(theIColumn, theIRow);
}
