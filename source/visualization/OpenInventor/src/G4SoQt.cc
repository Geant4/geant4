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
// F.W. Jones 05012018

#include <stdlib.h>
#include <string.h>

#include "G4ios.hh"

#include "G4SoQt.hh"
//#include "G4Qt.hh"
//#include "G4UImanager.hh"

#include <Inventor/Qt/SoQt.h>
#include <QMainWindow>

#include <qwidget.h>
#include <qapplication.h>

#ifndef G4GMAKE
#include "moc_G4SoQt.cpp"
#endif

#define G4warn G4cout

G4SoQt* G4SoQt::instance = NULL;

static G4bool QtInited = FALSE;


/***************************************************************************/
G4SoQt* G4SoQt::getInstance()
{
  if (instance==NULL) {
     instance = new G4SoQt();
  }
  return instance;
}


/***************************************************************************/
G4SoQt::G4SoQt()
// FWJ command-line arguments omitted for time being
//G4SoQt* G4SoQt(int a_argn, char** a_args, char*)
{
  externalApp = false;

  // FWJ in G4Qt this is used.  qApp is Qt global pointer
  // to the (one and only one) QApplication object:
  //      new QApplication (*p_argn, args);
  //      if(!qApp) {

  // FWJ detects existence of Qt UI or other running qApp
  if (qApp) externalApp = true;

  QWidget* mainWin = SoQt::init("Geant4");

  // FWJ Cf Xt:
  QtInited = TRUE;

  // FWJ DEBUG
  //  G4cout << "G4SoQt: mainWin=" << mainWin << G4endl;
  //  G4cout << "G4SoQt: toplevelwidget=" << SoQt::getTopLevelWidget() << G4endl;

// FWJ CAN'T GET MENUBAR THIS WAY OR BY CAST
  //  QWidget* toplevel = SoQt::getTopLevelWidget();
  //  QMenuBar* menubar = QMainWindow::menuBar();
  //  G4cout << "G4OpenInventorQtExaminerViewer menubar=" << menubar << G4endl;


  // FWJ will this work?
  SetMainInteractor(mainWin);
  //  SetMainInteractor(qApp);
  //  AddDispatcher     ((G4DispatchFunction)XtDispatchEvent);

// FWJ no locale for now (see G4Qt.cc)

}

/***************************************************************************/
G4SoQt::~G4SoQt()
{
  if(this==instance) {
    instance = NULL;
  }
}

G4bool G4SoQt::Inited()
{
  return QtInited;
}

/***************************************************************************/
void* G4SoQt::GetEvent()
{
  return 0;
}


/***************************************************************************/
void G4SoQt::FlushAndWaitExecution()
{
   // FWJ the following is used in G4Qt:
   //  if(!qApp) return;
   // This starts the Qt main loop:
   //  qApp->processEvents();

   // FWJ no, should be done in secondaryLoop()!
   //   SoQt::mainLoop();
}


/***************************************************************************/
void G4SoQt::SecondaryLoop()
{
   if (externalApp) return;

   // FWJ DEBUG
   //      G4cout <<
   //     "ENTERING OIQT VIEWER SECONDARY LOOP" << G4endl;
   //   else

   G4warn <<
      "ENTERING OIQT VIEWER SECONDARY LOOP... PRESS E KEY TO EXIT" << G4endl;

   SoQt::mainLoop();
}


/***************************************************************************/
void G4SoQt::ExitSecondaryLoop()
{
   // FWJ DEBUG
   //   G4cout << "G4SoQt: EXIT SECONDARY LOOP externalApp=" <<
   //      externalApp << G4endl;

   if (externalApp) return;   
   SoQt::exitMainLoop();
}
/***************************************************************************/
bool G4SoQt::IsExternalApp()
{
  return externalApp;
}
