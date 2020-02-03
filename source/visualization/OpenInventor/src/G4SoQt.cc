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

#if defined(G4VIS_BUILD_OIQT_DRIVER)
//#if defined(G4INTY_BUILD_QT) || defined(G4INTY_USE_QT)

#include <stdlib.h>
#include <string.h>

#include "G4ios.hh"

#include "G4SoQt.hh"
//#include "G4Qt.hh"
//#include "G4UImanager.hh"

#include <Inventor/Qt/SoQt.h>

#include <qwidget.h>
#include <qapplication.h>


G4SoQt* G4SoQt::instance    = NULL;

static G4bool QtInited  = FALSE;

/***************************************************************************/
G4SoQt* G4SoQt::getInstance() 
{
  if (instance==NULL) {
     instance = new G4SoQt();
  }
  return instance;
  //  return G4SoQt::getInstance(0, NULL, (char*)"Geant4");
}

/***************************************************************************/
//G4SoQt* G4SoQt::getInstance(int a_argn, char** a_args, char* a_class)
//{
//  if (instance==NULL) {
//    instance = new G4SoQt();
//  }
//  return instance;
//}

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
        
  QWidget* mainWin = SoQt::init("Geant4");

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
   G4cout << "G4SoQt: SECONDARY LOOP CALLED !!!!!!" << G4endl;
   SoQt::mainLoop();
}


/***************************************************************************/
bool G4SoQt::IsExternalApp()
{
  return externalApp;
}

#endif
