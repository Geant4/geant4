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
// $Id: G4Qt.cc,v 1.7 2007/11/15 18:24:28 lgarnier Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// L. Garnier

#if defined(G4INTY_BUILD_QT) || defined(G4INTY_USE_QT)

#include <stdlib.h>
#include <string.h>

#include "G4ios.hh"

#include "G4Qt.hh"

#include <qapplication.h>


#define NewString(str)  \
 ((str) != NULL ? (strcpy((char*)malloc((unsigned)strlen(str) + 1), str)) : NULL)

//static void XWidgetIconify                 (Widget);
//static void XWidgetUniconify               (Widget);
//static void XDisplaySetWindowToNormalState (Display*,Window);

G4Qt* G4Qt::instance    = NULL;

static G4bool QtInited  = FALSE;
//static int    argn      = 0;
//static char** args      = NULL;
// static QtAppContext appContext = NULL;
//static QApplication app = NULL;

/***************************************************************************/
G4Qt* G4Qt::getInstance (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return G4Qt::getInstance (0,NULL,(char*)"Geant4");
}
/***************************************************************************/
G4Qt* G4Qt::getInstance (
 int    a_argn
,char** a_args
,char*  a_class
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if (instance==NULL) {
    instance = new G4Qt(a_argn,a_args,a_class);
  }
  return instance;
}
/***************************************************************************/
G4Qt::G4Qt (
 int    a_argn
,char** a_args
,char*  a_class
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef GEANT4_QT_DEBUG
  printf("G4Qt::G4Qt try to inited Qt\n");
#endif
  if(QtInited==FALSE) {  //Qt should be Inited once !
#ifdef GEANT4_QT_DEBUG
    printf("G4Qt::G4Qt inited Qt\n");
#endif
#if QT_VERSION < 0x040000
    qApp = new QApplication (a_argn, a_args);
    //    QApplication qApp(a_argn, a_args);
    //    if(&qApp == NULL) {
#else
    new QApplication (a_argn, a_args);
#endif
    if(!qApp) {

      G4cout        << "G4Qt : Unable to init Qt." << G4endl;
    } else {
      QtInited  = TRUE;
      //#if QT_VERSION < 0x040000
      //      SetMainInteractor (&qApp);
      //#else
      SetMainInteractor (qApp);
      //#endif
      SetArguments      (a_argn,a_args);
#ifdef GEANT4_QT_DEBUG
      printf("G4Qt::G4Qt inited Qt END\n");
#endif
    }
  }
#ifdef GEANT4_QT_DEBUG
  if (qApp) {
    printf("G4Qt::qApp exist\n");
  }  else {
    printf("G4Qt::qApp not exist\n");
  }
#endif
  //  AddDispatcher     ((G4DispatchFunction)XtDispatchEvent);
}
/***************************************************************************/
G4Qt::~G4Qt (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(this==instance) {
    instance = NULL;
  }
}
/***************************************************************************/
G4bool G4Qt::Inited (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return QtInited;
}
/***************************************************************************/
/**
  Si j'ai bien compris, cette fonction ne sert à rien
 */
void* G4Qt::GetEvent (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
//FIXME
//   G4cout        << "G4Qt : Rien compris a cette fonction G4Qt::GetEvent." << G4endl;
//  static XEvent  event;
//  if(appContext==NULL) return NULL;
//  if(mainApp==NULL) return NULL;
//  QtAppNextEvent (appContext, &event);
//  return         &event;
  printf("*");
  return 0;
}
/***************************************************************************/
void G4Qt::FlushAndWaitExecution (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  //  printf("G4Qt::FlushAndWaitExecution ::  Flush ....\n");
  if(!qApp) return;
  qApp->processEvents();
}

#endif



