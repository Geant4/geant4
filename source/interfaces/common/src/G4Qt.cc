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
// $Id: G4Qt.cc 97174 2016-05-27 12:41:02Z gcosmo $
//
// L. Garnier

#if defined(G4INTY_BUILD_QT) || defined(G4INTY_USE_QT)

#include <stdlib.h>
#include <string.h>

#include "G4ios.hh"

#include "G4Qt.hh"
#include "G4UImanager.hh"
#include <qwidget.h>

#include <qapplication.h>


G4Qt* G4Qt::instance    = NULL;

static G4bool QtInited  = FALSE;

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
 ,char*  /*a_class */
 )
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  argn = 0;
  args = NULL;
  externalApp = false;

  // Check if Qt already init in another external app
  if(qApp) {
    externalApp = true;
    QtInited  = TRUE;
    SetMainInteractor (qApp);
    SetArguments      (a_argn,a_args);
    
  } else {
    
    if(QtInited==FALSE) {  //Qt should be Inited once !
      // Then two cases :
      // - It is the first time we create G4UI  (argc!=0)
      //   -> Inited and register
      // - It is the first time we create G4VIS  (argc == 0)
      //   -> Inited and NOT register
      
      if (a_argn != 0) {
        argn = a_argn;
        args = a_args;

      } else { //argc = 0

        // FIXME : That's not the good arguments, but I don't know how to get args from other Interactor.
        // Ex: How to get them from G4Xt ?
        argn = 1;
        args = (char **)malloc( 1 * sizeof(char *) );
        args[0] = (char *)malloc(10 * sizeof(char));
        strncpy(args[0], "my_app \0", 9);
      }

      int *p_argn = (int*)malloc(sizeof(int));
      *p_argn = argn;
      new QApplication (*p_argn, args);
      if(!qApp) {
        
        G4UImanager* UImanager = G4UImanager::GetUIpointer();
        G4int verbose = UImanager->GetVerboseLevel();
        if (verbose >= 2) {
          G4cout        << "G4Qt : Unable to init Qt." << G4endl;
        }
      } else {
        QtInited  = TRUE;
        if (a_argn != 0) {
          SetMainInteractor (qApp);
        }
        SetArguments      (a_argn,a_args);
      }
    }
  }
  //  AddDispatcher     ((G4DispatchFunction)XtDispatchEvent);
  
  /*
   * On some non-English locale, comma is used for the decimal separator instead of dot
   * bringing to weird behavior of strtod (string to double) function in user application.
   * This is "by design" from Qt, see https://bugreports.qt-project.org/browse/QTBUG-10994
   *
   *      $ LC_NUMERIC=fr_FR.UTF-8 ./qtstrtod
   *      strtod(0.1) = 0
   *      $ LC_NUMERIC=C ./qtstrtod
   *      strtod(0.1) = 0.1
   *
   * Jerome Suhard, jerome@suhard.fr
   */

  // explicitly set the LC_NUMBERIC locale to "C"
  setlocale (LC_NUMERIC, "C");
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
void* G4Qt::GetEvent (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return 0;
}
/***************************************************************************/
void G4Qt::FlushAndWaitExecution (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(!qApp) return;
  qApp->processEvents();
}

/***************************************************************************/
bool G4Qt::IsExternalApp (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return externalApp;
}

#endif


