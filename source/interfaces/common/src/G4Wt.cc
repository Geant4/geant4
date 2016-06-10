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
// $Id: G4Wt.cc,v 1.17 2010-05-20 07:01:03 lgarnier Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// L. Garnier

#if defined(G4INTY_BUILD_WT) || defined(G4INTY_USE_WT)

#include <stdlib.h>
#include <string.h>

#include "G4ios.hh"

#include "G4Wt.hh"

#include <Wt/WApplication>
#include <Wt/WEnvironment>


G4Wt* G4Wt::instance    = NULL;

static G4bool WtInited  = FALSE;

/***************************************************************************/
G4Wt* G4Wt::getInstance (
 ) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return G4Wt::getInstance (0,NULL,(char*)"Geant4");
}

/***************************************************************************/
G4Wt* G4Wt::getInstance (
 int    a_argn
,char** a_args
,char*  a_class
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if (instance==NULL) {
    instance = new G4Wt(a_argn,a_args,a_class);
  }
  return instance;
}
/***************************************************************************/
G4Wt::G4Wt (
 int    a_argn
,char** a_args
 ,char*  /*a_class */
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  argn = 0;
  args = NULL;

#ifdef G4DEBUG_INTERFACES_COMMON
  printf("G4Wt::G4Wt try to inited Wt\n");
#endif
  // Check if Wt already init in another external app
  if(wApp) {
    WtInited  = TRUE;
    SetMainInteractor (wServer);
    //#endif
    SetArguments      (a_argn,a_args);
    
#ifdef G4DEBUG_INTERFACES_COMMON
    printf("G4Wt::G4Wt alredy inited in external \n");
#endif
  } else {
    
    if(WtInited==FALSE) {  //Wt should be Inited once !
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
#ifdef G4DEBUG_INTERFACES_COMMON
      printf("G4Wt::G4Wt WAppl \n");
#endif
      Wt::WApplication:: instance();
      ////
      //WApplication (*p_argn, args);

      if(!wApp) {
        
        G4cout        << "G4Wt : Unable to init Wt." << G4endl;
      } else {
        WtInited  = TRUE;
        if (a_argn != 0) {
#ifdef G4DEBUG_INTERFACES_COMMON
          printf("G4Wt::G4Wt SetMainInteractor\n");
#endif
          //          SetMainInteractor ();
        }
        SetArguments      (a_argn,a_args);
#ifdef G4DEBUG_INTERFACES_COMMON
        printf("G4Wt::G4Wt inited Wt END\n");
#endif
      }
    }
  }
#ifdef G4DEBUG_INTERFACES_COMMON
  if (wApp) {
    printf("G4Wt::wApp already exist\n");
  }  else {
    printf("G4Wt::wApp not exist\n");
  }
#endif
  //  AddDispatcher     ((G4DispatchFunction)XtDispatchEvent);


  DisableSecondaryLoop ();
}
/***************************************************************************/
G4Wt::~G4Wt (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(this==instance) {
    instance = NULL;
  }
}
/***************************************************************************/
G4bool G4Wt::Inited (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return WtInited;
}
/***************************************************************************/
void* G4Wt::GetEvent (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
//FIXME
//   G4cout        << "G4Wt : Rien compris a cette fonction G4Wt::GetEvent." << G4endl;
//  static XEvent  event;
//  if(appContext==NULL) return NULL;
//  if(mainApp==NULL) return NULL;
//  WtAppNextEvent (appContext, &event);
//  return         &event;
  printf("*");
  return 0;
}
/***************************************************************************/
void G4Wt::FlushAndWaitExecution (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  //  printf("G4Wt::FlushAndWaitExecution ::  Flush ....\n");
  if(!wApp) return;
  //  wApp->processEvents();
}

#endif



