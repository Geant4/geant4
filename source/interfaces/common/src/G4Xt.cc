// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Xt.cc,v 1.5 1999-12-15 14:50:47 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G.Barrand

#if defined(G4INTY_BUILD_XT) || defined(G4INTY_USE_XT)

#include <stdlib.h>
#include <string.h>

#include <X11/Intrinsic.h>
#include <X11/Shell.h>

#include "G4ios.hh"

#include "G4Xt.hh"

#define NewString(str)  \
 ((str) != NULL ? (strcpy((char*)malloc((unsigned)strlen(str) + 1), str)) : NULL)

static void XWidgetIconify                 (Widget);
static void XWidgetUniconify               (Widget);
static void XDisplaySetWindowToNormalState (Display*,Window);

G4Xt* G4Xt::instance    = NULL;

static G4bool XtInited  = FALSE;
static int    argn      = 0;
static char** args      = NULL;
static XtAppContext appContext = NULL;
static Widget topWidget = NULL;
/***************************************************************************/
G4Xt* G4Xt::getInstance (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return G4Xt::getInstance (0,NULL,"Geant4");
}
/***************************************************************************/
G4Xt* G4Xt::getInstance (
 int    a_argn
,char** a_args
,char*  a_class
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if (instance==NULL) {
    instance = new G4Xt(a_argn,a_args,a_class);
  }
  return instance;
}
/***************************************************************************/
G4Xt::G4Xt (
 int    a_argn
,char** a_args
,char*  a_class
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(XtInited==FALSE) {  //Xt should be Inited once !
    if(a_argn!=0) {  //Save args.
      args = (char**)malloc(a_argn * sizeof(char*));
      if(args!=NULL) {
	argn = a_argn;
	for(int argi=0;argi<a_argn;argi++) {
	  args[argi] = (char*)NewString (a_args[argi]);
	}
      }
    }
#if XtSpecificationRelease == 4
    Cardinal     argc;
    argc         = (Cardinal)a_argn;
#else
    int          argc;
    argc         = a_argn;
#endif
    Arg          xargs[1];
    XtSetArg     (xargs[0],XtNgeometry,"100x100"); 
    topWidget    = XtAppInitialize (&appContext,a_class,
				    NULL,(Cardinal)0,
				    &argc,a_args,NULL,
				    xargs,1);
    if(topWidget==NULL) {
      G4cout        << "G4Xt : Unable to init Xt." << G4endl;
    }
    // Restore a_args. XtAppInitialize corrupts the given ones !!!
    if( (a_argn!=0) && (args!=NULL)) {
      for(int argi=0;argi<a_argn;argi++) {
	if(args[argi]!=NULL)
	  strcpy(a_args[argi],args[argi]);
	else
	  a_args[argi] = NULL;
      }
    }
    // If topWidget not realized, pbs with Inventor shells.
    XtSetMappedWhenManaged (topWidget,False);
    XtRealizeWidget (topWidget);
    XtInited = TRUE;
  }
  SetArguments      (argn,args);
  SetMainInteractor (topWidget);
  AddDispatcher     ((G4DispatchFunction)XtDispatchEvent);
}
/***************************************************************************/
G4Xt::~G4Xt (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(this==instance) {
    instance = NULL;
  }
}
/***************************************************************************/
G4bool G4Xt::Inited (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return XtInited;
}
/***************************************************************************/
void* G4Xt::GetEvent (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  static XEvent  event;
  if(appContext==NULL) return NULL;
  if(topWidget==NULL) return NULL;
  XtAppNextEvent (appContext, &event);
  return         &event;
}
/***************************************************************************/
void G4Xt::PutStringInResourceDatabase (
 char* a_string 
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(topWidget==NULL)  return;
  if(a_string==NULL)   return;
  Display*             dpy   = XtDisplay(topWidget);
  XrmDatabase          dbres = XrmGetStringDatabase (a_string);
  if(dbres==NULL)      return;
  XrmDatabase          database = XrmGetDatabase (dpy);
  if(database!=NULL)  {
    XrmMergeDatabases  (dbres,&database);
  } else {
    XrmSetDatabase     (dpy,dbres);
  }
}
/***************************************************************************/
void G4Xt::FlushAndWaitExecution (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(topWidget==NULL) return;
  XSync(XtDisplay(topWidget),False);
}
/***************************************************************************/
/******* From OPACS Xx/v3 package ******************************************/
/***************************************************************************/
/***************************************************************************/
void XWidgetIconify (
 Widget This 
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(This==NULL)          return;
  if(!XtIsRealized(This)) return;
  XIconifyWindow(XtDisplay(This),XtWindow(This),XScreenNumberOfScreen(XtScreen(This)));
}
/***************************************************************************/
void XWidgetUniconify (
 Widget This 
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(This==NULL)          return;
  if(!XtIsRealized(This)) return;
  XDisplaySetWindowToNormalState (XtDisplay(This),XtWindow(This));
}
/***************************************************************************/
void XDisplaySetWindowToNormalState (
 Display* This
,Window a_window
)
/***************************************************************************/
/*
  Used to deiconify a window.
*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  XWMHints         wh;
  if(This==NULL)   return;
  if(a_window==0L) return;
  wh.initial_state = NormalState; 
  wh.flags         = StateHint;
  XSetWMHints      (This,a_window,&wh);
  XMapWindow       (This,a_window);
}

#endif



