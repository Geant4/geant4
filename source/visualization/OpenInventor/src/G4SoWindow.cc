//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define SbWarn printf

#include <HEPVis/SbBasic.h>

// For backward compatibility :
#ifdef HEPVisWin
#define HEPVisWin32
#endif

#ifdef CDF
#ifndef LINUX
#define u_long unsigned long
#endif
#endif
#ifdef HEPVisXt
#include <X11/Intrinsic.h>
#include <X11/Shell.h>
#include <X11/StringDefs.h>
#include <X11/IntrinsicP.h>
#endif
#ifdef HEPVisWin32
#include <windowsx.h> //For GetWindowInstance.
#endif

#include <Inventor/SbString.h>
#include <Inventor/SbDict.h>

#include "G4SoWindow.hh"

#ifdef HEPVisXt
static XtTranslations trans_table = NULL;
static XtActionsRec closeAction = 
{(String)"CloseG4SoWindow",G4SoWindow::closeWindow};
//static void XWidgetSetDimension(Widget,Dimension,Dimension);
#endif
#ifdef HEPVisWin32
static LRESULT CALLBACK WindowProc (HWND,UINT,WPARAM,LPARAM);
static const char className[] = "G4SoWindow";
#endif

SbDict* G4SoWindow::fWidgetDictionary = NULL;

//////////////////////////////////////////////////////////////////////////////
void G4SoWindow::clearClass (
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fWidgetDictionary!=NULL) delete fWidgetDictionary;
  fWidgetDictionary = NULL;
}
//////////////////////////////////////////////////////////////////////////////
G4SoWindow::G4SoWindow(
 const char* aName
,SbBool aWMClose
)
:fWidget(NULL)
,fWMClose(aWMClose)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
#ifdef HEPVisXt
  if(trans_table==NULL) {
    trans_table = 
      XtParseTranslationTable("<Message>WM_PROTOCOLS:CloseG4SoWindow()");
    XtAppAddActions(SoXt::getAppContext(),&closeAction,1);
  }
  SbString shellName = aName;
  shellName += "_shell"; 
  Arg args[2];
  XtSetArg(args[0],XtNgeometry,XtNewString("400x400"));
  XtSetArg(args[1],XtNborderWidth,0);
  fWidget = XtAppCreateShell(shellName.getString(),"Inventor",
				  topLevelShellWidgetClass,
				  SoXt::getDisplay(),
				  args,2); 
  if(fWidget!=NULL) {
    XtOverrideTranslations(fWidget,trans_table);
    XtAddCallback(fWidget,XtNdestroyCallback,
		  (XtCallbackProc)G4SoWindow::widgetDestroyedCB,NULL);
  }
#endif
#ifdef HEPVisWin32
  static SbBool done = FALSE;
  if(done==FALSE) {
    WNDCLASS wc;
    wc.style = CS_HREDRAW | CS_VREDRAW;
    wc.lpfnWndProc = (WNDPROC)WindowProc;
    wc.cbClsExtra = 0;
    wc.cbWndExtra = 0;
    wc.hInstance = ::GetModuleHandle(NULL);
    wc.hIcon = ::LoadIcon(NULL, IDI_APPLICATION);
    wc.hCursor = ::LoadCursor(NULL, IDC_ARROW);
    wc.hbrBackground = NULL;
    wc.lpszMenuName = className;
    wc.lpszClassName = className;
    ::RegisterClass(&wc);
    done = TRUE;
  }
  //  Compell window to be created at 0,0 to bypass 
  // the 'black border' problem. 
  fWidget = ::CreateWindow(className, aName, 
			   WS_OVERLAPPEDWINDOW,
			   //CW_USEDEFAULT, CW_USEDEFAULT, 
			   0,0, 
			   400,400, 
			   NULL, NULL,
			   ::GetModuleHandle(NULL),
			   NULL);
  // Retreive window and client sizez :
  RECT wrect,crect;
  GetWindowRect((HWND)fWidget,&wrect);
  GetClientRect((HWND)fWidget,&crect);
  int ww = wrect.right-wrect.left;
  int wh = wrect.bottom-wrect.top;
  int cw = crect.right-crect.left;
  int ch = crect.bottom-crect.top;
  // Compell client rect to be 400 400 :
  MoveWindow((HWND)fWidget,wrect.left,wrect.top,400+ww-cw,400+wh-ch,TRUE);
  //
  ::SetWindowLong((HWND)fWidget,GWL_USERDATA,LONG(this));
#endif
  registerWidget(fWidget);
}
//////////////////////////////////////////////////////////////////////////////
G4SoWindow::~G4SoWindow(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fWidget==NULL) return;
  unregisterWidget(fWidget);
#ifdef HEPVisXt
  XtDestroyWidget(fWidget);
#endif
#ifdef HEPVisWin32
  ::SetWindowLong((HWND)fWidget,GWL_USERDATA,LONG(NULL));
  ::DestroyWindow((HWND)fWidget);
#endif
}
//////////////////////////////////////////////////////////////////////////////
Widget G4SoWindow::getWidget(
) const
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{ 
  return fWidget; 
}
//////////////////////////////////////////////////////////////////////////////
void G4SoWindow::setWMClose(
 SbBool aValue 			  
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{ 
  fWMClose = aValue;
}
//////////////////////////////////////////////////////////////////////////////
SbBool G4SoWindow::getWMClose(
) const
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{ 
  return fWMClose;
}
//////////////////////////////////////////////////////////////////////////////
void G4SoWindow::show(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fWidget==NULL) return;
  SoXt::show(fWidget);
}
//////////////////////////////////////////////////////////////////////////////
void G4SoWindow::hide(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fWidget==NULL) return;
  SoXt::hide(fWidget);
}
//////////////////////////////////////////////////////////////////////////////
void G4SoWindow::setTitle(
 const char* aString
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(aString==NULL) return;
#ifdef HEPVisXt
  Arg args[1];
  XtSetArg(args[0],XtNtitle,XtNewString(aString));
  XtSetValues(fWidget,args,1);
#endif
#ifdef HEPVisWin32
  ::SetWindowText((HWND)fWidget,aString);
#endif
}
//////////////////////////////////////////////////////////////////////////////
void G4SoWindow::setSize(
 const SbVec2s& aSize
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fWidget==NULL) return;
#ifdef HEPVisXt
  Arg args[1];
  char s[32];
  sprintf(s,"%dx%d",aSize[0],aSize[1]);
  XtSetArg(args[0],XtNgeometry,XtNewString(s));
  XtSetValues(fWidget,args,1);
  //XWidgetSetDimension(fWidget,aSize[0],aSize[1]);
#endif
#ifdef HEPVisWin32
  // Retreive window and client sizez :
  RECT wrect,crect;
  GetWindowRect((HWND)fWidget,&wrect);
  GetClientRect((HWND)fWidget,&crect);
  int ww = wrect.right-wrect.left;
  int wh = wrect.bottom-wrect.top;
  int cw = crect.right-crect.left;
  int ch = crect.bottom-crect.top;
  // Compell client rect to be 400 400 :
  MoveWindow((HWND)fWidget,
	       wrect.left,wrect.top,
	       aSize[0]+ww-cw,aSize[1]+wh-ch,TRUE);
#endif
}
//////////////////////////////////////////////////////////////////////////////
SbVec2s G4SoWindow::getSize(
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(fWidget==NULL) return SbVec2s();
#ifdef HEPVisXt
  Dimension width,height;
  Arg args[2];
  XtSetArg(args[0],XtNwidth  ,&width);
  XtSetArg(args[1],XtNheight,&height);
  XtGetValues(fWidget,args,2);
  return SbVec2s(width,height);
#endif
#ifdef HEPVisWin32
  RECT rect;
  ::GetClientRect((HWND)fWidget,&rect);
  return SbVec2s((short)(rect.right-rect.left),(short)(rect.bottom-rect.top));
#endif
}
//////////////////////////////////////////////////////////////////////////////
void G4SoWindow::registerWidget(
 Widget aWidget
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(aWidget==NULL) return;
  if(fWidgetDictionary==NULL) fWidgetDictionary = new SbDict;
  fWidgetDictionary->enter((unsigned long)aWidget,this);
}
//////////////////////////////////////////////////////////////////////////////
void G4SoWindow::unregisterWidget(
 Widget aWidget
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(aWidget==NULL) return;
  if(fWidgetDictionary==NULL) return;
  fWidgetDictionary->remove((unsigned long)aWidget);
}
//////////////////////////////////////////////////////////////////////////////
G4SoWindow* G4SoWindow::getWindow(
 Widget aWidget
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(aWidget==NULL) return NULL;
  if(fWidgetDictionary==NULL) return NULL;
  void*  value;
  SbBool found = fWidgetDictionary->find((unsigned long)aWidget,value);
  return (G4SoWindow*)(found==FALSE ? NULL : value);
}
#ifdef HEPVisXt
//////////////////////////////////////////////////////////////////////////////
void G4SoWindow::widgetDestroyedCB(
 Widget aWidget
,XtPointer
,XtPointer
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  G4SoWindow* This = G4SoWindow::getWindow(aWidget);
  if(This==NULL) return; // XtDestroyWidget called from ~G4SoWindow.
  //printf("debug : G4SoWindow::widgetDestroyCB\n");
  This->unregisterWidget(aWidget);
}
//////////////////////////////////////////////////////////////////////////////
void G4SoWindow::closeWindow(
 Widget aWidget
,XEvent*   /*aEvent*/
,char**    /*aArgs*/
,Cardinal* /*aNarg*/
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  G4SoWindow* This = G4SoWindow::getWindow(aWidget);
  if( (This!=NULL) && (This->getWMClose()==TRUE) ) {
    //SbWarn ("debug : CloseWindow : close\n");
#ifdef SoFreePackage
    SoXt::sendExit ();
#else
    exit (EXIT_SUCCESS);
#endif
  } else {
    // Should destroy G4SoWindow.
    XtUnmapWidget(aWidget);
  }
}
/*
#include <X11/IntrinsicP.h> //Dangerous.
void XWidgetSetDimension(Widget aWidget,Dimension a_width,Dimension a_height){
  if(aWidget==NULL)      return;
  if((a_width<=0)||(a_height<=0)) return;
  XtResizeWidget (aWidget,a_width,a_height,aWidget->core.border_width);
}
*/
#endif
#ifdef HEPVisWin32
//////////////////////////////////////////////////////////////////////////////
LRESULT CALLBACK WindowProc ( 
 HWND   aWindow
,UINT   aMessage
,WPARAM aWParam
,LPARAM aLParam
)
//////////////////////////////////////////////////////////////////////////////
//  Below treatment of WM_SIZE, WM_SETFOCUS not necessary 
// with TGS, but needed with SoFree. WM_DESTROY needed for 
// 'main top level window' so that 'Close window' induces
// the end of the task.
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  switch (aMessage) { 
  case WM_SIZE:{ // Assume one child window !
    int width = LOWORD(aLParam);
    int height = HIWORD(aLParam);
    //printf("debug : G4SoWindow : WMS_SIZE : %d %d\n",width,height);
    HWND hwnd = ::GetFirstChild(aWindow);
    if(hwnd!=NULL) {
      ::MoveWindow(hwnd,0,0,width,height,TRUE);
    }
  }return 0;
  case WM_SETFOCUS:{ // Assume one child window !
    HWND hwnd = ::GetFirstChild(aWindow);
    if(hwnd!=NULL) ::SetFocus(hwnd);
  }return 0;
  case WM_DESTROY:{
    G4SoWindow* This = (G4SoWindow*)::GetWindowLong(aWindow,GWL_USERDATA);
    if( (This!=NULL) && (This->getWMClose()==TRUE) ) {
      ::PostQuitMessage(0);
    }
  }return 0;
  default:
    return (::DefWindowProc(aWindow,aMessage,aWParam,aLParam));
  }
}
#endif




