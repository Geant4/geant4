// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmPushButton.cc,v 1.3 1999-12-15 14:54:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//Push button class. Inherits from G4OpenGLXmVWidgetComponent

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLXmVWidgetComponent.hh"
#include "G4OpenGLXmVWidgetContainer.hh"
#include "G4OpenGLXmPushButton.hh"
#include <X11/Intrinsic.h>
#include "globals.hh"

G4OpenGLXmPushButton::G4OpenGLXmPushButton (char* n,
					    XtCallbackRec* c) 
{
  name = n;
  callback = c;
}

G4OpenGLXmPushButton::~G4OpenGLXmPushButton ()
{}

void G4OpenGLXmPushButton::SetName (char* n) 
{
  name = n;
  XmString button_string = XmStringCreateLocalized (name);
  XtVaSetValues (button,
		 XmNlabelString, button_string,
		 NULL);
  XmStringFree (button_string);
}

char* G4OpenGLXmPushButton::GetName () 
{
  return name;
}

void G4OpenGLXmPushButton::AddYourselfTo (G4OpenGLXmVWidgetContainer* container)
{

  pView = container->GetView ();
  ProcesspView ();
  parent = container->GetPointerToWidget ();

  XmString button_str = XmStringCreateLocalized (name);
  button = XtVaCreateManagedWidget 
    (name,
     xmPushButtonWidgetClass,
     *parent,
     XmNlabelString, button_str,
     XmNalignment, XmALIGNMENT_CENTER,
     XmNuserData, pView,

     XtNvisual, visual, 
     XtNdepth, depth, 
     XtNcolormap, cmap, 
     XtNborderColor, borcol,
     XtNbackground, bgnd,
     
     NULL);
  
  XtAddCallbacks (button,
		  XmNarmCallback,
		  callback);
  
  XmStringFree (button_str);
}

Widget* G4OpenGLXmPushButton::GetPointerToParent ()
{
  return parent;
}

Widget* G4OpenGLXmPushButton::GetPointerToWidget () 
{
  return &button;
}

#endif


