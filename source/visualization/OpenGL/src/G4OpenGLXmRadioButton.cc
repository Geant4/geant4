// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmRadioButton.cc,v 1.2 1999-01-09 16:23:40 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//Radio button class. Inherits from G4OpenGLXmVWidgetComponent

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLXmVWidgetComponent.hh"
#include "G4OpenGLXmVWidgetContainer.hh"
#include "G4OpenGLXmRadioButton.hh"
#include <X11/Intrinsic.h>
#include "globals.hh"

G4OpenGLXmRadioButton::G4OpenGLXmRadioButton (char* n,
					      XtCallbackRec* c,
					      G4bool d,
					      G4int num) 
{
  name = n;
  callback = c;
  default_button = d;
  number = num;
}

G4OpenGLXmRadioButton::~G4OpenGLXmRadioButton ()
{}

void G4OpenGLXmRadioButton::SetName (char* n) 
{
  name = n;
  XmString button_string = XmStringCreateLocalized (name);
  XtVaSetValues (button,
		 XmNlabelString, button_string,
		 NULL);
  XmStringFree (button_string);
}

char* G4OpenGLXmRadioButton::GetName () 
{
  return name;
}

void G4OpenGLXmRadioButton::AddYourselfTo (G4OpenGLXmVWidgetContainer* container)
{

  pView = container->GetView ();
  ProcesspView ();
  parent = container->GetPointerToWidget ();
  XmString button_string = XmStringCreateLocalized (name);
  button = XtVaCreateManagedWidget (name,
				    xmToggleButtonWidgetClass,
				    *parent,
				    
				    XmNlabelString, button_string,
				    XmNset, default_button,
				    XmNuserData, number,
				    
				    XtNvisual, visual,
				    XtNdepth, depth,
				    XtNcolormap, cmap,
				    XtNborderColor, borcol,
				    XtNbackground, bgnd,
				    
				    NULL);
  
  XtAddCallbacks (button,
		  XmNarmCallback,
		  callback);

  XmStringFree (button_string);
}

Widget* G4OpenGLXmRadioButton::GetPointerToParent ()
{
  return parent;
}

Widget* G4OpenGLXmRadioButton::GetPointerToWidget () 
{
  return &button;
}

#endif
