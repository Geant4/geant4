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
//
// $Id: G4OpenGLXmRadioButton.cc,v 1.5 2001-07-11 10:08:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//Radio button class. Inherits from G4OpenGLXmVWidgetComponent

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLXmVWidgetComponent.hh"
#include "G4OpenGLXmVWidgetContainer.hh"
#include "G4OpenGLXmRadioButton.hh"
#include <X11/Intrinsic.h>
#include "globals.hh"

G4OpenGLXmRadioButton::G4OpenGLXmRadioButton (const char* n,
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

void G4OpenGLXmRadioButton::SetName (const char* n) 
{
  name = n;
  XmString button_string = XmStringCreateLocalized ((char*)name);
  XtVaSetValues (button,
		 XmNlabelString, button_string,
		 NULL);
  XmStringFree (button_string);
}

const char* G4OpenGLXmRadioButton::GetName () 
{
  return name;
}

void G4OpenGLXmRadioButton::AddYourselfTo (G4OpenGLXmVWidgetContainer* container)
{

  pView = container->GetView ();
  ProcesspView ();
  parent = container->GetPointerToWidget ();
  XmString button_string = XmStringCreateLocalized ((char*)name);
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
