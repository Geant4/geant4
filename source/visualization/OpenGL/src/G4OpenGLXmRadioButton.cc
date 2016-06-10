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
// $Id: G4OpenGLXmRadioButton.cc 68043 2013-03-13 14:27:49Z gcosmo $
//
//Radio button class. Inherits from G4OpenGLXmVWidgetComponent

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLXmVWidgetComponent.hh"
#include "G4OpenGLXmVWidgetContainer.hh"
#include "G4OpenGLXmRadioButton.hh"
#include <X11/Intrinsic.h>
#include <Xm/ToggleB.h>

#include "globals.hh"

G4OpenGLXmRadioButton::G4OpenGLXmRadioButton (const char* n,
					      XtCallbackRec* c,
					      G4bool d,
					      G4int num) 
: button(0)
, parent(0)
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
