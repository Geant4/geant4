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
// $Id: G4OpenGLXmPushButton.cc 68043 2013-03-13 14:27:49Z gcosmo $
//
//Push button class. Inherits from G4OpenGLXmVWidgetComponent

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLXmVWidgetComponent.hh"
#include "G4OpenGLXmVWidgetContainer.hh"
#include "G4OpenGLXmPushButton.hh"
#include <X11/Intrinsic.h>
#include <Xm/PushB.h>

#include "globals.hh"

G4OpenGLXmPushButton::G4OpenGLXmPushButton (const char* n,
					    XtCallbackRec* c)
: button(0)
, parent(0)
{
  name = n;
  callback = c;
}

G4OpenGLXmPushButton::~G4OpenGLXmPushButton ()
{}

void G4OpenGLXmPushButton::SetName (const char* n) 
{
  name = n;
  XmString button_string = XmStringCreateLocalized ((char*)name);
  XtVaSetValues (button,
		 XmNlabelString, button_string,
		 NULL);
  XmStringFree (button_string);
}

const char* G4OpenGLXmPushButton::GetName () 
{
  return name;
}

void G4OpenGLXmPushButton::AddYourselfTo (G4OpenGLXmVWidgetContainer* container)
{

  pView = container->GetView ();
  ProcesspView ();
  parent = container->GetPointerToWidget ();

  XmString button_str = XmStringCreateLocalized ((char*)name);
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


