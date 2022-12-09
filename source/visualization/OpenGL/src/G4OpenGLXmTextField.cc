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
//
//Text field class. Inherits from G4OpenGLXmVWidgetComponent

#include "G4OpenGLXmViewer.hh"
#include "G4OpenGLXmVWidgetComponent.hh"
#include "G4OpenGLXmVWidgetContainer.hh"
#include "G4OpenGLXmTextField.hh"

#include <X11/Intrinsic.h>
#include <Xm/Label.h>
#include <Xm/TextF.h>

#include "globals.hh"

G4OpenGLXmTextField::G4OpenGLXmTextField (const char* n,
					  G4double* val)
: text_label(0)
, text_field(0)
, parent(0)
{
  name = n;
  initial = new char[50];
  snprintf (initial, 50, "%6.2f", *val);
  value = (void*)val;
  text=false;
}

G4OpenGLXmTextField::G4OpenGLXmTextField (const char* n,
					  const char* val)
: text_label(0)
, text_field(0)
, parent(0)
{
  name = n;
  initial = new char[50];
  snprintf (initial, 50, "%s", val);
  value = (void*)val;
  text=true;
  //  strcpy (initial, val);
}

G4OpenGLXmTextField::~G4OpenGLXmTextField ()
{
  delete[] initial;
}

void G4OpenGLXmTextField::SetName (const char* n) 
{
  name = n;
  XmString text_string = XmStringCreateLocalized ((char*)name);
  XtVaSetValues (text_label,
		 XmNlabelString, text_string,
		 NULL);
  XmStringFree (text_string);
}

const char* G4OpenGLXmTextField::GetName () 
{
  return name;
}

void G4OpenGLXmTextField::SetValue (G4double val)
{
  snprintf (initial, 50, "%6.2f", val);
  
  XtVaSetValues (text_field,
		 XmNvalue, (String)initial,
		 NULL);
  
}

void G4OpenGLXmTextField::SetValue (const char* val)
{
  snprintf (initial, 50, "%s", val);
  //  strcpy (initial, val);

  XtVaSetValues (text_field,
		 XmNvalue, (String)initial,
		 NULL);
  
}

const char* G4OpenGLXmTextField::GetValue ()
{
  return initial;
}

void G4OpenGLXmTextField::AddYourselfTo (G4OpenGLXmVWidgetContainer* container)
{

  pView = container->GetView ();
  ProcesspView ();
  parent = container->GetPointerToWidget ();

  char local_w_text[50];
  strcpy (local_w_text, name);

  char label_name[50];
  strcpy (label_name, name);
  strcat (label_name, "_label");
  
  char text_field_name[50];
  strcpy (text_field_name, name);
  strcat (text_field_name, "_text_field");
  
  XmString local_text = XmStringCreateLocalized (local_w_text);
  text_label = XtVaCreateManagedWidget (label_name, 
					xmLabelWidgetClass,
					*parent,

					XmNlabelString, local_text,
					
					XtNvisual, visual, 
					XtNdepth, depth, 
					XtNcolormap, cmap, 
					XtNborderColor, borcol,
					XtNbackground, bgnd,
					
					NULL);
  XmStringFree (local_text);

  text_field = XtVaCreateManagedWidget (text_field_name,
					xmTextFieldWidgetClass,
					*parent,

					XmNvalue, (String)initial, 
					
					XtNvisual, visual, 
					XtNdepth, depth, 
					XtNcolormap, cmap, 
					XtNborderColor, borcol,
					XtNbackground, bgnd,
					
					NULL);

  if (!text) {
    XtAddCallback (text_field, 
		   XmNvalueChangedCallback,
		   G4OpenGLXmViewer::get_double_value_callback,
		   value);
  } else {
    XtAddCallback (text_field, 
		   XmNvalueChangedCallback,
		   G4OpenGLXmViewer::get_text_callback,
		   value);
  }
}

Widget* G4OpenGLXmTextField::GetPointerToParent ()
{
  return parent;
}

Widget* G4OpenGLXmTextField::GetPointerToWidget () 
{
  return &text_field;
}
