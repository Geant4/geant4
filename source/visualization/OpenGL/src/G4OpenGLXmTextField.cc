// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmTextField.cc,v 1.2 1999-01-09 16:23:44 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//Text field class. Inherits from G4OpenGLXmVWidgetComponent

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLXmVWidgetComponent.hh"
#include "G4OpenGLXmVWidgetContainer.hh"
#include "G4OpenGLXmTextField.hh"
#include <X11/Intrinsic.h>
#include "globals.hh"

G4OpenGLXmTextField::G4OpenGLXmTextField (char* n,
					  G4double* val)
{
  name = n;
  initial = new char[50];
  sprintf (initial, "%6.2f\0", *val);
  value = (void*)val;
  text=false;
}

G4OpenGLXmTextField::G4OpenGLXmTextField (char* n,
					  char* val)
{
  name = n;
  initial = new char[50];
  sprintf (initial, "%s", val);
  value = (void*)val;
  text=true;
  //  strcpy (initial, val);
}

G4OpenGLXmTextField::~G4OpenGLXmTextField ()
{
  delete[] initial;
}

void G4OpenGLXmTextField::SetName (char* n) 
{
  name = n;
  XmString text_string = XmStringCreateLocalized (name);
  XtVaSetValues (text_label,
		 XmNlabelString, text_string,
		 NULL);
  XmStringFree (text_string);
}

char* G4OpenGLXmTextField::GetName () 
{
  return name;
}

void G4OpenGLXmTextField::SetValue (G4double val)
{
  sprintf (initial, "%6.2f\0", val);
  
  XtVaSetValues (text_field,
		 XmNvalue, (String)initial,
		 NULL);
  
}

void G4OpenGLXmTextField::SetValue (char* val)
{
  sprintf (initial, "%s", val);
  //  strcpy (initial, val);

  XtVaSetValues (text_field,
		 XmNvalue, (String)initial,
		 NULL);
  
}

char* G4OpenGLXmTextField::GetValue ()
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

#endif
