// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmFramedBox.cc,v 1.3 1999-12-15 14:54:09 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//Framed box container class

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLXmBox.hh"
#include "G4OpenGLXmFramedBox.hh"
#include "G4OpenGLXmVWidgetComponent.hh"
#include "G4OpenGLXmVWidgetShell.hh"

G4OpenGLXmFramedBox::G4OpenGLXmFramedBox (char* n, 
					  G4bool r) :
G4OpenGLXmBox (n, r)
{
  frame = NULL;
}

G4OpenGLXmFramedBox::~G4OpenGLXmFramedBox () 
{}

void G4OpenGLXmFramedBox::AddChild (G4OpenGLXmVWidgetComponent* component)
{
  component->AddYourselfTo(this);
  Cardinal num_children;
  XtVaGetValues (box_row_col,
		 XmNnumChildren, &num_children,
		 NULL);
//  G4cout << name << " now parents " << num_children << " children." << G4endl;
}

void G4OpenGLXmFramedBox::AddYourselfTo (G4OpenGLXmVWidgetShell* window)
{

  pView = window->GetView ();
  ProcesspView ();
  char framename[50];
  strcpy (framename, name);
  strcat (framename, "_frame");

  parent = window->GetPointerToWidget ();
  frame = XtVaCreateManagedWidget (framename,
				    xmFrameWidgetClass,
				    *parent,
				    
				    XtNvisual, visual,
				    XtNdepth, depth,
				    XtNcolormap, cmap,
				    XtNborderColor, borcol,
				    XtNbackground, bgnd,
				    
				    NULL);
  
  
  
  box_row_col =  XtVaCreateManagedWidget (name,
					  xmRowColumnWidgetClass,
					  frame,
					  
					  XmNadjustMargin, True,
					  XmNisHomogeneous, False,
					  XmNlabelString, (XmString)name,
					  XmNradioAlwaysOne, radio,
					  XmNradioBehavior, radio,
					  
					  XtNvisual, visual,
					  XtNdepth, depth,
					  XtNcolormap, cmap,
					  XtNborderColor, borcol,
					  XtNbackground, bgnd,
					  
					  NULL);
  
}

#endif
