// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmSeparator.cc,v 1.1 1999-01-07 16:15:02 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//Separator class. Inherits from G4OpenGLXmVWidgetComponent

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLXmVWidgetComponent.hh"
#include "G4OpenGLXmVWidgetContainer.hh"
#include "G4OpenGLXmSeparator.hh"
#include <X11/Intrinsic.h>
#include "globals.hh"

G4OpenGLXmSeparator::G4OpenGLXmSeparator (unsigned char l) 
{
  line_type = l;
}

G4OpenGLXmSeparator::~G4OpenGLXmSeparator ()
{}

void G4OpenGLXmSeparator::AddYourselfTo (G4OpenGLXmVWidgetContainer* container)
{

  pView = container->GetView ();
  ProcesspView ();

  parent = container->GetPointerToWidget ();

  line = XtVaCreateManagedWidget ("sep",
				  xmSeparatorWidgetClass,
				  *parent,
				  XmNseparatorType, line_type,
				  XmNmargin, 1,
				  XmNorientation, XmHORIZONTAL,
				  
				  XtNvisual, visual, 
				  XtNdepth, depth, 
				  XtNcolormap, cmap, 
				  XtNborderColor, borcol,
				  XtNbackground, bgnd,
				  
				  NULL);
}

Widget* G4OpenGLXmSeparator::GetPointerToParent ()
{
  return parent;
}

Widget* G4OpenGLXmSeparator::GetPointerToWidget () 
{
  return &line;
}

#endif
