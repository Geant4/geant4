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
// $Id: G4OpenGLXmSeparator.cc,v 1.4 2001-07-11 10:08:58 gunter Exp $
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
