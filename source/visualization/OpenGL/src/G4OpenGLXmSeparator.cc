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
// $Id: G4OpenGLXmSeparator.cc 68043 2013-03-13 14:27:49Z gcosmo $
//
//Separator class. Inherits from G4OpenGLXmVWidgetComponent

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLXmVWidgetComponent.hh"
#include "G4OpenGLXmVWidgetContainer.hh"
#include "G4OpenGLXmSeparator.hh"
#include <X11/Intrinsic.h>
#include <Xm/Separator.h>

#include "globals.hh"

G4OpenGLXmSeparator::G4OpenGLXmSeparator (unsigned char l) 
: line(0)
, parent(0)
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
