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
// $Id: G4OpenGLXmBox.cc,v 1.5 2001-07-11 10:08:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//Box container class

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLXmFramedBox.hh"
#include "G4OpenGLXmBox.hh"
#include "G4OpenGLXmVWidgetComponent.hh"
#include "G4OpenGLXmVWidgetShell.hh"

G4OpenGLXmBox::G4OpenGLXmBox (const char* n, 
			      G4bool r)
{
  name = n;
  radio = r;
  parent = NULL;
  box_row_col = NULL;
}

G4OpenGLXmBox::~G4OpenGLXmBox () 
{}

void G4OpenGLXmBox::AddChild (G4OpenGLXmVWidgetComponent* component)
{
  component->AddYourselfTo(this);
  Cardinal num_children;
  XtVaGetValues (box_row_col,
		 XmNnumChildren, &num_children,
		 NULL);
//  G4cout << name << " now parents " << num_children << " children." << G4endl;
}

void G4OpenGLXmBox::AddYourselfTo (G4OpenGLXmVWidgetShell* window)
{

  pView = window->GetView ();
  ProcesspView ();
  parent = window->GetPointerToWidget ();
  
  box_row_col =  XtVaCreateManagedWidget (name,
					  xmRowColumnWidgetClass,
					  *parent,
					  
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

Widget* G4OpenGLXmBox::GetPointerToParent ()
{
  return parent;
}

Widget* G4OpenGLXmBox::GetPointerToWidget ()
{
  return &box_row_col;
}

const char* G4OpenGLXmBox::GetName ()
{
  return name;
}

void G4OpenGLXmBox::SetName (const char* n)
{
  name = n;
}

#endif
