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
// $Id: G4OpenGLXmTopLevelShell.cc,v 1.4 2001-07-11 10:08:58 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//Top level shell class

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLXmTopLevelShell.hh"
#include "G4OpenGLXmVWidgetContainer.hh"

G4OpenGLXmTopLevelShell::G4OpenGLXmTopLevelShell (G4OpenGLXmViewer* v,
						  char* n) 
{
  pView = v;
  ProcesspView ();
  name = n;
  toplevel = XtVaCreatePopupShell 
    (name,
     topLevelShellWidgetClass,
     top,
     
     XtNiconName, name,
     XtNtitle, name,
     XmNdeleteResponse, XmDO_NOTHING,
     XmNisHomogeneous, False,
     
     XtNvisual, visual, 
     XtNdepth, depth, 
     XtNcolormap, cmap, 
     XtNborderColor, borcol,
     XtNbackground, bgnd,
     NULL);

  frame = XtVaCreateManagedWidget (name,
				   xmFrameWidgetClass,
				   toplevel,
				   
				   XtNvisual, visual,
				   XtNdepth, depth,
				   XtNcolormap, cmap,
				   XtNborderColor, borcol,
				   XtNbackground, bgnd,
				   
				   NULL);
  
  
  
  top_box =  XtVaCreateManagedWidget (name,
				      xmRowColumnWidgetClass,
				      frame,
				      
				      XmNadjustMargin, True,
				      XmNisHomogeneous, False,
				      
				      XtNvisual, visual,
				      XtNdepth, depth,
				      XtNcolormap, cmap,
				      XtNborderColor, borcol,
				      XtNbackground, bgnd,
				      
				      NULL);  

}

G4OpenGLXmTopLevelShell::~G4OpenGLXmTopLevelShell ()
{
  XtDestroyWidget (toplevel);
}

void G4OpenGLXmTopLevelShell::AddChild (G4OpenGLXmVWidgetContainer* container)
{
  container->AddYourselfTo (this);
}

void G4OpenGLXmTopLevelShell::Realize () 
{
  Cardinal num_children;
  XtVaGetValues (toplevel,
		 XmNnumChildren, &num_children,
		 NULL);
//  G4cout << name << " now parents " << num_children << " children." << G4endl;
  XtManageChild (toplevel);
  XtRealizeWidget (toplevel);
  XtPopup (toplevel, XtGrabNonexclusive);
}

Widget* G4OpenGLXmTopLevelShell::GetPointerToWidget ()
{
  return &top_box;
}

char* G4OpenGLXmTopLevelShell::GetName ()
{
  return name;
}

#endif
