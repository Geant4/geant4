// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmTopLevelShell.cc,v 1.3 1999-12-15 14:54:10 gunter Exp $
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
