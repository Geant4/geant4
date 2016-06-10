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
// $Id: G4OpenGLXmTopLevelShell.cc 66373 2012-12-18 09:41:34Z gcosmo $
//
//Top level shell class

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLXmViewer.hh"
#include "G4OpenGLXmTopLevelShell.hh"
#include "G4OpenGLXmVWidgetContainer.hh"

#include <Xm/Frame.h>
#include <Xm/RowColumn.h>

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
