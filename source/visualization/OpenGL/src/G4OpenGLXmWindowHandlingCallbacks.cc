// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmWindowHandlingCallbacks.cc,v 1.2 1999-01-09 16:23:50 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  16th June 1997
// G4OpenGLXmWindowHandlingCallbacks : Callback functions for
//                                     (Motif) widgets to use.
//                                     in handling (Xm) windows

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLXmViewer.hh"

void G4OpenGLXmViewer::resize_callback (Widget w, 
				      XtPointer clientData, 
				      XtPointer) 
{
  Dimension width, height;
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*) clientData;
  
  XtVaGetValues (w, 
		 XmNwidth, &width, 
		 XmNheight, &height, 
		 NULL);
  
  pView->WinSize_x = (G4int) width;
  pView->WinSize_y = (G4int) height;
}



void G4OpenGLXmViewer::expose_callback (Widget w, 
				      XtPointer clientData, 
				      XtPointer) 
{
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*) clientData;
  Dimension width, height;

  XtVaGetValues (w, 
		 XmNwidth, &width, 
		 XmNheight, &height, 
		 NULL);

  pView->WinSize_x = (G4int) width;
  pView->WinSize_y = (G4int) height;

  glXMakeCurrent (pView->dpy, XtWindow(pView->glxarea), pView->cx);
  glViewport (0, 0, width, height);

  pView->ClearView ();
  pView->DrawView ();
}

void G4OpenGLXmViewer::print_callback (Widget, 
				    XtPointer clientData, 
				    XtPointer) 
{
  G4OpenGLXViewer* pView = (G4OpenGLXmViewer*) clientData;
  pView->print();
}

void G4OpenGLXmViewer::set_print_colour_callback (Widget w,
						XtPointer clientData,
						XtPointer) 
{
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*)clientData;
  
  G4int choice = G4OpenGLXmViewer::get_int_userData (w);
  
  pView->print_colour=(G4bool)choice;
  G4cout << "Print colour set to " << pView->print_colour;
  
}

void G4OpenGLXmViewer::set_print_style_callback (Widget w,
					       XtPointer clientData,
					       XtPointer) 
{
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*)clientData;
  
  G4int choice = G4OpenGLXmViewer::get_int_userData (w);
  
  pView->vectored_ps=(G4bool)choice;
  G4cout << "`Produce vectored PostScript ?' set to : " << pView->print_colour;
  
}

#endif
