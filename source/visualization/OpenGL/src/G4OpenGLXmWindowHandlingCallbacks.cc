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
// $Id: G4OpenGLXmWindowHandlingCallbacks.cc 98766 2016-08-09 14:17:17Z gcosmo $
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
				      XtPointer x) 
{
  expose_callback(w,clientData,x);
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

  pView->ResizeWindow(width,height);

//??????????????????????????? This might be a problem in MT mode.
//  glXMakeCurrent (pView->dpy, XtWindow(pView->glxarea), pView->cxMaster);
//  pView->SetView ();
//  pView->ClearView ();
//  pView->DrawView ();
//???????????????????????????? Commented out 14/06/16  JA
}

void G4OpenGLXmViewer::print_callback (Widget, 
				    XtPointer clientData, 
				    XtPointer) 
{
  G4OpenGLXViewer* pView = (G4OpenGLXmViewer*) clientData;
  pView->printEPS();
}

void G4OpenGLXmViewer::set_print_colour_callback (Widget w,
						XtPointer clientData,
						XtPointer) 
{
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*)clientData;
  
  G4int choice = get_int_userData (w);
  
  pView->fPrintColour=(G4bool)choice;
  G4cout << "Print colour set to " << pView->fPrintColour;
  
}

void G4OpenGLXmViewer::set_print_style_callback (Widget w,
					       XtPointer clientData,
					       XtPointer) 
{
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*)clientData;
  
  G4int choice = get_int_userData (w);
  
  pView->fVectoredPs=(G4bool)choice;
  G4cout << "`Produce vectored PostScript ?' set to : " << pView->fPrintColour;
  
}

#endif
