// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmStyleCallbacks.cc,v 1.3 1999-12-15 14:54:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  16th April 1997
// G4OpenGLXmStyleCallbacks : 
//                       Several callback functions used by
//                       elements of the control panel to
//                       determine how to visualize the view.

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLXmViewer.hh"

void G4OpenGLXmViewer::drawing_style_callback (Widget w, 
					     XtPointer clientData, 
					     XtPointer) 
{
  G4long choice = (G4long)clientData;
  G4OpenGLXmViewer* pView;
  XtVaGetValues (XtParent(w),
		 XmNuserData, &pView,
		 NULL);
  G4ViewParameters::DrawingStyle style;

  switch (choice) {
    
  case 0:
    style = G4ViewParameters::wireframe;
    break;

  case 1:
    style = G4ViewParameters::hlr;
    break;

  case 2:
    style = G4ViewParameters::hsr;
    break;

  case 3:
    style = G4ViewParameters::hlhsr;
    break;

  default:
    G4Exception("Unrecognised case in drawing_style_callback.");
  }

  pView->fVP.SetDrawingStyle (style);
  
  pView->SetView ();
  pView->ClearView ();
  pView->DrawView ();
}

void G4OpenGLXmViewer::rep_style_callback (Widget w, 
					 XtPointer clientData, 
					 XtPointer) 
{
  G4long choice = (G4long)clientData;
  G4OpenGLXmViewer* pView;
  XtVaGetValues (XtParent(w),
		 XmNuserData, &pView,
		 NULL);
  G4ViewParameters::RepStyle style;

  switch (choice) {
    
  case 0:
    style = G4ViewParameters::polyhedron;
    break;

  case 1:
    style = G4ViewParameters::nurbs;
    break;

  default:
    G4Exception("Unrecognised case in rep_style_callback.");
  }

  pView->fVP.SetRepStyle (style);

  pView->SetView ();
  pView->ClearView ();
  pView->DrawView ();
}

void G4OpenGLXmViewer::background_color_callback (Widget w, 
						XtPointer clientData, 
						XtPointer) 
{
  G4long choice = (G4long)clientData;
  G4OpenGLXmViewer* pView;
  XtVaGetValues (XtParent(w),
		 XmNuserData, &pView,
		 NULL);


  //I need to revisit the kernel if the background colour changes and hidden
  //line removal is enabled, because hlr drawing utilises the background
  //colour in its drawing...
  switch (choice) {
    
  case 0:
    if (!pView->white_background) {
      pView->white_background = true;
      if (pView->GetViewParameters().GetDrawingStyle() == G4ViewParameters::hlr) {
	pView->SetNeedKernelVisit ();
      }
    }
    break;

  case 1:
    if (pView->white_background) {
      pView->white_background = false;
      if (pView->GetViewParameters().GetDrawingStyle() == G4ViewParameters::hlr) {
	pView->SetNeedKernelVisit ();
      }
    }
    break;

  default:
    G4Exception("Unrecognised case in background_color_callback.");
  }

  pView->SetView ();
  pView->ClearView ();
  pView->DrawView ();
}

void G4OpenGLXmViewer::transparency_callback (Widget w, 
					    XtPointer clientData, 
					    XtPointer) 
{
  G4long choice = (G4long)clientData;
  G4OpenGLXmViewer* pView;
  XtVaGetValues (XtParent(w),
		 XmNuserData, &pView,
		 NULL);

  switch (choice) {
    
  case 0:
    pView->transparency_enabled = false;
    glDisable (GL_BLEND);
    break;

  case 1:
    pView->transparency_enabled = true;
    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glShadeModel (GL_FLAT);
    break;

  default:
    G4Exception("Unrecognised case in transparency_callback.");
  }

  pView->SetView ();
  pView->ClearView ();
  pView->DrawView ();
}

void G4OpenGLXmViewer::antialias_callback (Widget w, 
					 XtPointer clientData, 
					 XtPointer) 
{
  G4long choice = (G4long)clientData;
  G4OpenGLXmViewer* pView;
  XtVaGetValues (XtParent(w),
		 XmNuserData, &pView,
		 NULL);

  switch (choice) {
    
  case 0:
    pView->antialiasing_enabled = false;
    glDisable (GL_LINE_SMOOTH);
    glDisable (GL_POLYGON_SMOOTH);
    break;

  case 1:
    pView->antialiasing_enabled = true;
    glEnable (GL_LINE_SMOOTH);
    glHint (GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable (GL_POLYGON_SMOOTH);
    glHint (GL_POLYGON_SMOOTH_HINT, GL_NICEST);
    break;

  default:
    G4Exception("Unrecognised case in antialiasing_callback.");
  }

  pView->SetView ();
  pView->ClearView ();
  pView->DrawView ();
}

void G4OpenGLXmViewer::haloing_callback (Widget w, 
				       XtPointer clientData, 
				       XtPointer) 
{
  G4long choice = (G4long)clientData;
  G4OpenGLXmViewer* pView;
  XtVaGetValues (XtParent(w),
		 XmNuserData, &pView,
		 NULL);

  switch (choice) {
    
  case 0:
    pView->haloing_enabled = false;
    break;

  case 1:
    pView->haloing_enabled = true;
    break;

  default:
    G4Exception("Unrecognised case in haloing_callback.");
  }

  pView->SetView ();
  pView->ClearView ();
  pView->DrawView ();
}

void G4OpenGLXmViewer::projection_callback (Widget w, 
					  XtPointer clientData, 
					  XtPointer) 
{
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*)clientData;

  G4int choice = G4OpenGLXmViewer::get_int_userData (w);

  switch (choice) {
  case 0:
    {
      pView->fVP.SetFieldHalfAngle (0.);
      break;
    }

  case 1:
    {
      if (pView->fov > 89.5 || pView->fov <= 0.0) {
	G4cout << "Field half angle should be 0 < angle <= 89.5 degrees.";
	G4cout << G4endl;
      }
      else {
	pView->fVP.SetFieldHalfAngle (pView->fov * deg);
      }
      break;
    }
  default:
    {
      G4Exception("Unrecognised choice made in projection_callback");
    }
  }

  pView->SetView ();
  pView->ClearView ();
  pView->DrawView ();
}  

#endif

