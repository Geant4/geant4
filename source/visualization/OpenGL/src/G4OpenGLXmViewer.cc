// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmViewer.cc,v 1.6 2001-03-07 15:29:50 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  10th February 1997
// G4OpenGLXmViewer : Class derived from G4OpenGLXViewer, to provide
//                  (Motif) widget OpenGL functionality for GEANT4.

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "globals.hh"

#include "G4OpenGLXmViewer.hh"

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

#include "G4ios.hh"
#include <assert.h>
#include <unistd.h>
#include <stdlib.h>

#include "G4VisExtent.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Point3D.hh"
#include "G4Normal3D.hh"

#include "G4VisManager.hh"
#include "G4OpenGLXmSliderBar.hh"

#include "G4Xt.hh"
#include <X11/Shell.h>

static G4ViewParameters viewingParameters;
void G4OpenGLXmViewerSecondaryLoopPostAction ();

void G4OpenGLXmViewer::ShowView () {

  viewingParameters = fVP;
  G4Xt::getInstance () -> SecondaryLoop ();

}

void G4OpenGLXmViewerSecondaryLoopPostAction ()
{
  if(G4Xt::getInstance () -> GetExitSecondaryLoopCode ()==OGL_EXIT_CODE)
    {
      G4VisManager::GetInstance() -> SetCurrentViewParameters() = viewingParameters;
    }
}


void G4OpenGLXmViewer::GetXmConnection () {
  
  G4Xt* interactorManager = G4Xt::getInstance ();
  toplevel = (Widget)interactorManager->GetMainInteractor();
  app      = XtWidgetToApplicationContext(toplevel);

  if (!toplevel) {
    fViewId = -1;  // This flags an error.
    G4cerr << "G4OpenGLXmViewer::GetXmConnection unable to Initialize"
      " application context." << G4endl;
    return;
  }

  // Better to put this in an X11 resource file !!!
  interactorManager->PutStringInResourceDatabase ((char*)"\
*glxarea*width: 500\n\
*glxarea*height: 500\n\
*frame*x: 10\n\
*frame*y: 10\n\
*frame*topOffset: 10\n\
*frame*bottomOffset: 10\n\
*frame*rightOffset: 10\n\
*frame*leftOffset: 10\n\
*frame*shadowType: SHADOW_IN\n\
*frame*useColorObj: False\n\
*frame*primaryColorSetId: 3\n\
*frame*secondaryColorSetId: 3\n\
*menubar*useColorObj: False\n\
*menubar*primaryColorSetId: 3\n\
*menubar*secondaryColorSetId: 3\n\
*toplevel*useColorObj: False\n\
*toplevel*primaryColorSetId: 3\n\
*toplevel*secondaryColorSetId: 3\n\
");
  interactorManager->AddSecondaryLoopPostAction ((G4SecondaryLoopAction)G4OpenGLXmViewerSecondaryLoopPostAction);
  
  shell = XtAppCreateShell ((String)fName.data(),(String)fName.data(),topLevelShellWidgetClass,XtDisplay(toplevel),NULL,0); 
  interactorManager->AddShell (shell);

  dpy = XtDisplay (shell);

  if (!dpy) {
    fViewId = -1;  // This flags an error.
    G4cerr << "G4OpenGLXmViewer::GetXmConnection unable to connect to display."
	 << G4endl;
    return;
  }

  if (!glXQueryExtension (dpy, &errorBase, &eventBase)) {
    fViewId = -1;  // This flags an error.
    G4cerr << "G4OpenGLXmViewer::GetXmConnection. X Server has no GLX extension."
	 << G4endl;;
    return;
  }
}

void G4OpenGLXmViewer::CreateMainWindow () {

  bgnd = XWhitePixelOfScreen (XtScreen(shell));
  borcol = XBlackPixelOfScreen (XtScreen(shell));
  
  XtVaSetValues (shell, 
		 XtNvisual, vi -> visual, 
       		 XtNdepth, vi -> depth, 
       		 XtNcolormap, cmap, 
		 XtNborderColor, &borcol,
		 XtNbackground, &bgnd,
		 XmNtitle, fName.data(),
		 NULL);

  main_win = XtVaCreateManagedWidget ("main_win", 
				      xmMainWindowWidgetClass,
				      shell,
				      XtNvisual, vi -> visual, 
				      XtNdepth, vi -> depth, 
				      XtNcolormap, cmap, 
				      XtNborderColor, borcol,
				      XtNbackground, bgnd,
				      NULL);

  //*********Create a menu bar for the window********
  style_str = XmStringCreateLocalized ((char*)"Style");
  actions_str = XmStringCreateLocalized ((char*)"Actions");
  misc_str = XmStringCreateLocalized ((char*)"Miscellany");
  spec_str = XmStringCreateLocalized ((char*)"Special");

  menubar = XmVaCreateSimpleMenuBar (main_win,
				     (char*)"menubar",
				     XmVaCASCADEBUTTON, style_str, 'S',
				     XmVaCASCADEBUTTON, actions_str, 'A',
				     XmVaCASCADEBUTTON, misc_str, 'M',
				     XmVaCASCADEBUTTON, spec_str, 'p',
				     XtNvisual, vi -> visual, 
				     XtNdepth, vi -> depth, 
				     XtNcolormap, cmap, 
				     XtNborderColor, borcol,
				     XtNbackground, bgnd,
				     NULL);

  XmStringFree (style_str);
  XmStringFree (actions_str);
  XmStringFree (misc_str);
  XmStringFree (spec_str);

  G4cout << "Created menubar" << G4endl;


  //*********Create style pulldown menu on menubar*********
  rep_str = XmStringCreateLocalized ((char*)"Representation");
  draw_str = XmStringCreateLocalized ((char*)"Drawing");
  bgnd_str = XmStringCreateLocalized ((char*)"Background color");

  style_cascade = XmVaCreateSimplePulldownMenu
    (menubar,
     (char*)"style",
     0,
     NULL,
     XmVaCASCADEBUTTON, rep_str, 'R',
     XmVaCASCADEBUTTON, draw_str, 'D',
     XmVaCASCADEBUTTON, bgnd_str, 'B',
     XtNvisual, vi -> visual, 
     XtNdepth, vi -> depth, 
     XtNcolormap, cmap, 
     XtNborderColor, borcol,
     XtNbackground, bgnd,
     NULL);
  
  XmStringFree (rep_str);
  XmStringFree (draw_str);
  XmStringFree (bgnd_str);

  //  G4cout << "Created Style pulldown menu" << G4endl;

  //Add Representation pullright menu to style cascade...
  polyhedron_str = XmStringCreateLocalized ((char*)"Polyhedron");
  nurbs_str = XmStringCreateLocalized ((char*)"NURBS");

  rep_style_pullright = XmVaCreateSimplePulldownMenu 
    (style_cascade,
     (char*)"rep_style",
     0,
     G4OpenGLXmViewer::rep_style_callback,
     XmVaRADIOBUTTON, polyhedron_str, 'P', NULL, NULL,
     XmVaRADIOBUTTON, nurbs_str, 'N', NULL, NULL,
     XmNradioBehavior, True, 
     XmNradioAlwaysOne, True, 
     XmNuserData, this, 
     XtNvisual, vi -> visual, 
     XtNdepth, vi -> depth, 
     XtNcolormap, cmap, 
     XtNborderColor, borcol,
     XtNbackground, bgnd,
     NULL);
  
  Widget special_widget;

  G4ViewParameters::RepStyle style;
  style = fVP.GetRepStyle();
  
  if (style == G4ViewParameters::polyhedron) {
    if(special_widget = XtNameToWidget(rep_style_pullright, "button_0")) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else if (style == G4ViewParameters::nurbs) {
    if(special_widget = XtNameToWidget(rep_style_pullright, "button_1")) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else {
    G4Exception("Invalid Representation style in G4OpenGLXmViewer::CreateContext");
  }
  XmStringFree (polyhedron_str);
  XmStringFree (nurbs_str);
  
  //  G4cout << "Created Representation pulldown menu" << G4endl;

  //Add Drawing pullright menu to style cascade...
  wireframe_str = XmStringCreateLocalized ((char*)"Wireframe");
  hlr_str = XmStringCreateLocalized ((char*)"Hidden line removal");
  hsr_str = XmStringCreateLocalized ((char*)"Hidden surface removal");
  hlhsr_str = XmStringCreateLocalized ((char*)"Hidden line and surface removal");

  drawing_style_pullright = XmVaCreateSimplePulldownMenu 
    (style_cascade,
     (char*)"drawing_style",
     1,
     G4OpenGLXmViewer::drawing_style_callback,
     XmVaRADIOBUTTON, wireframe_str, 'W', NULL, NULL,
     XmVaRADIOBUTTON, hlr_str, 'L', NULL, NULL,
     XmVaRADIOBUTTON, hsr_str, 'S', NULL, NULL,
     XmVaRADIOBUTTON, hlhsr_str, 'H', NULL, NULL,
     XmNradioBehavior, True, 
     XmNradioAlwaysOne, True, 
     XmNuserData, this,
     XtNvisual, vi -> visual, 
     XtNdepth, vi -> depth, 
     XtNcolormap, cmap, 
     XtNborderColor, borcol,
     XtNbackground, bgnd,
     NULL);
  
  G4ViewParameters::DrawingStyle d_style;
  d_style = fVP.GetDrawingStyle();
  
  if (d_style == G4ViewParameters::wireframe) {
    if(special_widget = XtNameToWidget(drawing_style_pullright, "button_0")) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else if (d_style == G4ViewParameters::hlr) {
    if(special_widget = XtNameToWidget(drawing_style_pullright, "button_1")) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else if (d_style == G4ViewParameters::hsr) {
    if(special_widget = XtNameToWidget(drawing_style_pullright, "button_2")) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else if (d_style == G4ViewParameters::hlhsr) {
    if(special_widget = XtNameToWidget(drawing_style_pullright, "button_3")) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else {
    G4Exception("Invalid Drawing style in G4OpenGLXmViewer::CreateContext");
  }

  XmStringFree (wireframe_str);
  XmStringFree (hlr_str);
  XmStringFree (hsr_str);
  XmStringFree (hlhsr_str);

  //  G4cout << "Created Drawing pullright menu" << G4endl;

  //Add Drawing pullright menu to style cascade...
  white_str = XmStringCreateLocalized ((char*)"White");
  black_str = XmStringCreateLocalized ((char*)"Black");

  background_color_pullright = XmVaCreateSimplePulldownMenu 
    (style_cascade,
     (char*)"background_color",
     2,
     G4OpenGLXmViewer::background_color_callback,
     XmVaRADIOBUTTON, white_str, 'W', NULL, NULL,
     XmVaRADIOBUTTON, black_str, 'B', NULL, NULL,
     XmNradioBehavior, True, 
     XmNradioAlwaysOne, True, 
     XmNuserData, this,
     XtNvisual, vi -> visual, 
     XtNdepth, vi -> depth, 
     XtNcolormap, cmap, 
     XtNborderColor, borcol,
     XtNbackground, bgnd,
     NULL);
  
  if (white_background == true) {
    if(special_widget = XtNameToWidget(background_color_pullright, 
				       "button_0")) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else if (white_background == false) {
    if(special_widget = XtNameToWidget(background_color_pullright,
				       "button_1")) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else {
    G4Exception("white_background in G4OpenGLXmViewer is neither true nor false!!");
  }

  XmStringFree (white_str);
  XmStringFree (black_str);

  //  G4cout << "Created Background color pullright menu" << G4endl;

  //*********Create actions pulldown menu on menubar*********
  rot_str = XmStringCreateLocalized ((char*)"Rotation control panel");
  pan_str = XmStringCreateLocalized ((char*)"Panning control panel");
  set_str = XmStringCreateLocalized ((char*)"Set control panel limits");

  actions_cascade = XmVaCreateSimplePulldownMenu
    (menubar,
     (char*)"actions",
     1,
     G4OpenGLXmViewer::actions_callback,
     XmVaPUSHBUTTON, rot_str, 'R', NULL, NULL,
     XmVaPUSHBUTTON, pan_str, 'P', NULL, NULL,
     XmVaPUSHBUTTON, set_str, 'S', NULL, NULL,
     XmNuserData, this, 
     XtNvisual, vi -> visual, 
     XtNdepth, vi -> depth, 
     XtNcolormap, cmap, 
     XtNborderColor, borcol,
     XtNbackground, bgnd,
     NULL);
  
  XmStringFree (rot_str);
  XmStringFree (pan_str);
  XmStringFree (set_str);
  G4cout << "Created Actions pulldown menu" << G4endl;

  misc_str = XmStringCreateLocalized ((char*)"Miscellany control panel");
  exit_str = XmStringCreateLocalized ((char*)"Exit to G4Vis>");
  print_str = XmStringCreateLocalized ((char*)"Create .eps file");

  //*********Create miscellany pulldown menu on menubar*********
  misc_cascade = XmVaCreateSimplePulldownMenu
    (menubar,
     (char*)"miscellany",
     2,
     G4OpenGLXmViewer::misc_callback,
     XmVaPUSHBUTTON, misc_str, 'M', NULL, NULL,
     XmVaPUSHBUTTON, exit_str, 'E', NULL, NULL,
     XmVaPUSHBUTTON, print_str, 'P', NULL, NULL,
     XmNuserData, this,
     XtNvisual, vi -> visual, 
     XtNdepth, vi -> depth, 
     XtNcolormap, cmap, 
     XtNborderColor, borcol,
     XtNbackground, bgnd,
     NULL);
  
  XmStringFree (misc_str);
  XmStringFree (exit_str);
  XmStringFree (print_str);
  G4cout << "Created Miscellany pulldown menu" << G4endl;

  trans_str = XmStringCreateLocalized ((char*)"Transparency");
  anti_str = XmStringCreateLocalized ((char*)"Antialiasing");
  halo_str = XmStringCreateLocalized ((char*)"Haloing");

  //*********Create special pulldown menu on menubar*********
  spec_cascade = XmVaCreateSimplePulldownMenu
    (menubar,
     (char*)"special",
     3,
     NULL,
     XmVaCASCADEBUTTON, trans_str, 'T',
     XmVaCASCADEBUTTON, anti_str, 'A',
     XmVaCASCADEBUTTON, halo_str, 'H',
     XtNvisual, vi -> visual, 
     XtNdepth, vi -> depth, 
     XtNcolormap, cmap, 
     XtNborderColor, borcol,
     XtNbackground, bgnd,
     NULL);
  
  XmStringFree (trans_str);
  XmStringFree (anti_str);
  XmStringFree (halo_str);

  //  G4cout << "Created Special pulldown menu" << G4endl;

  //Add Transparency pullright menu to special cascade...
  off_str = XmStringCreateLocalized ((char*)"Off");
  on_str = XmStringCreateLocalized ((char*)"On");

  transparency_pullright = XmVaCreateSimplePulldownMenu 
    (spec_cascade,
     (char*)"transparency",
     0,
     G4OpenGLXmViewer::transparency_callback,
     XmVaRADIOBUTTON, off_str, 'f', NULL, NULL,
     XmVaRADIOBUTTON, on_str, 'n', NULL, NULL,
     XmNradioBehavior, True, 
     XmNradioAlwaysOne, True, 
     XmNuserData, this,
     XtNvisual, vi -> visual, 
     XtNdepth, vi -> depth, 
     XtNcolormap, cmap, 
     XtNborderColor, borcol,
     XtNbackground, bgnd,
     NULL);
  
  if (transparency_enabled == false) {
    if(special_widget = XtNameToWidget(transparency_pullright, 
				       "button_0")) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else if (transparency_enabled == true) {
    if(special_widget = XtNameToWidget(transparency_pullright,
				       "button_1")) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else {
    G4Exception("transparency_enabled in G4OpenGLXmViewer is neither true nor false!!");
  }

  //Add antialias pullright menu to special cascade...
  antialias_pullright = XmVaCreateSimplePulldownMenu 
    (spec_cascade,
     (char*)"antialias",
     1,
     G4OpenGLXmViewer::antialias_callback,
     XmVaRADIOBUTTON, off_str, 'f', NULL, NULL,
     XmVaRADIOBUTTON, on_str, 'n', NULL, NULL,
     XmNradioBehavior, True, 
     XmNradioAlwaysOne, True, 
     XmNuserData, this,
     XtNvisual, vi -> visual, 
     XtNdepth, vi -> depth, 
     XtNcolormap, cmap, 
     XtNborderColor, borcol,
     XtNbackground, bgnd,
     NULL);
  
  if (antialiasing_enabled == false) {
    if(special_widget = XtNameToWidget(antialias_pullright, 
				       "button_0")) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else if (antialiasing_enabled == true) {
    if(special_widget = XtNameToWidget(antialias_pullright,
				       "button_1")) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else {
    G4Exception("antialiasing_enabled in G4OpenGLXmViewer is neither true nor false!!");
  }

  //Add Haloing pullright menu to special cascade...
  haloing_pullright = XmVaCreateSimplePulldownMenu 
    (spec_cascade,
     (char*)"haloing",
     2,
     G4OpenGLXmViewer::haloing_callback,
     XmVaRADIOBUTTON, off_str, 'f', NULL, NULL,
     XmVaRADIOBUTTON, on_str, 'n', NULL, NULL,
     XmNradioBehavior, True, 
     XmNradioAlwaysOne, True, 
     XmNuserData, this,
     XtNvisual, vi -> visual, 
     XtNdepth, vi -> depth, 
     XtNcolormap, cmap, 
     XtNborderColor, borcol,
     XtNbackground, bgnd,
     NULL);
  
  if (haloing_enabled == false) {
    if(special_widget = XtNameToWidget(haloing_pullright, 
				       "button_0")) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else if (haloing_enabled == true) {
    if(special_widget = XtNameToWidget(haloing_pullright,
				       "button_1")) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else {
    G4Exception("haloing_enabled in G4OpenGLXmViewer is neither true nor false!!");
  }

  XtManageChild (menubar);
  frame = XtVaCreateManagedWidget ((char*)"frame",
				   xmFrameWidgetClass, main_win,
				   XtNvisual, vi -> visual, 
				   XtNdepth, vi -> depth, 
				   XtNcolormap, cmap, 
				   XtNborderColor, borcol,
				   XtNbackground, bgnd,
				   NULL);

  glxarea = XtVaCreateManagedWidget ((char*)"glxarea", 
				     xmDrawingAreaWidgetClass,
				     frame,
				     XtNvisual, vi -> visual, 
				     XtNdepth, vi -> depth, 
				     XtNcolormap, cmap, 
				     XtNborderColor, borcol,
				     XtNbackground, bgnd,
				     NULL);
  
  XtAddCallback (glxarea, 
		 XmNexposeCallback, 
		 G4OpenGLXmViewer::expose_callback, 
		 this);

  XtAddCallback (glxarea, 
		 XmNresizeCallback, 
		 G4OpenGLXmViewer::resize_callback, 
		 this);

  XmMainWindowSetAreas (main_win,  // main widget, children are specified 
			menubar,   // widget to use as menu bar
			NULL,      // widget to use as command window
			NULL,      // widget for horizontal scroll bar
			NULL,      // widget for vertical scroll bar
			frame      // widget to be used for work window
			);

  XtRealizeWidget(shell);
  
  // Once widget is realized (ie, associated with a created X window), we
  // can bind the OpenGL rendering context to the window.

  Dimension width, height;
  XtVaGetValues (glxarea,XmNwidth,&width,XmNheight,&height,NULL);
  WinSize_x = (G4int) width;
  WinSize_y = (G4int) height;

  win = XtWindow (glxarea);

  glXMakeCurrent (dpy, win, cx);

}

G4OpenGLXmViewer::G4OpenGLXmViewer (G4OpenGLSceneHandler& scene):
G4VViewer (scene, -1),
G4OpenGLViewer (scene),
G4OpenGLXViewer (scene),
zoom_low (fVP.GetZoomFactor() / 10.0),
zoom_high (fVP.GetZoomFactor() * 10.0),
dolly_low (fVP.GetDolly() - 1000.0),
dolly_high (fVP.GetDolly() + 1000.0),
rot_sens (4.),
frameNo (0),
wob_sens (20.),
rot_sens_limit (90.),
pan_sens_limit (100.),
fov (0.0),
original_vp(fVP.GetViewpointDirection()),
fppanning_top (NULL),
fprotation_top (NULL),
fpmiscellany_top (NULL),
fpsetting_top (NULL),
fpprint_top (NULL),
fppanning_slider (NULL),
fprotation_slider (NULL),
fpzoom_slider (NULL),
fpdolly_slider (NULL)
{

  WinSize_x = 100;
  WinSize_y = 100;
  
  GetXmConnection ();
  if (fViewId < 0) return;

}

G4OpenGLXmViewer::~G4OpenGLXmViewer ()
{
  XtDestroyWidget  (shell);
  win = 0; // ...to avoid XDestroyWindow in G4OpenGLXViewer base class
  // because XtDestroyWidget has already destroyed it.
  G4Xt::getInstance () ->RemoveShell (shell);

/******************************
  if (fprotation_top) {
    delete fprotation_top;
  }

  if (fppanning_top) {
    delete fppanning_top;
  }

  if (fpsetting_top) {
    delete fpsetting_top;
  }

  if (fpmiscellany_top) {
    delete fpmiscellany_top;
  }
******************************/
}

#endif
