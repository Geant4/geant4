// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmView.cc,v 1.1 1999-01-07 16:15:04 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  10th February 1997
// G4OpenGLXmView : Class derived from G4OpenGLXView, to provide
//                  (Motif) widget OpenGL functionality for GEANT4.

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "globals.hh"

#include "G4OpenGLXmView.hh"

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
static void SecondaryLoopPostAction ();

void G4OpenGLXmView::ShowView () {

  viewingParameters = fVP;
  G4Xt::getInstance () -> SecondaryLoop ();

}

static void SecondaryLoopPostAction ()
{
  if(G4Xt::getInstance () -> GetExitSecondaryLoopCode ()==OGL_EXIT_CODE)
    {
      G4VisManager::GetInstance() -> SetCurrentViewParameters() = viewingParameters;
    }
}


void G4OpenGLXmView::GetXmConnection () {
  
  G4Xt* interactorManager = G4Xt::getInstance ();
  toplevel = (Widget)interactorManager->GetMainInteractor();
  app      = XtWidgetToApplicationContext(toplevel);

  if (!toplevel) {
    fViewId = -1;  // This flags an error.
    G4cerr << "G4OpenGLXmView::GetXmConnection unable to Initialize"
      " application context." << endl;
    return;
  }

  // Better to put this in an X11 resource file !!!
  interactorManager->PutStringInResourceDatabase ("\
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
  interactorManager->AddSecondaryLoopPostAction ((G4SecondaryLoopAction)SecondaryLoopPostAction);
  
  shell = XtAppCreateShell ((String)fName.data(),(String)fName.data(),topLevelShellWidgetClass,XtDisplay(toplevel),NULL,0); 
  interactorManager->AddShell (shell);

  dpy = XtDisplay (shell);

  if (!dpy) {
    fViewId = -1;  // This flags an error.
    G4cerr << "G4OpenGLXmView::GetXmConnection unable to connect to display."
	 << endl;
    return;
  }

  if (!glXQueryExtension (dpy, &errorBase, &eventBase)) {
    fViewId = -1;  // This flags an error.
    G4cerr << "G4OpenGLXmView::GetXmConnection. X Server has no GLX extension."
	 << endl;;
    return;
  }
}

void G4OpenGLXmView::CreateMainWindow () {

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
  style_str = XmStringCreateLocalized ("Style");
  actions_str = XmStringCreateLocalized ("Actions");
  misc_str = XmStringCreateLocalized ("Miscellany");
  spec_str = XmStringCreateLocalized ("Special");

  menubar = XmVaCreateSimpleMenuBar (main_win,
				     "menubar",
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

  G4cout << "Created menubar" << endl;


  //*********Create style pulldown menu on menubar*********
  rep_str = XmStringCreateLocalized ("Representation");
  draw_str = XmStringCreateLocalized ("Drawing");
  bgnd_str = XmStringCreateLocalized ("Background color");

  style_cascade = XmVaCreateSimplePulldownMenu
    (menubar,
     "style",
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

  //  G4cout << "Created Style pulldown menu" << endl;

  //Add Representation pullright menu to style cascade...
  polyhedron_str = XmStringCreateLocalized ("Polyhedron");
  nurbs_str = XmStringCreateLocalized ("NURBS");

  rep_style_pullright = XmVaCreateSimplePulldownMenu 
    (style_cascade,
     "rep_style",
     0,
     G4OpenGLXmView::rep_style_callback,
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
    G4Exception("Invalid Representation style in G4OpenGLXmView::CreateContext");
  }
  XmStringFree (polyhedron_str);
  XmStringFree (nurbs_str);
  
  //  G4cout << "Created Representation pulldown menu" << endl;

  //Add Drawing pullright menu to style cascade...
  wireframe_str = XmStringCreateLocalized ("Wireframe");
  hlr_str = XmStringCreateLocalized ("Hidden line removal");
  hsr_str = XmStringCreateLocalized ("Hidden surface removal");
  hlhsr_str = XmStringCreateLocalized ("Hidden line and surface removal");

  drawing_style_pullright = XmVaCreateSimplePulldownMenu 
    (style_cascade,
     "drawing_style",
     1,
     G4OpenGLXmView::drawing_style_callback,
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
    G4Exception("Invalid Drawing style in G4OpenGLXmView::CreateContext");
  }

  XmStringFree (wireframe_str);
  XmStringFree (hlr_str);
  XmStringFree (hsr_str);
  XmStringFree (hlhsr_str);

  //  G4cout << "Created Drawing pullright menu" << endl;

  //Add Drawing pullright menu to style cascade...
  white_str = XmStringCreateLocalized ("White");
  black_str = XmStringCreateLocalized ("Black");

  background_color_pullright = XmVaCreateSimplePulldownMenu 
    (style_cascade,
     "background_color",
     2,
     G4OpenGLXmView::background_color_callback,
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
    G4Exception("white_background in G4OpenGLXmView is neither true nor false!!");
  }

  XmStringFree (white_str);
  XmStringFree (black_str);

  //  G4cout << "Created Background color pullright menu" << endl;

  //*********Create actions pulldown menu on menubar*********
  rot_str = XmStringCreateLocalized ("Rotation control panel");
  pan_str = XmStringCreateLocalized ("Panning control panel");
  set_str = XmStringCreateLocalized ("Set control panel limits");

  actions_cascade = XmVaCreateSimplePulldownMenu
    (menubar,
     "actions",
     1,
     G4OpenGLXmView::actions_callback,
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
  G4cout << "Created Actions pulldown menu" << endl;

  misc_str = XmStringCreateLocalized ("Miscellany control panel");
  exit_str = XmStringCreateLocalized ("Exit to G4Vis>");
  print_str = XmStringCreateLocalized ("Create .eps file");

  //*********Create miscellany pulldown menu on menubar*********
  misc_cascade = XmVaCreateSimplePulldownMenu
    (menubar,
     "miscellany",
     2,
     G4OpenGLXmView::misc_callback,
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
  G4cout << "Created Miscellany pulldown menu" << endl;

  trans_str = XmStringCreateLocalized ("Transparency");
  anti_str = XmStringCreateLocalized ("Antialiasing");
  halo_str = XmStringCreateLocalized ("Haloing");

  //*********Create special pulldown menu on menubar*********
  spec_cascade = XmVaCreateSimplePulldownMenu
    (menubar,
     "special",
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

  //  G4cout << "Created Special pulldown menu" << endl;

  //Add Transparency pullright menu to special cascade...
  off_str = XmStringCreateLocalized ("Off");
  on_str = XmStringCreateLocalized ("On");

  transparency_pullright = XmVaCreateSimplePulldownMenu 
    (spec_cascade,
     "transparency",
     0,
     G4OpenGLXmView::transparency_callback,
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
    G4Exception("transparency_enabled in G4OpenGLXmView is neither true nor false!!");
  }

  //Add antialias pullright menu to special cascade...
  antialias_pullright = XmVaCreateSimplePulldownMenu 
    (spec_cascade,
     "antialias",
     1,
     G4OpenGLXmView::antialias_callback,
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
    G4Exception("antialiasing_enabled in G4OpenGLXmView is neither true nor false!!");
  }

  //Add Haloing pullright menu to special cascade...
  haloing_pullright = XmVaCreateSimplePulldownMenu 
    (spec_cascade,
     "haloing",
     2,
     G4OpenGLXmView::haloing_callback,
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
    G4Exception("haloing_enabled in G4OpenGLXmView is neither true nor false!!");
  }

  XtManageChild (menubar);
  frame = XtVaCreateManagedWidget ("frame", xmFrameWidgetClass, main_win,
				   XtNvisual, vi -> visual, 
				   XtNdepth, vi -> depth, 
				   XtNcolormap, cmap, 
				   XtNborderColor, borcol,
				   XtNbackground, bgnd,
				   NULL);

  glxarea = XtVaCreateManagedWidget ("glxarea", 
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
		 G4OpenGLXmView::expose_callback, 
		 this);

  XtAddCallback (glxarea, 
		 XmNresizeCallback, 
		 G4OpenGLXmView::resize_callback, 
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

G4OpenGLXmView::G4OpenGLXmView (G4OpenGLScene& scene):
G4VView (scene, -1),
G4OpenGLView (scene),
G4OpenGLXView (scene),
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
  vi_immediate = 0;
  vi_stored = 0;
  
  GetXmConnection ();
  if (fViewId < 0) return;

  // Try for a visual suitable for OpenGLImmediate..
  // first try for a single buffered RGB window
  if (vi_immediate =
      glXChooseVisual (dpy, XDefaultScreen (dpy), snglBuf_RGBA)) {
    attributeList = snglBuf_RGBA;
    doublebuffer = false;
  }

  if (!vi_immediate){
    // next try for a double buffered RGB, but Draw to top buffer
    if (vi_immediate =
	glXChooseVisual (dpy, XDefaultScreen (dpy), dblBuf_RGBA)) {
      attributeList = dblBuf_RGBA;
      doublebuffer = true;
    }
  }

  // Now try for a visual suitable for OpenGLStored...
  // Try for a double buffered RGB window
  if (vi_stored = glXChooseVisual (dpy, XDefaultScreen (dpy), dblBuf_RGBA)) {
    attributeList = dblBuf_RGBA;
    doublebuffer = true;
  }

}
G4OpenGLXmView::~G4OpenGLXmView ()
{
  XtDestroyWidget  (shell);
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
