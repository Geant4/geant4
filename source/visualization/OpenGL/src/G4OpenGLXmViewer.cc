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
//
// 
// Andrew Walkden  10th February 1997
// G4OpenGLXmViewer : Class derived from G4OpenGLXViewer, to provide
//                  (Motif) widget OpenGL functionality for GEANT4.

#include "globals.hh"

#include "G4OpenGLXmViewer.hh"
#include "G4OpenGLSceneHandler.hh"

#include "G4VisExtent.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Point3D.hh"
#include "G4Normal3D.hh"

#include "G4Scene.hh"

#include "G4OpenGLXmSliderBar.hh"
#include "G4OpenGLXmTextField.hh"

#include "G4Xt.hh"
#include <X11/Shell.h>
#include <Xm/MainW.h>
#include <Xm/Frame.h>
#include <Xm/DrawingA.h>

#include <sstream>

void G4OpenGLXmViewer::ShowView () {

//  glXWaitGL (); //Wait for effects of all previous OpenGL commands to
                //be propagated before progressing.
// JA: Commented out July 2021 - slows rendering down in some cases and I
// don't see any adverse effects.

  glFlush ();

  G4Xt::getInstance () -> SecondaryLoop ();

}

void G4OpenGLXmViewer::ResetView () {
  // reset global parameters
  G4OpenGLViewer::ResetView();

  //reset Xm parameteres
  zoom_high  = fVP.GetZoomFactor() * 10.0;
  zoom_low = fVP.GetZoomFactor() / 10.0;
  rot_sens_limit = 90.;
  wob_low = 0.;
  wob_high = 50.;
  wob_sens = 20.;
  
  bool firstInit = true;
  if (GetSceneHandler() != NULL) {
    if (GetSceneHandler()->GetScene() != NULL) {
      firstInit = false;
    }
  }
  if (firstInit) {
    pan_sens_limit = 100.;
    fPan_sens = pan_sens_limit / 10.0;
    dolly_low  = fVP.GetDolly() - 1000.0;
    dolly_high = fVP.GetDolly() + 1000.0;
  } else {
    fPan_sens = GetSceneHandler()->GetScene()->GetExtent().GetExtentRadius() / 10.0;
    pan_sens_limit = GetSceneHandler()->GetScene()->GetExtent().GetExtentRadius();
    
    dolly_high = GetSceneHandler()->GetScene()->GetExtent().GetExtentRadius();
    dolly_low = -(GetSceneHandler()->GetScene()->GetExtent().GetExtentRadius());
  }

  UpdateControlPanel ();


  // FIXME : L.Garnier 12 Oct 2011
  // Has also to change the Camera/Object, but tricky to do...  

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
  std::ostringstream oss;
  oss <<
    "*glxarea*width: " << fVP.GetWindowSizeHintX() << "\n"
    "*glxarea*height: " << fVP.GetWindowSizeHintY() << "\n"
    /*
    // Tried this as a replacement for the above two lines, but
    // sub-windows (rotation, etc.) came same size!!
    "*geometry: " << fVP.GetXGeometryString() << "\n"
    */
    "*frame*x: 10\n"
    "*frame*y: 10\n"
    "*frame*topOffset: 10\n"
    "*frame*bottomOffset: 10\n"
    "*frame*rightOffset: 10\n"
    "*frame*leftOffset: 10\n"
    "*frame*shadowType: SHADOW_IN\n"
    "*frame*useColorObj: False\n"
    "*frame*primaryColorSetId: 3\n"
    "*frame*secondaryColorSetId: 3\n"
    "*menubar*useColorObj: False\n"
    "*menubar*primaryColorSetId: 3\n"
    "*menubar*secondaryColorSetId: 3\n"
    "*toplevel*useColorObj: False\n"
    "*toplevel*primaryColorSetId: 3\n"
    "*toplevel*secondaryColorSetId: 3\n";
  interactorManager->PutStringInResourceDatabase ((char*)oss.str().c_str());

  //  interactorManager->AddSecondaryLoopPostAction ((G4SecondaryLoopAction)G4OpenGLXmViewerSecondaryLoopPostAction);
  
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
  
  ResizeWindow(fVP.GetWindowSizeHintX(),fVP.GetWindowSizeHintY());

  G4int x_origin = fVP.GetWindowAbsoluteLocationHintX(DisplayWidth(dpy, vi -> screen));

  // FIXME,  screen size != window size on MAC, but I don't know have to get the menuBar
  // size on MAC. L.Garnier 01/2009
  G4int y_origin = fVP.GetWindowAbsoluteLocationHintY(DisplayHeight(dpy, vi -> screen));

  if (fVP.IsWindowSizeHintX () && fVP.IsWindowLocationHintX () && fVP.IsWindowLocationHintY ()) {
    XtVaSetValues (shell, 
                   XtNvisual, vi -> visual, 
                   XtNdepth, vi -> depth,
                   XtNcolormap, cmap, 
                   XtNwidth, getWinWidth(),
                   XtNheight, getWinHeight(),
                   XtNx, x_origin,
                   XtNy, y_origin,
                   XtNborderColor, &borcol,
                   XtNbackground, &bgnd,
                   XmNtitle, fName.data(),
                   NULL);
  } else if (fVP.IsWindowSizeHintX () && !(fVP.IsWindowLocationHintX () || fVP.IsWindowLocationHintY ())) {
    XtVaSetValues (shell, 
                   XtNvisual, vi -> visual, 
                   XtNdepth, vi -> depth,
                   XtNcolormap, cmap, 
                   XtNwidth, getWinWidth(),
                   XtNheight, getWinHeight(),
                   XtNborderColor, &borcol,
                   XtNbackground, &bgnd,
                   XmNtitle, fName.data(),
                   NULL);
  } else if ((!fVP.IsWindowSizeHintX ()) && fVP.IsWindowLocationHintX () && fVP.IsWindowLocationHintY ()) {
    XtVaSetValues (shell, 
                   XtNvisual, vi -> visual, 
                   XtNdepth, vi -> depth,
                   XtNcolormap, cmap, 
                   XtNx, x_origin,
                   XtNy, y_origin,
                   XtNborderColor, &borcol,
                   XtNbackground, &bgnd,
                   XmNtitle, fName.data(),
                   NULL);
  } else {
    XtVaSetValues (shell, 
                   XtNvisual, vi -> visual, 
                   XtNdepth, vi -> depth,
                   XtNcolormap, cmap, 
                   XtNborderColor, &borcol,
                   XtNbackground, &bgnd,
                   XmNtitle, fName.data(),
                   NULL);
  }


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
  style_str = XmStringCreateLocalized ((char*)menu_str[0].c_str());
  actions_str = XmStringCreateLocalized ((char*)menu_str[2].c_str());
  misc_str = XmStringCreateLocalized ((char*)menu_str[4].c_str());
  spec_str = XmStringCreateLocalized ((char*)menu_str[6].c_str());

  menubar = XmVaCreateSimpleMenuBar (main_win,
				     (char*)menu_str[8].c_str(),
				     XmVaCASCADEBUTTON, style_str,   (KeySym)XK_S,  /*G.Barrand : cast to KeySym and use XK_*/
				     XmVaCASCADEBUTTON, actions_str, (KeySym)XK_A,
				     XmVaCASCADEBUTTON, misc_str,    (KeySym)XK_M,
				     XmVaCASCADEBUTTON, spec_str,    (KeySym)XK_p,
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
  draw_str = XmStringCreateLocalized ((char*)menu_str[9].c_str());
  bgnd_str = XmStringCreateLocalized ((char*)menu_str[10].c_str());

  style_cascade = XmVaCreateSimplePulldownMenu
    (menubar,
     (char*)menu_str[1].c_str(),
     0,
     NULL,
     XmVaCASCADEBUTTON, draw_str, (KeySym)XK_D,
     XmVaCASCADEBUTTON, bgnd_str, (KeySym)XK_B,
     XtNvisual, vi -> visual, 
     XtNdepth, vi -> depth, 
     XtNcolormap, cmap, 
     XtNborderColor, borcol,
     XtNbackground, bgnd,
     NULL);
  
  XmStringFree (draw_str);
  XmStringFree (bgnd_str);

  //  G4cout << "Created Style pulldown menu" << G4endl;

  //Add Drawing pullright menu to style cascade...
  wireframe_str = XmStringCreateLocalized ((char*)menu_str[11].c_str());
  hlr_str = XmStringCreateLocalized ((char*)menu_str[12].c_str());
  hsr_str = XmStringCreateLocalized ((char*)menu_str[13].c_str());
  hlhsr_str = XmStringCreateLocalized ((char*)menu_str[14].c_str());

  drawing_style_pullright = XmVaCreateSimplePulldownMenu 
    (style_cascade,
     (char*)menu_str[15].c_str(),
     1,
     drawing_style_callback,
     XmVaRADIOBUTTON, wireframe_str, (KeySym)XK_W, NULL, NULL,
     XmVaRADIOBUTTON, hlr_str,       (KeySym)XK_L, NULL, NULL,
     XmVaRADIOBUTTON, hsr_str,       (KeySym)XK_S, NULL, NULL,
     XmVaRADIOBUTTON, hlhsr_str,     (KeySym)XK_H, NULL, NULL,
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

  G4ViewParameters::DrawingStyle d_style;
  d_style = fVP.GetDrawingStyle();
  
  if (d_style == G4ViewParameters::wireframe) {
    special_widget = XtNameToWidget(drawing_style_pullright, "button_0");
    if(special_widget) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else if (d_style == G4ViewParameters::hlr) {
    special_widget = XtNameToWidget(drawing_style_pullright, "button_1");
    if(special_widget) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else if (d_style == G4ViewParameters::hsr) {
    special_widget = XtNameToWidget(drawing_style_pullright, "button_2");
    if(special_widget) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else if (d_style == G4ViewParameters::hlhsr) {
    special_widget = XtNameToWidget(drawing_style_pullright, "button_3");
    if(special_widget) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else {
    G4Exception
      ("G4OpenGLXmViewer::CreateMainWindow",
       "opengl2015", FatalException,
       "Invalid Drawing style in G4OpenGLXmViewer::CreateContext");
  }

  XmStringFree (wireframe_str);
  XmStringFree (hlr_str);
  XmStringFree (hsr_str);
  XmStringFree (hlhsr_str);

  //  G4cout << "Created Drawing pullright menu" << G4endl;

  //Add Drawing pullright menu to style cascade...
  white_str = XmStringCreateLocalized ((char*)menu_str[16].c_str());
  black_str = XmStringCreateLocalized ((char*)menu_str[17].c_str());

  background_color_pullright = XmVaCreateSimplePulldownMenu 
    (style_cascade,
     (char*)menu_str[18].c_str(),
     2,
     background_color_callback,
     XmVaRADIOBUTTON, white_str, (KeySym)XK_W, NULL, NULL,
     XmVaRADIOBUTTON, black_str, (KeySym)XK_B, NULL, NULL,
     XmNradioBehavior, True, 
     XmNradioAlwaysOne, True, 
     XmNuserData, this,
     XtNvisual, vi -> visual, 
     XtNdepth, vi -> depth, 
     XtNcolormap, cmap, 
     XtNborderColor, borcol,
     XtNbackground, bgnd,
     NULL);
  
  if (background.GetRed() == 1. &&
      background.GetGreen() == 1. &&
      background.GetBlue() == 1.) {
    special_widget = XtNameToWidget(background_color_pullright, "button_0");
    if(special_widget) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else {
    special_widget = XtNameToWidget(background_color_pullright, "button_1");
    if(special_widget) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  }

  XmStringFree (white_str);
  XmStringFree (black_str);

  //  G4cout << "Created Background color pullright menu" << G4endl;

  //*********Create actions pulldown menu on menubar*********
  rot_str = XmStringCreateLocalized ((char*)menu_str[19].c_str());
  pan_str = XmStringCreateLocalized ((char*)menu_str[20].c_str());
  set_str = XmStringCreateLocalized ((char*)menu_str[21].c_str());

  actions_cascade = XmVaCreateSimplePulldownMenu
    (menubar,
     (char*)menu_str[3].c_str(),
     1,
     actions_callback,
     XmVaPUSHBUTTON, rot_str, (KeySym)XK_R, NULL, NULL,
     XmVaPUSHBUTTON, pan_str, (KeySym)XK_P, NULL, NULL,
     XmVaPUSHBUTTON, set_str, (KeySym)XK_S, NULL, NULL,
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

  misc_str = XmStringCreateLocalized ((char*)menu_str[22].c_str());
  exit_str = XmStringCreateLocalized ((char*)menu_str[23].c_str());
  print_str = XmStringCreateLocalized ((char*)menu_str[24].c_str());

  //*********Create miscellany pulldown menu on menubar*********
  misc_cascade = XmVaCreateSimplePulldownMenu
    (menubar,
     (char*)menu_str[5].c_str(),
     2,
     misc_callback,
     XmVaPUSHBUTTON, misc_str,  (KeySym)XK_M, NULL, NULL,
     XmVaPUSHBUTTON, exit_str,  (KeySym)XK_E, NULL, NULL,
     XmVaPUSHBUTTON, print_str, (KeySym)XK_P, NULL, NULL,
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

  trans_str = XmStringCreateLocalized ((char*)menu_str[25].c_str());
  anti_str = XmStringCreateLocalized ((char*)menu_str[27].c_str());
  halo_str = XmStringCreateLocalized ((char*)menu_str[29].c_str());
  aux_edge_str = XmStringCreateLocalized ((char*)menu_str[31].c_str());

  //*********Create special pulldown menu on menubar*********
  spec_cascade = XmVaCreateSimplePulldownMenu
    (menubar,
     (char*)menu_str[7].c_str(),
     3,
     NULL,
     XmVaCASCADEBUTTON, trans_str,    (KeySym)XK_T,
     XmVaCASCADEBUTTON, anti_str,     (KeySym)XK_A,
     XmVaCASCADEBUTTON, halo_str,     (KeySym)XK_H,
     XmVaCASCADEBUTTON, aux_edge_str, (KeySym)XK_E,
     XtNvisual, vi -> visual, 
     XtNdepth, vi -> depth, 
     XtNcolormap, cmap, 
     XtNborderColor, borcol,
     XtNbackground, bgnd,
     NULL);
  
  XmStringFree (trans_str);
  XmStringFree (anti_str);
  XmStringFree (halo_str);
  XmStringFree (aux_edge_str);

  //  G4cout << "Created Special pulldown menu" << G4endl;

  //Add Transparency pullright menu to special cascade...
  off_str = XmStringCreateLocalized ((char*)menu_str[33].c_str());
  on_str = XmStringCreateLocalized ((char*)menu_str[34].c_str());

  transparency_pullright = XmVaCreateSimplePulldownMenu 
    (spec_cascade,
     (char*)menu_str[26].c_str(),
     0,
     transparency_callback,
     XmVaRADIOBUTTON, off_str, (KeySym)XK_f, NULL, NULL,
     XmVaRADIOBUTTON, on_str,  (KeySym)XK_n, NULL, NULL,
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
    special_widget = XtNameToWidget(transparency_pullright, "button_0");
    if(special_widget) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else if (transparency_enabled == true) {
    special_widget = XtNameToWidget(transparency_pullright, "button_1");
    if(special_widget) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else {
    G4Exception
      ("G4OpenGLXmViewer::CreateMainWindow",
       "opengl2016", FatalException,
       "transparency_enabled in G4OpenGLXmViewer is neither true nor false!!");
  }

  //Add antialias pullright menu to special cascade...
  antialias_pullright = XmVaCreateSimplePulldownMenu 
    (spec_cascade,
     (char*)menu_str[28].c_str(),
     1,
     antialias_callback,
     XmVaRADIOBUTTON, off_str, (KeySym)XK_f, NULL, NULL,
     XmVaRADIOBUTTON, on_str,  (KeySym)XK_n, NULL, NULL,
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
    special_widget = XtNameToWidget(antialias_pullright, "button_0");
    if(special_widget) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else if (antialiasing_enabled == true) {
    special_widget = XtNameToWidget(antialias_pullright, "button_1");
    if(special_widget) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else {
    G4Exception
      ("G4OpenGLXmViewer::CreateMainWindow",
       "opengl2017", FatalException,
       "antialiasing_enabled in G4OpenGLXmViewer is neither true nor false!!");
  }

  //Add Haloing pullright menu to special cascade...
  haloing_pullright = XmVaCreateSimplePulldownMenu 
    (spec_cascade,
     (char*)menu_str[30].c_str(),
     2,
     haloing_callback,
     XmVaRADIOBUTTON, off_str, (KeySym)XK_f, NULL, NULL,
     XmVaRADIOBUTTON, on_str,  (KeySym)XK_n, NULL, NULL,
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
    special_widget = XtNameToWidget(haloing_pullright, "button_0");
    if(special_widget) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else if (haloing_enabled == true) {
    special_widget = XtNameToWidget(haloing_pullright, "button_1");
    if(special_widget) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else {
    G4Exception
      ("G4OpenGLXmViewer::CreateMainWindow",
       "opengl2018", FatalException,
       "haloing_enabled in G4OpenGLXmViewer is neither true nor false!!");
  }

  //Add Aux_Edge pullright menu to special cascade...
  aux_edge_pullright = XmVaCreateSimplePulldownMenu 
    (spec_cascade,
     (char*)menu_str[32].c_str(),
     3,
     aux_edge_callback,
     XmVaRADIOBUTTON, off_str, (KeySym)XK_f, NULL, NULL,
     XmVaRADIOBUTTON, on_str,  (KeySym)XK_n, NULL, NULL,
     XmNradioBehavior, True, 
     XmNradioAlwaysOne, True, 
     XmNuserData, this,
     XtNvisual, vi -> visual, 
     XtNdepth, vi -> depth, 
     XtNcolormap, cmap, 
     XtNborderColor, borcol,
     XtNbackground, bgnd,
     NULL);
  
  if (!fVP.IsAuxEdgeVisible()) {
    special_widget = XtNameToWidget(aux_edge_pullright, "button_0");
    if(special_widget) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  } else {
    special_widget = XtNameToWidget(aux_edge_pullright, "button_1");
    if(special_widget) {
      XtVaSetValues (special_widget, XmNset, True, NULL);
    }
  }

  XtManageChild (menubar);
  frame = XtVaCreateManagedWidget ((char*)menu_str[35].c_str(),
				   xmFrameWidgetClass, main_win,
				   XtNvisual, vi -> visual, 
				   XtNdepth, vi -> depth, 
				   XtNcolormap, cmap, 
				   XtNborderColor, borcol,
				   XtNbackground, bgnd,
				   NULL);

  glxarea = XtVaCreateManagedWidget ((char*)menu_str[36].c_str(), 
				     xmDrawingAreaWidgetClass,
				     frame,
				     XtNvisual, vi -> visual, 
				     XtNdepth, vi -> depth, 
				     XtNcolormap, cmap, 
				     XtNborderColor, borcol,
				     XtNbackground, bgnd,
				     NULL);
  

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

  win = XtWindow (glxarea);

  glXMakeCurrent (dpy, win, cxMaster);

  // This should be add AFTER glXMakeCurrent done because it will fire a resizeCallback
  XtAddCallback (glxarea, 
 		 XmNresizeCallback, 
 		 resize_callback, 
 		 this);

  XtAddCallback (glxarea, 
		 XmNexposeCallback, 
		 expose_callback, 
		 this);

}

G4OpenGLXmViewer::G4OpenGLXmViewer (G4OpenGLSceneHandler& scene):
G4VViewer (scene, -1),
G4OpenGLViewer (scene),
G4OpenGLXViewer (scene),
toplevel (0),
shell (0),
main_win (0),
menubar (0),
style_cascade (0),
actions_cascade (0),
misc_cascade (0),
spec_cascade (0),
drawing_style_pullright (0),
background_color_pullright (0),
transparency_pullright (0),
antialias_pullright (0),
haloing_pullright (0),
aux_edge_pullright (0),
frame (0),
glxarea (0),
style_str (0),
actions_str (0),
misc_str (0),
spec_str (0),
draw_str (0),
polyhedron_str (0),
wireframe_str (0),
hlr_str (0),
hsr_str (0),
hlhsr_str (0),
set_str (0),
rot_str (0),
pan_str (0),
exit_str (0),
quit_str (0),
print_str (0),
white_str (0),
black_str (0),
anti_str (0),
trans_str (0),
halo_str (0),
aux_edge_str (0),
bgnd_str (0),
off_str (0),
on_str (0),
zoom_high (0.0),
zoom_low (0.0),
pan_low (0.0),
pan_high (0.0),
dolly_low (0.0),
dolly_high (0.0),
fov (0.0),
rot_sens_limit (0.0),
pan_sens_limit (0.0),
wob_high (0.0),
wob_low (0.0),
wob_sens (0.0),
pan_right (false),
rotate_right (false),
pan_up (false),
rotate_up (false),
original_vp(fVP.GetViewpointDirection()),
frameNo (0),
fprotation_top (0),
fprotation_button_box (0),
fprotation_button1 (0),
fprotation_button2 (0),
fprotation_slider_box (0),
fprotation_slider (0),
fprotation_arrow_box (0),
fprotation_arrow (0),
fppanning_top (0),
fppanning_box (0),
fppanning_arrows (0),
fppanning_slider (0),
fpzoom_box (0),
fpzoom_slider (0),
fpdolly_box (0),
fpdolly_slider (0),
fpsetting_top (0),
fpsetting_box (0),
fppan_set (0),
fprot_set (0),
fpzoom_upper (0),
fpzoom_lower (0),
fpdolly_upper (0),
fpdolly_lower (0),
fpok_button (0),
fpmiscellany_top (0),
fpwobble_box (0),
fpwobble_button (0),
fpwobble_slider (0),
fpreset_box (0),
fpreset_button (0),
fpproj_style_box (0),
fporthogonal_button (0),
fpperspective_button (0),
fpfov_text (0),
fpprint_top (0),
fpprint_box (0),
fpprint_col_box (0),
fpprint_style_box (0),
fpprint_text (0),
fpprint_button (0),
fpprint_line (0),
fpprint_col_radio1 (0),
fpprint_col_radio2 (0),
fpprint_style_radio1 (0),
fpprint_style_radio2 (0)
{
  GetXmConnection ();
  ResetView();
  if (fViewId < 0) return;
}


void G4OpenGLXmViewer::UpdateControlPanel () {

  // set new values

  if (fprotation_slider) {
    fprotation_slider->SetInitialValue(fRot_sens);
    fprotation_slider->SetMaxValue(rot_sens_limit);
    fprotation_slider->SetMinValue(0);
  }
  if (fppanning_slider) {
    fppanning_slider->SetInitialValue(fPan_sens);
    fppanning_slider->SetMaxValue(pan_sens_limit);
    fppanning_slider->SetMinValue(0);
  }
  if (fpzoom_slider) {
    fpzoom_slider->SetInitialValue(fVP.GetZoomFactor());
    fpzoom_slider->SetMinValue(zoom_low);
    fpzoom_slider->SetMaxValue(zoom_high);
  }
  if (fpdolly_slider) {
    fpdolly_slider->SetInitialValue(fVP.GetDolly());
    fpdolly_slider->SetMinValue(dolly_low);
    fpdolly_slider->SetMaxValue(dolly_high);
  }

  if (fpwobble_slider) {
    fpwobble_slider->SetInitialValue(fVP.GetDolly());
  }

  if (fppan_set) {
    fppan_set->SetValue(pan_sens_limit);
  }

  if (fprot_set) {
    fprot_set->SetValue(rot_sens_limit);
  }

  if (fpzoom_upper) {
    fpzoom_upper->SetValue(zoom_high);
  }

  if (fpzoom_lower) {
    fpzoom_lower->SetValue(zoom_low);
  }
  if (fpdolly_upper) {
    fpdolly_upper->SetValue(dolly_high);
  }

  if (fpdolly_lower) {
    fpdolly_lower->SetValue(dolly_low);
  }


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
