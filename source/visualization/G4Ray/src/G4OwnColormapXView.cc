// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OwnColormapXView.cc,v 2.4 1998/10/23 17:58:22 allison Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// John Allison  1st June 1997
// Loosely based on G4OpenGLXView by Andrew Walkden  7th February 1997
// G4OwnColormapXView : Class to provide X Windows with its own colormap
//   for a GEANT4 view.  The derived view must define the colormap.

#ifdef G4VIS_BUILD_RAYX_DRIVER

#include "G4OwnColormapXView.hh"

#include <X11/Xatom.h>
#include <X11/Xutil.h>
#include <iomanip.h>

#include <stdlib.h>
#include <string.h>

#define NewString(str) \
 ((str) != NULL ? (strcpy((char*)malloc((unsigned)strlen(str) + 1), str)) : (char*)NULL)

static Bool WaitForNotify (Display* d, XEvent* e, char* arg) {
  return (e->type == MapNotify) && (e->xmap.window == (Window) arg);
}

G4OwnColormapXView::G4OwnColormapXView (G4VScene& scene):
G4VView (scene, -1),
fRedShift (0),
fGreenShift (0),
fBlueShift (0),
fRedMax (1),
fGreenMax (1),
fBlueMax (1)
{
  fpDisplay = XOpenDisplay (0);
  if (!fpDisplay) {
    G4cerr << "G4OwnColormapXView::G4OwnColormapXView couldn't open display." << endl;
    fViewId = -1;  // This flags an error.
    return;
  }

  fScreen   = XDefaultScreen   (fpDisplay);

  fDepth    = XDefaultDepth    (fpDisplay, fScreen);
  if (fDepth < 8) {
    G4cerr << "G4OwnColormapXView::G4OwnColormapXView: too few pixel"
      "\nplanes - need at least 8 for own colormap." << endl;
    fViewId = -1;  // This flags an error.
    return;
  }

  XVisualInfo vinfo;
  if (XMatchVisualInfo (fpDisplay, fScreen, fDepth, PseudoColor, &vinfo)) {
    G4cout << "G4OwnColormapXView::G4OwnColormapXView: PseudoColor." << endl;
    fVisualName = "PseudoColor";
    fTrueColor = false;
  }
  else if
    (XMatchVisualInfo (fpDisplay, fScreen, fDepth, DirectColor, &vinfo)) {
    G4cout << "G4OwnColormapXView::G4OwnColormapXView: DirectColor." << endl;
    fVisualName = "DirectColor";
    fTrueColor = false;
  }
  else if
    (XMatchVisualInfo (fpDisplay, fScreen, fDepth, TrueColor, &vinfo)) {
    G4cout << "G4OwnColormapXView::G4OwnColormapXView: TrueColor." << endl;
    fVisualName = "TrueColor";
    fTrueColor = true;
  }
  else {
    G4cerr << "G4OwnColormapXView::G4OwnColormapXView: neither a PseudoColor nor"
      "\na DirectColor nor a TrueColor visual -  needs to be one of these"
      "\nfor own colormap."
	 << endl;
    fViewId = -1;  // This flags an error.
    return;
  }
 
  fpVisual = vinfo.visual;
  if (!fpVisual) {
    G4cerr << "G4OwnColormapXView::G4OwnColormapXView couldn't get a visual."
	 << endl;
    fViewId = -1;  // This flags an error.
    return;
  }

  fMapEntries = fpVisual -> map_entries;

  fRootW = XRootWindow (fpDisplay, fScreen);

  if (fTrueColor) {
    //  Set up parameters for masks etc.
    unsigned long oneBit;
    oneBit = 1;
    while (!(fpVisual->red_mask & oneBit)) {oneBit <<= 1; fRedShift++;}
    while   (fpVisual->red_mask & oneBit)  {oneBit <<= 1; fRedMax *= 2.;}
    fRedMax -= 1.;
    oneBit = 1;
    while (!(fpVisual->green_mask & oneBit)) {oneBit <<= 1; fGreenShift++;}
    while   (fpVisual->green_mask & oneBit)  {oneBit <<= 1; fGreenMax *= 2.;}
    fGreenMax -= 1.;
    oneBit = 1;
    while (!(fpVisual->blue_mask & oneBit)) {oneBit <<= 1; fBlueShift++;}
    while   (fpVisual->blue_mask & oneBit)  {oneBit <<= 1; fBlueMax *= 2.;}
    fBlueMax -= 1.;
  }
  else {
    // Own colour map.
    fColormap = XCreateColormap (fpDisplay, fRootW, fpVisual, AllocAll);
    if (!fColormap || fColormap == XDefaultColormap (fpDisplay, fScreen)) {
      G4cerr << "G4OwnColormapXView::G4OwnColormapXView failed to create"
	"\na own Colormap." << endl;
      fViewId = -1;  // This flags an error.
      return;
    }
  }

  fGC = XDefaultGC (fpDisplay, fScreen);    // Graphics Context for XDraw...

  XGetWindowAttributes (fpDisplay, fRootW, &fRootWA);  // Get root window atts.

  //Choose a preferred window size.
  G4int requestWinSize_x;
  G4int requestWinSize_y;
  if (fRootWA.width > fRootWA.height) {
    requestWinSize_x = (fRootWA.height)/2;
    requestWinSize_y = (fRootWA.height)/2;
  }
  else {
    requestWinSize_x = (fRootWA.width)/2;
    requestWinSize_y = (fRootWA.width)/2;
  }
  if (requestWinSize_x < fVP.GetWindowSizeHintX ())
    requestWinSize_x = fVP.GetWindowSizeHintX ();
  if (requestWinSize_y < fVP.GetWindowSizeHintY ())
    requestWinSize_y = fVP.GetWindowSizeHintY ();

  // Choose an origin (the window manager will ignore it!)
  G4int x_origin;
  G4int y_origin;
  x_origin = fRootWA.x;
  y_origin = fRootWA.y;

  if (!fTrueColor) {
    fSWA.colormap = fColormap;
  }
  fSWA.background_pixel = 0;
  fSWA.border_pixel = 0;
  fSWA.event_mask = ExposureMask | ButtonPressMask | StructureNotifyMask;
  fSWA.backing_store = WhenMapped;
  unsigned int mask = CWBackPixel | CWBorderPixel |
    CWEventMask | CWBackingStore;
  if (!fTrueColor) mask = mask | CWColormap;
  fWindow = XCreateWindow (fpDisplay, fRootW, x_origin, 
			   y_origin, requestWinSize_x, requestWinSize_y, 0,
			   fDepth, InputOutput, fpVisual,
			   mask,
			   &fSWA);

  G4cout << "Window name: " << fName << endl;
  fcharViewName = new char [100];
  strncpy (fcharViewName, fName, 100);
  XStringListToTextProperty (&fcharViewName, 1, &fWindowName);
  XStringListToTextProperty (&fcharViewName, 1, &fIconName);

  fSize_hints = XAllocSizeHints();
  fSize_hints -> flags = PPosition | PSize | PMinSize;
  fSize_hints -> min_width = 300;
  fSize_hints -> min_height = 300;

  fWM_hints = XAllocWMHints();
  fWM_hints -> initial_state = NormalState;
  fWM_hints -> input = True;
  fWM_hints -> icon_pixmap = fIcon_pixmap;
  fWM_hints -> flags = StateHint | IconPixmapHint | InputHint;

  fClass_hints = XAllocClassHint();
  fClass_hints -> res_name  = NewString("G4SimpleX");
  fClass_hints -> res_class = NewString("G4SimpleX");

  XSetWMProperties (fpDisplay, fWindow, &fWindowName, &fIconName, 0, 0, 
                    fSize_hints, fWM_hints, fClass_hints);
  
// request X to Draw window on screen.
  XMapWindow (fpDisplay, fWindow);

// Wait for window to appear (wait for an "expose" event).
  XIfEvent (fpDisplay, &fEvent, WaitForNotify, (char*) fWindow);

  XGetWindowAttributes (fpDisplay, fWindow, &fWA);
  fWidth  = fWA.width;
  fHeight = fWA.height;
  G4cout << "G4OwnColormapXView instantiated."
    "\n A simple X window of visual type "
       << fVisualName
       << ",\n size " << fWidth << 'x' << fHeight << ", depth " << fDepth;
  if (fTrueColor) {
    G4cout << "\nRed shift  : " << setw(2) << fRedShift
	 << ", red max  : " << fRedMax;
    G4cout << "\nGreen shift: " << setw(2) << fGreenShift
	 << ", green max: " << fGreenMax;
    G4cout << "\nBlue shift : " << setw(2) << fBlueShift
	 << ", blue max : " << fBlueMax;
  }
  else {
    G4cout << ",\n with " << fMapEntries << " colormap entries";
  }
  G4cout << ", has been opened." << endl;
}

G4OwnColormapXView::~G4OwnColormapXView () {
  if (fViewId >= 0) {
    XFree (fSize_hints);
    XFree (fWM_hints);
    if(fClass_hints -> res_name!=NULL)  free(fClass_hints -> res_name);
    if(fClass_hints -> res_class!=NULL) free(fClass_hints -> res_class);
    XFree (fClass_hints);
    //Close a window from here
    XUnmapWindow   (fpDisplay, fWindow);
    XDestroyWindow (fpDisplay, fWindow);
  }
  if (fpDisplay) {
    XFlush         (fpDisplay);
    XCloseDisplay  (fpDisplay);
  }
}

void G4OwnColormapXView::ClearView () {
  XClearWindow (fpDisplay, fWindow);
}

void G4OwnColormapXView::FinishView () {
  XFlush (fpDisplay);
}

#endif
