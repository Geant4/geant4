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
// $Id: G4RayTracerXViewer.cc,v 1.1 2005-07-17 13:59:24 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4RayTracerXViewer.hh"

#include "G4VSceneHandler.hh"
#include "G4Scene.hh"
#include "G4RayTracerX.hh"
#include "G4UImanager.hh"

#include <cstdlib>

//#include <X11/Xlib.h>
#include <X11/Xutil.h>
//#include <X11/Xos.h>
#include <X11/Xatom.h>

extern "C" {
  Bool G4RayTracerXViewerWaitForNotify (Display*, XEvent* e, char* arg) {
    return (e->type == MapNotify) && (e->xmap.window == (Window) arg);
  }
}

G4RayTracerXViewer::G4RayTracerXViewer
(G4VSceneHandler& sceneHandler, const G4String& name):
  G4RayTracerViewer(sceneHandler, name)
{
  G4RayTracerX* theTracer = 
    (G4RayTracerX*) fSceneHandler.GetGraphicsSystem();

  Display*& display = theTracer->display;
  Window& win = theTracer->win;
  GC& gc = theTracer->gc;
  XStandardColormap*& scmap = theTracer->scmap;

  display = XOpenDisplay(0);  // Use display defined by DISPLAY environment.
  if (!display) {
    fViewId = -1;  // This flags an error.
    G4cerr << "G4RayTracerXViewer::G4RayTracerXViewer: cannot get display."
	   << G4endl;
    return;
  }

  int screen_num = DefaultScreen(display);

  win = XCreateSimpleWindow
    (display, RootWindow(display, screen_num),
     0, 0,                              // Corner position.
     fVP.GetWindowSizeHintX(),          // Width.
     fVP.GetWindowSizeHintY(),          // Height.
     0,                                 // Border width.
     WhitePixel(display, screen_num),   // Border colour.
     BlackPixel(display, screen_num));  // Background colour.

  XGCValues values;
  gc = XCreateGC(display, win, 0, &values);

  int nMaps;
  Status status = XGetRGBColormaps
    (display, RootWindow(display, screen_num),
     &scmap, &nMaps, XA_RGB_BEST_MAP);
  if (!status) {
    system("xstdcmap -best");  // ...and try again...
    Status status = XGetRGBColormaps
      (display, RootWindow(display, screen_num),
       &scmap, &nMaps, XA_RGB_BEST_MAP);
    if (!status) {
      fViewId = -1;  // This flags an error.
      G4cerr <<
	"G4RayTracerXViewer::G4RayTracerXViewer: cannot get color map."
	"\n  Perhaps your system does not support RGB_BEST_MAP."
	     << G4endl;
      return;
    }
  }
  if (!scmap->colormap) {
    fViewId = -1;  // This flags an error.
    G4cerr << "G4RayTracerXViewer::G4RayTracerXViewer: color map empty."
	   << G4endl;
    return;
  }

  XSizeHints* size_hints = XAllocSizeHints();
  XWMHints* wm_hints = XAllocWMHints();
  XClassHint* class_hint = XAllocClassHint();
  const char* window_name = fName.c_str();
  XTextProperty windowName;
  XStringListToTextProperty((char**)&window_name, 1, &windowName);
  XSetWMProperties(display, win, &windowName, &windowName,
		   0, 0, size_hints, wm_hints, class_hint);

  XMapWindow(display, win);

  // Wait for window to appear (wait for an "map notify" event).
  XSelectInput(display, win, StructureNotifyMask);
  XEvent event;
  XIfEvent (display, &event, G4RayTracerXViewerWaitForNotify, (char*) win);

}

G4RayTracerXViewer::~G4RayTracerXViewer() {}
