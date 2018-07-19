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
// $Id: G4RTXScanner.cc 103626 2017-04-19 13:29:18Z gcosmo $
//
//

#ifdef G4VIS_BUILD_RAYTRACERX_DRIVER

#include "G4RTXScanner.hh"

#include "G4TheRayTracer.hh"
#include "G4RayTracerXViewer.hh"
#include "G4ViewParameters.hh"
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>

extern "C" {
  Bool G4RayTracerXScannerWaitForNotify (Display*, XEvent* e, char* arg) {
    return (e->type == MapNotify) && (e->xmap.window == (Window) arg);
  }
}

G4RTXScanner::G4RTXScanner():
  G4VRTScanner(), theNRow(0), theNColumn(0), theStep(0)
  ,theIRow(0), theIColumn(0)
  ,display(0), win(0), scmap(0)
{
  theGSName = "RayTracerX";
  theGSNickname = "RayTracerX";
}

G4RTXScanner::~G4RTXScanner() {}

const G4String& G4RTXScanner::GetGSName() const
{return theGSName;}

const G4String& G4RTXScanner::GetGSNickname() const
{return theGSNickname;}

void G4RTXScanner::Initialize(G4int nRow, G4int nColumn) {
  theNRow = nRow;
  theNColumn = nColumn;
  G4int nMax = std::max (nRow, nColumn);
  theStep = 1;
  if (nMax > 3) {
    for (;;) {
      theStep *= 3;
      if (theStep > nMax) break;
    }
  }
  theIRow = theStep / 2;
  theIColumn = theStep / 2 - theStep;
}

G4bool G4RTXScanner::Coords(G4int& iRow, G4int& iColumn)
{
  // Increment column...
  theIColumn += theStep;

  // Skip coordinates covered in the previous scan...
  if ((theIColumn + (3 * theStep) / 2 + 1)%(3 * theStep) == 0 &&
      (theIRow + (3 * theStep) / 2 + 1)%(3 * theStep) == 0)
    theIColumn += theStep;

  //  If necessary, increment row...
  if (theIColumn >= theNColumn) {
    theIColumn = theStep / 2;
    theIRow += theStep;
  }

  // Return if finished...
  if (theIRow >= theNRow && theStep <= 1) return false;

  // Start next scan if necessary...
  if (theIRow >= theNRow) {
    theStep /= 3;
    theIRow = theStep / 2;
    theIColumn = theStep / 2;
  }

  // Return current row and column...
  iRow = theIRow;
  iColumn = theIColumn;
  return true;
}

G4bool G4RTXScanner::GetXWindow(const G4String& name, G4ViewParameters& vp)
{
  display = XOpenDisplay(0);  // Use display defined by DISPLAY environment.
  if (!display) {
    G4cerr << "G4RTXScanner::Initialize(): cannot get display."
	   << G4endl;
    return false;
  }

  int screen_num = DefaultScreen(display);

  // Window size and position...
  int xOffset = 0, yOffset = 0;
  XSizeHints* size_hints = XAllocSizeHints();
  unsigned int width, height;
  const G4String& XGeometryString = vp.GetXGeometryString();
  if (!XGeometryString.empty()) {
    G4int geometryResultMask = XParseGeometry
      ((char*)XGeometryString.c_str(),
       &xOffset, &yOffset, &width, &height);
    if (geometryResultMask & (WidthValue | HeightValue)) {
      if (geometryResultMask & XValue) {
	if (geometryResultMask & XNegative) {
	  xOffset = DisplayWidth(display, screen_num) + xOffset - width;
	}
	size_hints->flags |= PPosition;
	size_hints->x = xOffset;
      }
      if (geometryResultMask & YValue) {
	if (geometryResultMask & YNegative) {
	  yOffset = DisplayHeight(display, screen_num) + yOffset - height;
	}
	size_hints->flags |= PPosition;
	size_hints->y = yOffset;
      }
    } else {
      G4cout << "ERROR: Geometry string \""
	     << XGeometryString
	     << "\" invalid.  Using \"600x600\"."
	     << G4endl;
      width = 600;
      height = 600;
    }
  } else {
    G4cout << "ERROR: Geometry string \""
	   << XGeometryString
	   << "\" is empty.  Using \"600x600\"."
	   << G4endl;
    width = 600;
    height = 600;
  }
  size_hints->width = width;
  size_hints->height = height;
  size_hints->flags |= PSize;

  win = XCreateSimpleWindow
    (display, RootWindow(display, screen_num),
     xOffset, yOffset, width, height,
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
    status = XGetRGBColormaps
      (display, RootWindow(display, screen_num),
       &scmap, &nMaps, XA_RGB_BEST_MAP);
    if (!status) {
      G4cerr <<
	"G4RTXScanner::Initialize(): cannot get color map."
	"\n  Perhaps your system does not support XA_RGB_BEST_MAP."
	     << G4endl;
      return false;
    }
  }
  if (!scmap->colormap) {
    G4cerr << "G4RTXScanner::Initialize(): color map empty."
	   << G4endl;
    return false;
  }

  XWMHints* wm_hints = XAllocWMHints();
  XClassHint* class_hint = XAllocClassHint();
  const char* window_name = name.c_str();
  XTextProperty windowName;
  XStringListToTextProperty((char**)&window_name, 1, &windowName);

  XSetWMProperties(display, win, &windowName, &windowName,
		   0, 0, size_hints, wm_hints, class_hint);

  XMapWindow(display, win);

  // Wait for window to appear (wait for an "map notify" event).
  XSelectInput(display, win, StructureNotifyMask);
  XEvent event;
  XIfEvent (display, &event, G4RayTracerXScannerWaitForNotify, (char*) win);

  return true;
}

void G4RTXScanner::Draw
(unsigned char red, unsigned char green, unsigned char blue)
// Draw coloured square at current position.
{
  unsigned long pixel_value = scmap->base_pixel +
    ((unsigned long) ((red * scmap->red_max) / 256.) * scmap->red_mult) +
    ((unsigned long) ((green * scmap->green_max) / 256.) * scmap->green_mult) +
    ((unsigned long) ((blue * scmap->blue_max) / 256.) * scmap->blue_mult);
  XSetForeground(display, gc, pixel_value);

  if (theStep > 1) {
    XFillRectangle(display, win, gc,
		   theIColumn - theStep / 2,
		   theIRow - theStep / 2,
		   theStep, theStep);
  } else {
    XDrawPoint(display, win, gc, theIColumn, theIRow);
  }

  XFlush(display);
}

#endif
