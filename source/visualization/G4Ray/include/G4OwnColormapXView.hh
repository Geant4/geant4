// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OwnColormapXView.hh,v 2.1 1998/10/23 17:58:21 allison Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// John Allison  1st June 1997
// Loosely based on G4OpenGLXView by Andrew Walkden  7th February 1997
// G4OwnColormapXView : Class to provide X Windows with its own colormap
//   for a GEANT4 view.  The derived view must define the colormap.

#ifdef G4VIS_BUILD_RAYX_DRIVER

#ifndef G4OWNCOLORMAPXVIEW_HH
#define G4OWNCOLORMAPXVIEW_HH

#include "G4VView.hh"

#include <X11/Xlib.h>
#include <X11/Intrinsic.h>

class G4VScene;

class G4OwnColormapXView: virtual public G4VView {

public:
  G4OwnColormapXView (G4VScene& scene);
  ~G4OwnColormapXView ();
  void ClearView ();
  void FinishView ();

protected:
  char*                fcharViewName;
  Display*             fpDisplay;
  int                  fScreen;
  int                  fDepth;
  G4String             fVisualName;
  G4bool               fTrueColor;
  Visual*              fpVisual;
  int                  fMapEntries;
  Window               fRootW;
  int                  fRedShift;
  int                  fGreenShift;
  int                  fBlueShift;
  float                fRedMax;
  float                fGreenMax;
  float                fBlueMax;
  Colormap             fColormap;
  XWindowAttributes    fRootWA;
  GC                   fGC;
  XSetWindowAttributes fSWA;
  Window               fWindow;
  XTextProperty        fWindowName;
  XTextProperty        fIconName;
  XSizeHints*          fSize_hints;
  XWMHints*            fWM_hints;
  XClassHint*          fClass_hints;
  XWindowAttributes    fWA;
  int                  fWidth;
  int                  fHeight;
  Pixmap               fIcon_pixmap;
  XEvent               fEvent;
};

#endif

#endif
