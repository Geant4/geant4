// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RayXView.cc,v 2.2 1998/11/06 13:41:58 allison Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// Nikos Savvas  1st June 1997
// GEANT4 Ray Tracing view.

#ifdef G4VIS_BUILD_RAYX_DRIVER

#include "G4RayXView.hh"

#include "G4RayScene.hh"

G4RayXView::G4RayXView (G4RayScene& scene, const G4String& name):
G4RayView (scene),
G4OwnColormapXView (scene),
G4VView (scene, scene.IncrementViewCount (), name) {

  if (fViewId < 0) return;  // Traps possible error in base class creation.

  if (fTrueColor) {
    fMaxImageMapEntries = 256;  // make a 256 colour map anyway.
  }

  if (fMapEntries != 256 ) {
    G4cerr << "G4RayXView::G4RayXView: no of colormap entries is " << ncells
	 << "\nWe cannot deal with this case for ray tracing yet."
	 << endl;
    fViewId = -1;
    return;
  }

  /************************************
  // Fill 3x3x2 colormap.
  // (In future, I'll write something to evolve the colormap on the fly.
  if (fMapEntries != 256 ) {
    G4cerr << "G4RayXView::G4RayXView: no of colormap entries is " << ncells
	 << "\nWe cannot deal with this case for ray tracing yet."
	 << endl;
    fViewId = -1;
    return;
  }
  XColor* cells = new XColor [fMapEntries];
  red_max = 7;
  green_max = 7;
  blue_max = 4;
  int red_max_1   = red_max   + 1;
  int green_max_1 = green_max + 1;
  int blue_max_1  = blue_max  + 1;
  float red_scale   = 65535 * (float (red_max_1)   / red_max);
  float green_scale = 65535 * (float (green_max_1) / green_max);
  float blue_scale  = 65535 * (float (blue_max_1)  / blue_max);
  red_mult   = 32;
  green_mult = 4;
  blue_mult  = 1;
  for (int iPixel  = 0; iPixel < ncells; iPixel++) {
    cells [iPixel].pixel = iPixel;
    cells [iPixel].red   = (short unsigned int)
      (red_scale   * ((iPixel / red_mult)   % red_max_1));
    cells [iPixel].green = (short unsigned int)
      (green_scale * ((iPixel / green_mult) % green_max_1));
    cells [iPixel].blue = (short unsigned int)
      (blue_scale  * ((iPixel / blue_mult)  % blue_max_1));
    cells [iPixel].flags = DoRed | DoGreen | DoBlue;
  }
  XStoreColors (fpDisplay, fColormap, cells, ncells);
  ************************************/
}

void G4RayXView::DrawView () {

  // If view parameters have changed, start again from scratch.
  // But if not, be smart, and start from private map, depending on
  //   the value of resolution.

  SetView ();

  // Decide whether to initialise image...
  if (fImageViewParameters != fVP         ||
      fMaxImageMapEntries  != fMapEntries ||
      fImageWidth          != fWidth      ||
      fImageHeight         != fHeight) {  // Clear and reset image stores.
    G4cout
      << "G4RayXView::DrawView: image re-initialised because...\n  ";
    if (fImageViewParameters != fVP)
      G4cout << "view parameters, ";
    if (fMaxImageMapEntries  != fMapEntries)
      G4cout << "no. of entries in colormap, ";
    if (fImageWidth          != fWidth)
      G4cout << "width, ";
    if (fImageHeight         != fHeight)
      G4cout << "  height, ";
    G4cout << "changed." << endl;
    // Reset above parameters...
    fImageViewParameters = fVP;
    fMaxImageMapEntries = fMapEntries;
    fImageWidth = fWidth;
    fImageHeight = fHeight;
    // Now things that need initialising...
    InitialiseImage ();
  }
  // Now things that need resetting for every Draw.
  TriggerRedraw ();

  XGCValues xgcValues;
  XColor cell;
  XColor* cells = new XColor [fMapEntries];

  int i, j, k, iPixVal;
  float red, green, blue;
  // Square of side k whose upper left corner is (i, j).  The above
  // are set by the function Scan, inherited from G4RayView.  You
  // don't have to use it - you can call Trace, also inherited from
  // G4RayView, directly.  iPixVal is the pixel value for this colour,
  // i.e., the key in fRayHitColourIndex.
  while (Scan (i, j, k, red, green, blue, iPixVal)) {

    if(fTrueColor) {
      unsigned long int trueColorPixval =
	((unsigned long int) (red   * fRedMax))   << fRedShift   |
	((unsigned long int) (green * fGreenMax)) << fGreenShift |
	((unsigned long int) (blue  * fBlueMax))  << fBlueShift;
      xgcValues.foreground = trueColorPixval;
      XChangeGC (fpDisplay, fGC, GCForeground, &xgcValues);
    }

    else {

      G4int nCells = fRayHitColourMap.entries ();
      if (nCells > fMapEntries) {
	G4Exception ("G4RayXView::DrawView: fRayHitColourMap overflow!!!");
      }

      if (fImageMapRejigged && ++fWaitCount > fWaitLimit) {
	// Note: ++fWaitCount only executed if fImageMapRejigged is true.
	fWaitCount = 0;
	fImageMapRejigged = false;
	fTriggerRescan = true;
	//G4cout << "G4RayXView::DrawView: remaking whole colormap." << endl;
	for (int iCells = 0; iCells < nCells; iCells++) {
	  cells [iCells].pixel = iCells;
	  cells [iCells].red   = (short unsigned int)
	    (65535 * fRayHitColourMap (iCells).red);
	  cells [iCells].green = (short unsigned int)
	    (65535 * fRayHitColourMap (iCells).green);
	  cells [iCells].blue  = (short unsigned int)
	    (65535 * fRayHitColourMap (iCells).blue);
	  cells [iCells].flags = DoRed | DoGreen | DoBlue;
	}
	XStoreColors (fpDisplay, fColormap, cells, nCells);
      }
      else {
	if (iPixVal == nCells - 1) {
	  // Just do the last one, the latest addition.
	  cell.pixel = iPixVal;
	  cell.red   = (short unsigned int)
	    (65535 * fRayHitColourMap (iPixVal).red);
	  cell.green = (short unsigned int)
	    (65535 * fRayHitColourMap (iPixVal).green);
	  cell.blue  = (short unsigned int)
	    (65535 * fRayHitColourMap (iPixVal).blue);
	  cell.flags = DoRed | DoGreen | DoBlue;
	  XStoreColor (fpDisplay, fColormap, &cell);
	}
	// (Otherwise assume an old pixel value is being used,
	// i.e., do nothing.)
      }
      
      XSetForeground (fpDisplay, fGC, iPixVal);
    }

    if (k <= 1) {
      XDrawPoint (fpDisplay, fWindow, fGC, i, j);
    }
    else {
      XFillRectangle (fpDisplay, fWindow, fGC, i, j, k, k);
    }

    XFlush (fpDisplay);

    //for (int ii = 0; ii < 1000000; ii++);  //?? 1000000 = about 0.5 secs!!
  }

  delete cells;

  FinishView ();
}

#endif
