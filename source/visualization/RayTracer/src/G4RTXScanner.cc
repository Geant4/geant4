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
// $Id: G4RTXScanner.cc,v 1.2 2005/07/20 20:39:02 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//

#ifdef G4VIS_BUILD_RAYTRACERX_DRIVER

#include "G4RTXScanner.hh"

#include "G4RayTracerX.hh"

#include <X11/Xlib.h>
#include <X11/Xutil.h>

G4RTXScanner::G4RTXScanner():
  theGSName("RayTracerX"), theGSNickname("RayTracerX"),
  theNRow(0), theNColumn(0), theIRow(0), theIColumn(0) {}

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

void G4RTXScanner::Draw
(unsigned char red, unsigned char green, unsigned char blue,
 G4RayTracer* aTracer)
// Draw coloured square at current position.
{
  G4RayTracerX* theTracer = (G4RayTracerX*) aTracer;
  Display*& display = theTracer->display;
  Window& win = theTracer->win;
  GC& gc = theTracer->gc;
  XStandardColormap*& scmap = theTracer->scmap;

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
