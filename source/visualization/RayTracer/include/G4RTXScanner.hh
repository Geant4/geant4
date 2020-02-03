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

#ifdef G4VIS_BUILD_RAYTRACERX_DRIVER

#ifndef G4RTXScanner_H
#define G4RTXScanner_H 1

// class description:
//
// G4RTXScanner
// Provides a sequence of window coordinates suitable for a visible
// window of ever increasing resolution.

#include "G4VRTScanner.hh"

#include <X11/Xlib.h>
#include <X11/Xutil.h>

class G4ViewParameters;

class G4RTXScanner: public G4VRTScanner {

public: // with description

  G4RTXScanner();
  virtual ~G4RTXScanner();

  // Compiler defaults for copy constructor and assignmemt.

  virtual const G4String& GetGSName() const;
  // Get name that acts as graphics system name.

  virtual const G4String& GetGSNickname() const;
  // Get name that acts as graphics system nickname.  It is this that
  // the user specifies on the /vis/open and /vis/sceneHandler/create
  // commands.

  virtual void Initialize(G4int nRow, G4int nColumn);
  // Intialises scanner for window with nRow rows and nColumn columns.

  virtual G4bool Coords(G4int& iRow, G4int& iColumn);
  // Supplies coordinate (iRow,iColumn) and returns false when the
  // sequence has finished, i.e., on the call *after* suplying the
  // last valid coordinate.

  virtual void Draw
  (unsigned char red, unsigned char green, unsigned char blue);
  // Draw coloured square at current position.

  G4bool GetXWindow(const G4String& name, G4ViewParameters&);

protected:
  G4String theGSName, theGSNickname;
  G4int theNRow, theNColumn, theStep, theIRow, theIColumn;
  // X Window variables...
  Display* display;
  Window win;
  GC gc;
  XStandardColormap *scmap;
};

#endif

#endif
