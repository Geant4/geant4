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
// $Id: G4VRTScanner.hh,v 1.4 2006/01/11 18:01:33 allison Exp $
// GEANT4 tag $Name: geant4-08-00-patch-01 $
//
//

#ifndef G4VRTScanner_H
#define G4VRTScanner_H 1

// class description:
//
// G4VRTScanner
// Interface class for provider of a sequence of window coordinates.

#include "globals.hh"

class G4VRTScanner {

public: // with description

  virtual const G4String& GetGSName() const = 0;
  // Get name that acts as graphics system name.

  virtual const G4String& GetGSNickname() const = 0;
  // Get name that acts as graphics system nickname.  It is this that
  // the user specifies on the /vis/open and /vis/sceneHandler/create
  // commands.

  virtual void Initialize(G4int nRow, G4int nColumn) = 0;
  // Intialises scanner for window with nRow rows and nColumn columns.

  virtual G4bool Coords(G4int& iRow, G4int& iColumn) = 0;
  // Supplies coordinate (iRow,iColumn) and returns false when the
  // sequence has finished, i.e., on the call *after* suplying the
  // last valid coordinate.

  virtual void Draw
  (unsigned char red, unsigned char green, unsigned char blue);
  // Draw coloured square at current position.

};

inline void G4VRTScanner::Draw
(unsigned char, unsigned char, unsigned char)
{}

#endif
