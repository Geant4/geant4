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
// $Id: G4VRTScanner.hh 66373 2012-12-18 09:41:34Z gcosmo $
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

  G4VRTScanner();
  virtual ~G4VRTScanner();
  
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

inline
void G4VRTScanner::Draw(unsigned char, unsigned char, unsigned char) {}

#endif
