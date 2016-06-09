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
// $Id: G4VRTScanner.hh,v 1.1 2003/09/18 11:14:03 johna Exp $
// GEANT4 tag $Name: geant4-06-00 $
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

  virtual void Initialize(G4int nRow, G4int nColumn) = 0;
  // Intialises scanner for window with nRow rows and nColumn columns.

  virtual G4bool Coords(G4int& iRow, G4int& iColumn) = 0;
  // Supplies coordinate (iRow,iColumn) and returns false when the
  // sequence has finished, i.e., on the call *after* suplying the
  // last valid coordinate.

};

#endif
