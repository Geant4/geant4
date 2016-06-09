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
// $Id: G4RTSimpleScanner.cc,v 1.1 2003/09/18 11:14:04 johna Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//

#include "G4RTSimpleScanner.hh"

G4RTSimpleScanner::G4RTSimpleScanner():
  theNRow(0), theNColumn(0), theIRow(0), theIColumn(0) {}

void G4RTSimpleScanner::Initialize(G4int nRow, G4int nColumn) {
  theNRow = nRow;
  theNColumn = nColumn;
  theIRow = 0;
  theIColumn = 0;
}

G4bool G4RTSimpleScanner::Coords(G4int& iRow, G4int& iColumn) {
  if (theIRow >= theNRow) return false;
  iRow = theIRow;
  iColumn = theIColumn;
  if (++theIColumn >= theNColumn) {
    theIColumn = 0;
    ++theIRow;
  }
  return true;
}
