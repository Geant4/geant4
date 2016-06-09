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
// $Id: G4RTSimpleScanner.cc,v 1.2 2005/07/17 13:59:24 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//

#include "G4RTSimpleScanner.hh"

G4RTSimpleScanner::G4RTSimpleScanner():
  theGSName("RayTracer"), theGSNickname("RayTracer"),
  theNRow(0), theNColumn(0), theIRow(0), theIColumn(0) {}

const G4String& G4RTSimpleScanner::GetGSName() const
{return theGSName;}

const G4String& G4RTSimpleScanner::GetGSNickname() const
{return theGSNickname;}

void G4RTSimpleScanner::Initialize(G4int nRow, G4int nColumn) {
  theNRow = nRow;
  theNColumn = nColumn;
  theIRow = 0;
  theIColumn = -1;
}

G4bool G4RTSimpleScanner::Coords(G4int& iRow, G4int& iColumn)
{
  // Increment column and, if necessary, increment row...
  ++theIColumn;
  if (theIColumn >= theNColumn) {
    theIColumn = 0;
    ++theIRow;
  }

  // Return if finished...
  if (theIRow >= theNRow) return false;

  // Return current row and column...
  iRow = theIRow;
  iColumn = theIColumn;
  return true;
}
