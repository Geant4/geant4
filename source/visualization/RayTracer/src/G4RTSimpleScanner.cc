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
// $Id: G4RTSimpleScanner.cc 66373 2012-12-18 09:41:34Z gcosmo $
//
//

#include "G4RTSimpleScanner.hh"

G4RTSimpleScanner::G4RTSimpleScanner():
  G4VRTScanner(), theNRow(0), theNColumn(0), theIRow(0), theIColumn(0)
{
  theGSName = "RayTracer";
  theGSNickname = "RayTracer";
}

G4RTSimpleScanner::~G4RTSimpleScanner(){}

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
