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
// $Id: G4UIArrayString.hh 66892 2013-01-17 10:57:59Z gunter $
//

#ifndef G4UIArrayString_h
#define G4UIArrayString_h

#include "globals.hh"

//
//   Description: 
//   This class is a utility class for listing commands like 
//   UNIX "ls" command.
//

class G4UIArrayString {
private:
  G4String* stringArray;
  G4int nElement;
  G4int nColumn;

  G4String* GetElement(G4int icol, G4int irow) const;
  // #rows in a column
  G4int GetNRow(G4int icol) const;
  // max field width of a column
  G4int GetNField(G4int icol) const;
  // calculate total column width
  G4int CalculateColumnWidth() const;

public:
  G4UIArrayString(const G4String& stream);
  ~G4UIArrayString();

  void Show(G4int ncol);
};

#endif

