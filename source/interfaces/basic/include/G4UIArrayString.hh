// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIArrayString.hh,v 1.1 2000-03-26 23:02:33 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

