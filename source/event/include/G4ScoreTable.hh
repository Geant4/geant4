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
// $Id: G4ScoreTable.hh,v 1.2 2003/06/16 16:50:29 gunter Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// ----------------------------------------------------------------------
// Class G4ScoreTable
//
// Class description:
// 
// This class creates a table from a given G4CellScoreStore
// and prints it to a given output stream.
//

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
//

#ifndef G4ScoreTable_hh
#define G4ScoreTable_hh G4ScoreTable_hh

#include <iostream>
#include <iomanip>
#include <vector>

#include "globals.hh"
#include "G4CellScorerStore.hh"
#include "G4CellScoreValues.hh"
class G4VIStore;


class G4ScoreTable
{
public: // with description
  explicit G4ScoreTable(const G4VIStore *aIStore = 0);

  ~G4ScoreTable();

  void Print(const G4MapGeometryCellCellScorer &cs, 
	     std::ostream *out = 0);
    // create the table and print it to the ouput stream.

  void PrintHeader(std::ostream *out);
    // print the table header, done by the above Print() function

  void PrintTable(const G4MapGeometryCellCellScorer &cs,
		  std::ostream *out);
    // print the table contend, done by the above Print() function

  std::string FillString(const std::string &name, 
			   char c, G4int n, G4bool back = true);
    // create a string of length n, by filling up
    // name with the char c.

private:
  G4ScoreTable(const G4ScoreTable &);
  G4ScoreTable &operator=(const G4ScoreTable &);

  G4String CreateName(const G4GeometryCell &gCell);
  void PrintLine(const G4String &name,
		 const G4CellScoreValues &sc_scores,
		 std::ostream *out);
  const G4VIStore *fIStore;
  G4int FieldName;
  G4int FieldValue;
};

#endif
