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
// $Id: G4StandardScoreTable.hh,v 1.1 2002-07-10 15:51:04 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4StandardScoreTable
//
// Class description:
// 
// This class creates a table from a given G4MapPtkStandardCellScorer
// and prints it to a given output stream.
//
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
//

#ifndef G4StandardScoreTable_hh
#define G4StandardScoreTable_hh G4StandardScoreTable_hh

#include "g4std/iostream"
#include "g4std/iomanip"
#include "g4std/vector"

#include "globals.hh"
#include "G4MapPtkStandardCellScorer.hh"

class G4VIStore;


class G4StandardScoreTable
{
public:
  G4StandardScoreTable(const G4VIStore *aIStore = 0);
  ~G4StandardScoreTable(){}
  void Print(const G4MapPtkStandardCellScorer &mptkscorer, 
	     G4std::ostream *out);
    // create the table and print it to the ouput stream.

  void PrintHeader(G4std::ostream *out);
    // print the table header, done by the above Print() function

  void PrintTable(const G4MapPtkStandardCellScorer &mptkscorer,
		  G4std::ostream *out);
    // print the table contend, done by the above Print() function

  G4std::string FillString(const G4std::string &name, 
			   char c, G4int n, G4bool back = true);
    // create a string of length n, by filling up
    // name with the char c.

private:
  G4String CreateName(G4PTouchableKey ptk);
  void PrintLine(G4String &name,
		 G4StandardCellScoreValues sc_scores,
		 G4std::ostream *out);
  const G4VIStore *fIStore;
  G4int FieldName;
  G4int FieldValue;
};

#endif
