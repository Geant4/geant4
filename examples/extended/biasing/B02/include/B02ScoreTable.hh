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
/// \file biasing/B02/include/B02ScoreTable.hh
/// \brief Definition of the B02ScoreTable class
//
// $Id$
//
// ----------------------------------------------------------------------
// Class B02ScoreTable
//
// Class description:
// 
// This class creates a table from a given G4CellScoreStore
// and prints it to a given output stream.
//

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
//

#ifndef B02ScoreTable_hh
#define B02ScoreTable_hh B02ScoreTable_hh

#include <iostream>
#include <iomanip>
#include <vector>

#include "globals.hh"
#include "G4CellScorerStore.hh"
#include "G4CellScoreValues.hh"
class G4VIStore;


class B02ScoreTable
{
public: // with description
  explicit B02ScoreTable(const G4VIStore *aIStore = 0);

  ~B02ScoreTable();

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
  B02ScoreTable(const B02ScoreTable &);
  B02ScoreTable &operator=(const B02ScoreTable &);

  G4String CreateName(const G4GeometryCell &gCell);
  void PrintLine(const G4String &name,
                 const G4CellScoreValues &sc_scores,
                 std::ostream *out);
  const G4VIStore *fIStore;
  G4int FieldName;
  G4int FieldValue;
};

#endif
