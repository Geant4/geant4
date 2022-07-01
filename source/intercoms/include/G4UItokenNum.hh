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
// G4UItokenNum
//
// Description:
//
// Namespace with enumerator of tokens

// Author: Makoto Asai, 1998
// --------------------------------------------------------------------
#ifndef G4UItokenNum_hh
#define G4UItokenNum_hh 1

#include "globals.hh"

namespace G4UItokenNum
{
  enum tokenNum
  {
    NONE        = 0,
    IDENTIFIER  = 257,
    CONSTINT    = 258,
    CONSTDOUBLE = 259,
    CONSTCHAR   = 260,
    CONSTSTRING = 261,
    GT          = 262,
    GE          = 263,
    LT          = 264,
    LE          = 265,
    EQ          = 266,
    NE          = 267,
    // LOGICALNOT = 268,
    CONSTLONG   = 268,
    LOGICALOR  = 269,
    LOGICALAND = 270,
    SCAREAMER  = 33,
    LPAREN     = 40,
    PLUS       = 43,
    MINUS      = 45
  };

  using yystype = struct yystype
  {
    tokenNum type{ tokenNum::NONE };
    G4double D{ 0.0 };
    G4int I{ 0 };
    G4long L{ 0 };
    char C{ ' ' };
    G4String S;

    yystype()
      : S("")
    {}
    G4bool operator==(const yystype& right) const { return this == &right; }
    yystype& operator=(const yystype& right)
    {
      if(&right == this)
      {
        return *this;
      }
      type = right.type;
      D    = right.D;
      I    = right.I;
      L    = right.L;
      C    = right.C;
      S    = right.S;
      return *this;
    }
    yystype(const yystype& right) { *this = right; }
  };
}  // namespace G4UItokenNum

#endif
