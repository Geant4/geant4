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
// $Id: G4UItokenNum.hh 67965 2013-03-13 09:35:29Z gcosmo $
//
// G4UItokenNum.hh

#ifndef G4UItokenNum_hh
#define G4UItokenNum_hh 1
#include "globals.hh"


enum  tokenNum
{
  NONE       = 0,
  IDENTIFIER = 257,
  CONSTINT   = 258,
  CONSTDOUBLE= 259,
  CONSTCHAR  = 260,
  CONSTSTRING= 261,
  GT         = 262,
  GE         = 263,
  LT         = 264,
  LE         = 265,
  EQ         = 266,
  NE         = 267,
  //LOGICALNOT = 268,
  LOGICALOR  = 269,
  LOGICALAND = 270,
  SCAREAMER  = 33,
  LPAREN     = 40,
  PLUS       = 43,
  MINUS      = 45
};


typedef struct yystype
{
    tokenNum type;
    G4double D;
    G4int    I;
    char     C;
    G4String S;

    yystype() : type(NONE), D(0.0), I(0), C(' '), S("")
    {
    }
    G4int operator==(const yystype& right) const
    {
      return (this == &right)?1:0;
    }
    yystype& operator=(const yystype& right)
    {
      if (&right==this) return *this;
      type = right.type;
      D = right.D;
      I = right.I;
      C = right.C;
      S = right.S;
      return *this;
    }
    yystype(const yystype& right)
    {
      *this=right;
    }
} yystype;

#endif
