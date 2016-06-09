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
// $Id: G4UItokenNum.hh,v 1.7 2003/06/07 16:40:48 gcosmo Exp $
// GEANT4 tag $Name: geant4-05-02-patch-01 $
//
// G4UItokenNum.hh

#ifndef G4UItokenNum_hh
#define G4UItokenNum_hh 1
#include "globals.hh"


enum  tokenNum
{
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
  LOGICALAND = 270
};


typedef struct yystype
{
    tokenNum      type;
    G4double D;
    G4int    I;
    char     C;
    G4String S;

    yystype() : D(0.0), I(0), C(' '), S("")
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
