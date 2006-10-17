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
//    **************************************
//    *                                    *
//    *    CellSurvivalRT112LQLambin.hh    *
//    *                                    *
//    **************************************
//
// Author: Barbara Mascialino (Barbara.Mascialino@ge.infn.it)
//
//
// History:
// -----------
// 12 October 2006 B. Mascialino      first implementation
// -------------------------------------------------------------------


#ifndef CellSurvivalRT112LQLambin_h
#define CellSurvivalRT112LQLambin_h 1

#include "globals.hh"

class CellSurvivalRT112LQLambin 
{
  public:
    CellSurvivalRT112LQLambin();
   ~CellSurvivalRT112LQLambin();

  public:
  void  SurvivalFormula(G4double dose);
  G4double GetSurvival();

private:
  G4double probability;
};
#endif
