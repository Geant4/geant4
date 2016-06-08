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
#ifndef G4BaryonSplitter_h
#define G4BaryonSplitter_h

// HPW Feb1999 based on prototype, needs urgent clean-up of data structures.
// Also needs clean-up of interfaces.
// @@@@@@@@@@@@@@@@@
// clean-up of data structures

#include "globals.hh"
#include "G4SPBaryonTable.hh"

class G4BaryonSplitter
{
public:
  G4BaryonSplitter();
  G4bool SplitBarion(G4int PDGCode, G4int* q_or_qqbar, G4int* qbar_or_qq);
  G4bool FindDiquark(G4int PDGCode, G4int Quark, G4int* Diquark);
  const G4SPBaryon & GetSPBaryon(G4int PDGCode);

private:

  G4SPBaryonTable theBaryons;
};

#endif
