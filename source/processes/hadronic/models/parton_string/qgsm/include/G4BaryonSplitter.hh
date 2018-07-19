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

