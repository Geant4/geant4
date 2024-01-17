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
// G4MuonicAtomHelper

// Author: K.Lynch, 01.07.2016 - First implementation
// Revision:
// - 12.06.2017, K L Genser
// --------------------------------------------------------------------
#ifndef G4MuonicAtomHelper_hh
#define G4MuonicAtomHelper_hh 1

#include "G4Ions.hh"
#include "G4MuonicAtom.hh"

class G4MuonicAtomHelper
{
  public:
    static G4MuonicAtom* ConstructMuonicAtom(const G4String& name, G4int encoding,
                                             G4Ions const* baseion);

    static G4double GetMuonCaptureRate(G4int Z, G4int A);

    static G4double GetMuonDecayRate(G4int Z);

    static G4double GetMuonZeff(G4int Z);

    static G4double GetKShellEnergy(G4double A);

    static G4double GetLinApprox(G4int N, const G4double* const X, const G4double* const Y,
                                 G4double Xuser);
};

#endif
