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
// Hadronic Process: Nuclear De-excitations
// by V. Lara
//
// Modified:
// 21.03.2013 V.Ivanchenko redesigned and cleaned up

#ifndef G4CameronShellPlusPairingCorrections_h
#define G4CameronShellPlusPairingCorrections_h 1

#include "globals.hh"

class G4CameronShellPlusPairingCorrections
{
public:

  explicit G4CameronShellPlusPairingCorrections();

  ~G4CameronShellPlusPairingCorrections() = default;

  inline G4bool GetPairingCorrection(G4int N, G4int Z, G4double& result) const
  {
    G4bool res = false;
    if(Z <= TableSize && N <= TableSize) { 
      result = SPZTable[Z-1] + SPNTable[N-1];
      res = true; 
    }
    return res;
  }
  
  enum  { TableSize = 200 };
  
  G4CameronShellPlusPairingCorrections(const G4CameronShellPlusPairingCorrections & right) = delete;
  const G4CameronShellPlusPairingCorrections & operator=
  (const G4CameronShellPlusPairingCorrections & right) = delete;

private:

  static G4double SPZTable[TableSize];
  static G4double SPNTable[TableSize];
  
};
#endif
