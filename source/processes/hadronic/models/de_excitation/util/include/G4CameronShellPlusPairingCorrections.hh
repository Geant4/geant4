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
// $Id: G4CameronShellPlusPairingCorrections.hh 68724 2013-04-05 09:26:32Z gcosmo $
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

  G4CameronShellPlusPairingCorrections();

  ~G4CameronShellPlusPairingCorrections();
  
  inline
  G4double GetShellPlusPairingZ(G4int Z) const 
  {
    G4double res = 0.0;
    if (Z <= TableSize && Z > 1) { res = SPZTable[Z-1]; }
    return res;
  }
  
  inline
  G4double GetShellPlusPairingN(G4int N) const 
  {
    G4double res = 0.0;
    if (N <= TableSize && N > 0) { res = SPNTable[N-1]; }
    return res;
  }
  
  enum  { TableSize = 200 };
  
private:
  
  static G4double SPZTable[TableSize];
  static G4double SPNTable[TableSize];
  
};
#endif
