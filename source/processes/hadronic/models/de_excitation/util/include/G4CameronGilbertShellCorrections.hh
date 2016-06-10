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
// $Id: G4CameronGilbertShellCorrections.hh 68724 2013-04-05 09:26:32Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara
//
// Modified:
// 21.03.2013 V.Ivanchenko redesigned and cleaned up

#ifndef G4CameronGilbertShellCorrections_h
#define G4CameronGilbertShellCorrections_h 1

#include "globals.hh"

class G4CameronGilbertShellCorrections
{
public:

  G4CameronGilbertShellCorrections();

  ~G4CameronGilbertShellCorrections();
  
  inline 
  G4double GetShellCorrection(G4int A, G4int Z) const 
  {
    return GetShellZ(Z) + GetShellN(A-Z);
  }

  inline 
  G4double GetShellZ(G4int Z) const
  {
    G4double res = 0.0;
    if (IsInTableThisZ(Z)) { res = ShellZTable[Z-ZTableMin]; }
    return res;
  }

  inline 
  G4bool IsInTableThisZ(G4int Z) const 
  {
    return ( Z >= ZTableMin && Z <= ZTableMax );
  }
	
  inline 
  G4double GetShellN(G4int N) const 
  {
    G4double res = 0.0;
    if (IsInTableThisN(N)) { res = ShellNTable[N-NTableMin]; }
    return res;
  }
  
  inline 
  G4bool IsInTableThisN(G4int N) const 
  {
    return ( N >= NTableMin && N <= NTableMax );
  }

  enum  { ZTableSize = 88, NTableSize = 140, ZTableMin = 11, ZTableMax = 98,
	  NTableMin = 11, NTableMax = 150 };

private:

  static G4double ShellZTable[ZTableSize];  
  static G4double ShellNTable[NTableSize];
	
};
#endif
