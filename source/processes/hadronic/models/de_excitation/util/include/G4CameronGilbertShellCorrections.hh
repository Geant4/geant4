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

#ifndef G4CameronGilbertShellCorrections_h
#define G4CameronGilbertShellCorrections_h 1

#include "globals.hh"

class G4CameronGilbertShellCorrections
{
public:

  G4CameronGilbertShellCorrections();

  ~G4CameronGilbertShellCorrections() = default;

  inline G4bool GetShellCorrection(G4int N, G4int Z, G4double& result) const
  {
    G4bool res = false;
    if(Z >= ZTableMin && Z <= ZTableMax && N >= NTableMin && N <= NTableMax) { 
      result = ShellZTable[Z-ZTableMin] + ShellNTable[N-NTableMin];
      res = true; 
    }
    return res;
  }
  
  enum  { ZTableSize = 88, NTableSize = 140, ZTableMin = 11, ZTableMax = 98,
	  NTableMin = 11, NTableMax = 150 };

  G4CameronGilbertShellCorrections(const G4CameronGilbertShellCorrections & right) = delete;
  const G4CameronGilbertShellCorrections & operator=
  (const G4CameronGilbertShellCorrections & right) = delete;

private:

  static G4double ShellZTable[ZTableSize];  
  static G4double ShellNTable[NTableSize];
	
};
#endif
