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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4CameronShellPlusPairingCorrections.hh,v 1.7 2002/01/15 12:03:00 vlara Exp $
// GEANT4 tag $Name: geant4-04-00-patch-02 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4CameronShellPlusPairingCorrections_h
#define G4CameronShellPlusPairingCorrections_h 1

#include "globals.hh"

//#define verbose 1

class G4CameronShellPlusPairingCorrections
{
private:

  // Dummy constructor
  G4CameronShellPlusPairingCorrections();

  static G4CameronShellPlusPairingCorrections* theInstance;
	
public:
	
  static G4CameronShellPlusPairingCorrections* GetInstance();

  ~G4CameronShellPlusPairingCorrections() {};
  
  G4double GetShellPlusPairingZ(const G4int Z) const 
  {
    if (Z <= TableSize && Z > 1) return SPZTable[Z-1]*MeV;
    else {
#ifdef verbose
      G4cerr << "G4CameronShellPlusPairingCorrections: out of table for Z = " << Z << G4endl;
#endif
      return 0.0;
    }
  }
  
  G4double GetShellPlusPairingN(const G4int N) const 
  {
    if (N <= TableSize && N > 0) return SPNTable[N-1]*MeV;
    else {
#ifdef verbose
      G4cerr << "G4CameronShellPlusPairingCorrections: out of table for N = " << N << G4endl;
#endif
      return 0.0;
    }
  }
  
  
  enum  { TableSize = 200 };
  
private:
  
  const static G4double SPZTable[TableSize];
  
  const static G4double SPNTable[TableSize];
  
};
#endif
