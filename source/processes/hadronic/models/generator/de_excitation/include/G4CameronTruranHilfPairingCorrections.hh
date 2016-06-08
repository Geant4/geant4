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
// $Id: G4CameronTruranHilfPairingCorrections.hh,v 1.7 2002/01/15 12:03:12 vlara Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4CameronTruranHilfPairingCorrections_h
#define G4CameronTruranHilfPairingCorrections_h 1

#include "globals.hh"

//#define verbose 1

class G4CameronTruranHilfPairingCorrections
{
private:

  G4CameronTruranHilfPairingCorrections();
	
  static G4CameronTruranHilfPairingCorrections* theInstance;


public:
  static G4CameronTruranHilfPairingCorrections* GetInstance();

  ~G4CameronTruranHilfPairingCorrections() {};

  G4double GetParingCorrection(const G4int A, const  G4int Z) const
  {
    return GetPairingZ(Z) + GetPairingN(A-Z);
  }


  G4double GetPairingZ(const G4int Z) const 
  {
    if (IsInTableThisZ(Z)) return -PairingZTable[Z-ZTableMin]*MeV; // Notice the sign 
    else {
#ifdef verbose
      G4cerr << "G4CameronTruranHilfPairingCorrections: out of table for Z = " << Z << G4endl;
#endif
      return 0.0;
    }
  }
  
  G4bool IsInTableThisZ(const G4int Z) const 
  {
    if ( Z >= ZTableMin && Z <= ZTableMax ) return true;
    else return false;
  }

	
  G4double GetPairingN(const G4int N)  const 
  {
    if (IsInTableThisN(N)) return -PairingNTable[N-NTableMin]*MeV; // Notice the sign
    else {
#ifdef verbose
      G4cerr << "G4CameronTruranHilfPairingCorrections: out of table for N = " << N << G4endl;
#endif
      return 0.0;
    }
  }
  
  G4bool IsInTableThisN(const G4int N)  const 
  {
    if ( N >= NTableMin && N <= NTableMax ) return true;
    else return false;
  }


  enum  { ZTableSize = 93, NTableSize = 146, ZTableMin = 10, ZTableMax = 102,
	  NTableMin = 10, NTableMax = 155 };
private:



    static const G4double PairingZTable[ZTableSize];

    static const G4double PairingNTable[NTableSize];
	
};
#endif
