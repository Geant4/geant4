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
//
// $Id: G4CameronGilbertPairingCorrections.hh,v 1.8 2002/12/12 19:17:03 gunter Exp $
// GEANT4 tag $Name: geant4-05-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara


#ifndef G4CameronGilbertPairingCorrections_h
#define G4CameronGilbertPairingCorrections_h 1

#include "globals.hh"


//#define verbose 1

class G4CameronGilbertPairingCorrections
{
private:
  // Dummy constructor
  G4CameronGilbertPairingCorrections();
	
  static G4CameronGilbertPairingCorrections* theInstance;


public:
  static G4CameronGilbertPairingCorrections* GetInstance();

  ~G4CameronGilbertPairingCorrections() {};

  G4double GetPairingCorrection(const G4int A, const G4int Z) const
  {
    return GetPairingZ(Z) + GetPairingN(A-Z);
  }

  G4double GetPairingZ(const G4int Z) const 
  {
    if (IsInTableThisZ(Z)) return PairingZTable[Z-ZTableMin]*MeV;
    else {
#ifdef verbose
      G4cerr << "G4CameronGilbertPairingCorrections: out of table for Z = " << Z << G4endl;
#endif
      return 0.0;
    }
  }

  G4bool IsInTableThisZ(const G4int Z) const 
  {
    if ( Z >= ZTableMin && Z <= ZTableMax ) return true;
    else return false;
  }

	
  G4double GetPairingN(const G4int N) const 
  {
   if (IsInTableThisN(N)) return PairingNTable[N-NTableMin]*MeV;
    else {
#ifdef verbose
      G4cerr << "G4CameronGilbertPairingCorrections: out of table for N = " << N << G4endl;
#endif
      return 0.0;
    }
  }

  G4bool IsInTableThisN(const G4int N) const 
  {
    if ( N >= NTableMin && N <= NTableMax ) return true;
    else return false;
  }


  enum  { ZTableSize = 88, NTableSize = 140, ZTableMin = 11, ZTableMax = 98,
	  NTableMin = 11, NTableMax = 150 };

private:



  static const G4double PairingZTable[ZTableSize];

  static const G4double PairingNTable[NTableSize];
	
};
#endif
