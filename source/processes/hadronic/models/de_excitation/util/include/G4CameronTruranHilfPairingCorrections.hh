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
// $Id: G4CameronTruranHilfPairingCorrections.hh,v 1.5 2009-03-04 11:05:02 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

  ~G4CameronTruranHilfPairingCorrections();

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
