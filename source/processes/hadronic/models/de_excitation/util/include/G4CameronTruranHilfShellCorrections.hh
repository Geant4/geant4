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
// $Id: G4CameronTruranHilfShellCorrections.hh,v 1.1 2003/08/26 18:50:11 lara Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4CameronTruranHilfShellCorrections_h
#define G4CameronTruranHilfShellCorrections_h 1

#include "globals.hh"

//#define verbose 1

class G4CameronTruranHilfShellCorrections
{
private:
	
  G4CameronTruranHilfShellCorrections();
	
  static G4CameronTruranHilfShellCorrections* theInstance;


public:
  static G4CameronTruranHilfShellCorrections* GetInstance();
  
  ~G4CameronTruranHilfShellCorrections() {};

  G4double GetShellCorrection(const G4int A, const G4int Z) const 
  {
    return GetShellZ(Z) + GetShellN(A-Z);
  }

  G4double GetShellZ(const G4int Z) const 
  {
    if (IsInTableThisZ(Z)) return ShellZTable[Z-ZTableMin]*MeV;
    else {
#ifdef verbose
      G4cerr << "G4CameronTruranHilfShellCorrections: out of table for Z = " << Z << G4endl;
#endif
      return 0.0;
    }
  }
  
  G4bool IsInTableThisZ(const G4int Z)  const 
  {
    if ( Z >= ZTableMin && Z <= ZTableMax ) return true;
    else return false;
  }
  
  G4double GetShellN(const G4int N) const 
  {
    if (IsInTableThisN(N)) return ShellNTable[N-NTableMin]*MeV;
    else {
#ifdef verbose
      G4cerr << "G4CameronTruranHilfShellCorrections: out of table for N = " << N << G4endl;
#endif
      return 0.0;
    }
  }
  
  G4bool IsInTableThisN(const G4int N) const 
  {
    if ( N >= NTableMin && N <= NTableMax ) return true;
    else return false;
  }
  
  enum  { ZTableSize = 93, NTableSize = 146, ZTableMin = 10, ZTableMax = 102,
	  NTableMin = 10, NTableMax = 155 };

private:



  static const G4double ShellZTable[ZTableSize];
  
  static const G4double ShellNTable[NTableSize];
	
};
#endif
