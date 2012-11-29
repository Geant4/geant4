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
// $Id$
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara


#ifndef G4CookShellCorrections_h
#define G4CookShellCorrections_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh" 

//#define verbose 1

class G4CookShellCorrections
{
private:

  G4CookShellCorrections();
	
  static G4CookShellCorrections* theInstance;


public:
	
  static G4CookShellCorrections* GetInstance();

  ~G4CookShellCorrections();

  G4double GetShellCorrection(G4int A, G4int Z) const 
  {
    return GetShellZ(Z) + GetShellN(A-Z); 
  }

  G4double GetShellZ(const G4int Z) const 
  {
    if ( this->IsInTableThisZ(Z) ) return ShellZTable[Z-ZTableMin]*CLHEP::MeV;
    else {
#ifdef verbose
      G4cerr << "G4CookShellCorrections: out of table for Z = " << Z << G4endl;
#endif
      return 0.0;
    }
  }

  G4bool IsInTableThisZ(const G4int Z) const 
  {
    if ( Z >= ZTableMin && Z <= ZTableMax ) return true;
    else return false;
  }
  
  G4double GetShellN(const G4int N) const 
  {
    if ( this->IsInTableThisN(N) ) return ShellNTable[N-NTableMin]*CLHEP::MeV;
    else {
#ifdef verbose
      G4cerr << "G4CookShellCorrections: out of table for N = " << N << G4endl;
#endif
      return 0.0;
    }
  }
    
  G4bool IsInTableThisN(const G4int N) const 
  {
    if ( N >= NTableMin && N <= NTableMax ) return true;
    else return false;
  }
  
  enum  { ZTableSize = 68, NTableSize = 118, ZTableMin = 28, ZTableMax = 95,
	  NTableMin = 33, NTableMax = 150 };
private:
  
  
  
  static const G4double ShellZTable[ZTableSize];
  
  static const G4double ShellNTable[NTableSize];
  
};
#endif
