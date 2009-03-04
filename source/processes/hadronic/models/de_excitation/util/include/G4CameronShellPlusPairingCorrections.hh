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
// $Id: G4CameronShellPlusPairingCorrections.hh,v 1.5 2009-03-04 11:05:02 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

  ~G4CameronShellPlusPairingCorrections();
  
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
  
  static const G4double SPZTable[TableSize];
  
  static const G4double SPNTable[TableSize];
  
};
#endif
