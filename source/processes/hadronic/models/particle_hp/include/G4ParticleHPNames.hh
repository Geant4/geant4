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
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#ifndef G4ParticleHPNames_h
#define G4ParticleHPNames_h 1

#include "G4ios.hh"
#include <fstream>
// #include <strstream>
#include <stdlib.h>
#include "globals.hh"
#include "G4ParticleHPDataUsed.hh"

class G4ParticleHPNames
{
  public:
  
  G4ParticleHPNames(){theMaxOffSet = 5;}
  G4ParticleHPNames(G4int maxOffSet){theMaxOffSet = maxOffSet;}
  ~G4ParticleHPNames(){}
  
  //G4ParticleHPDataUsed GetName(G4int A, G4int Z, G4String base, G4String rest, G4bool & active);
  G4ParticleHPDataUsed GetName(G4int A, G4int Z, G4String base, G4String rest, G4bool & active) { G4int M = 0; return GetName( A, Z, M, base, rest, active); };
  G4ParticleHPDataUsed GetName(G4int A, G4int Z, G4int M, G4String base, G4String rest, G4bool & active);
  G4String GetName(G4int i);
  void SetMaxOffSet(G4int anOffset) { theMaxOffSet = anOffset; }
  
  public:
  
  static const G4String theString[100];
  G4int theMaxOffSet;
  G4String itoa(int current)
  {
    const char theDigits[11] = "0123456789";
    G4String result;
    int digit;
    do
    {
      digit = current-10*(current/10);
      result=theDigits[digit]+result;
      current/=10;
    }
    while(current!=0); // Loop checking, 11.05.2015, T. Koi
    return result;
  }
};

#endif
