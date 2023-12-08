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
// 070625 add natural abundance (nat) flag by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
// V. Ivanchenko, July-2023 Basic revision of particle HP classes
//

#ifndef G4ParticleHPDataUsed_h
#define G4ParticleHPDataUsed_h 1

#include "globals.hh"

class G4ParticleHPDataUsed
{
public:
  G4ParticleHPDataUsed() = default;

  void SetA(G4double anA) { theA = G4lrint(anA); }
  void SetA(G4int anA) { theA = anA; }
  void SetZ(G4int aZ) { theZ = aZ; }
  void SetM(G4int aM) { theM = aM; }
  void SetName(const G4String& aName) { theName = aName; }

  G4int GetZ() { return theZ; }
  G4int GetA() { return theA; }
  G4int GetM() { return theM; }
  G4String& GetName() { return theName; }

  G4bool IsThisNaturalAbundance() { return nat; };
  void SetNaturalAbundanceFlag() { nat = TRUE; };

private:
  G4int theA{0};
  G4int theZ{0};
  G4int theM{0};
  G4bool nat{true};
  G4String theName{""};
};

#endif
