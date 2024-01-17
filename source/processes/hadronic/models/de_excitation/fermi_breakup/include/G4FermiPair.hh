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
// FermiBreakUp de-excitation model
// by V. Ivanchenko (July 2016)
//

#ifndef G4FermiPair_h
#define G4FermiPair_h 1

#include "globals.hh"
#include "G4FermiFragment.hh"

class G4FermiPair 
{
public:

  explicit G4FermiPair(const G4FermiFragment* f1, const G4FermiFragment* f2);

  ~G4FermiPair() = default;

  G4int GetA() const { return totalA; }
  G4int GetZ() const { return totalZ; }
  G4double GetMass() const { return mass; } 
  G4double GetExcitationEnergy() const { return excitEnergy; }
  G4double GetTotalEnergy() const { return mass + excitEnergy; }
  const G4FermiFragment* GetFragment1() const { return fragment1; }
  const G4FermiFragment* GetFragment2() const { return fragment2; }

  G4double GetMinMass(G4double Eex) const;

  void SetProbability(const G4double p) { prob = p; }
  G4double Probability() const { return prob; }
  
  inline G4FermiPair(const G4FermiPair &) = delete;
  inline const G4FermiPair & operator=(const G4FermiPair &) = delete;
  inline G4bool operator==(const G4FermiPair &) const = delete;
  inline G4bool operator!=(const G4FermiPair &) const = delete;
  
private:

  G4int totalZ;
  G4int totalA;

  G4double mass;
  G4double excitEnergy;
  G4double prob{1.0};

  const G4FermiFragment* fragment1;
  const G4FermiFragment* fragment2;
};

#endif


