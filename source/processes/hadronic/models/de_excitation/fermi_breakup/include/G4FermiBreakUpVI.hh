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

#ifndef G4FermiBreakUpVI_h
#define G4FermiBreakUpVI_h 1

#include "G4VFermiBreakUp.hh"
#include "globals.hh"
#include "G4FermiFragment.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include <vector>

class G4FermiFragmentsPoolVI;

class G4FermiBreakUpVI final: public G4VFermiBreakUp 
{
public:

  G4FermiBreakUpVI();
  ~G4FermiBreakUpVI() override;

  void Initialise() override;

  // check if the Fermi Break Up model can be used 
  // mass is an effective mass of a fragment
  G4bool IsApplicable(G4int ZZ, G4int AA, G4double eexc) const override;

  // new interface - vector of products is added to the provided vector
  // primary fragment is deleted or is modified and added to the list
  // of products 
  void BreakFragment(G4FragmentVector*, G4Fragment* theNucleus) override;

  G4FermiBreakUpVI(const G4FermiBreakUpVI &right) = delete;  
  const G4FermiBreakUpVI & operator=(const G4FermiBreakUpVI &right) = delete;
  G4bool operator==(const G4FermiBreakUpVI &right) const = delete;
  G4bool operator!=(const G4FermiBreakUpVI &right) const = delete;
  
private:

  G4bool SampleDecay(const G4int z, const G4int a, const G4double mass,
                     const G4double ext, G4LorentzVector&);

  static G4FermiFragmentsPoolVI* fPool;

  const G4int maxZ{9};
  const G4int maxA{17};
  G4int secID;  // Creator model ID for the secondaries created by this model
  
  G4double fTolerance{0.0};
  G4double fElim{0.0};
  G4double fTimeLim{1.0}; // in ns

  G4bool isFirst{false};

  std::vector<G4double> prob;
  std::vector<const G4FermiFragment*> frag;
  std::vector<G4LorentzVector> lvect;

};

#endif
