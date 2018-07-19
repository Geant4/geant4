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
// $Id: G4FermiBreakUpVI.hh 96711 2016-05-02 10:01:24Z vnivanch $
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
#include "G4Threading.hh"
#include <CLHEP/Random/RandomEngine.h>
#include <vector>

class G4FermiFragmentsPoolVI;
class G4FermiDecayProbability;

class G4FermiBreakUpVI : public G4VFermiBreakUp 
{
public:

  explicit G4FermiBreakUpVI();
  virtual ~G4FermiBreakUpVI();

  virtual void Initialise() final;

  // check if the Fermi Break Up model can be used 
  // mass is an effective mass of a fragment
  virtual G4bool IsApplicable(G4int ZZ, G4int AA, G4double etot) const final;

  // new interface - vector of products is added to the provided vector
  // primary fragment is deleted or is modified and added to the list
  // of products 
  virtual void BreakFragment(G4FragmentVector*, G4Fragment* theNucleus) final;
  
private:

  virtual void InitialisePool() final;

  G4bool SampleDecay();

  G4FermiBreakUpVI(const G4FermiBreakUpVI &right) = delete;  
  const G4FermiBreakUpVI & operator=(const G4FermiBreakUpVI &right) = delete;
  G4bool operator==(const G4FermiBreakUpVI &right) const = delete;
  G4bool operator!=(const G4FermiBreakUpVI &right) const = delete;

  static G4FermiFragmentsPoolVI* thePool;
  const  G4FermiDecayProbability* theDecay;

  CLHEP::HepRandomEngine* rndmEngine;

  G4int verbose;
  G4int maxZ;
  G4int maxA;

  G4int Z;
  G4int A;
  G4int spin;

  G4double mass;
  G4double excitation;
  G4double tolerance;
  G4double elim;

  const G4FermiFragment* frag1;
  const G4FermiFragment* frag2;

  G4LorentzVector lv0;
  G4ThreeVector boostVector;

  std::vector<G4double> prob;
  std::vector<const G4FermiFragment*> frag;
  std::vector<G4LorentzVector> lvect;

#ifdef G4MULTITHREADED
  static G4Mutex FermiBreakUpVIMutex;
#endif
};

#endif
