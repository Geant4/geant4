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
// Hadronic Process: Nuclear De-excitations
// by V. Lara
//
// Modification: 13.08.2025 V.Ivanchenko rewrite

#ifndef G4StatMFMicroCanonical_h
#define G4StatMFMicroCanonical_h 1

#include <vector>

#include "globals.hh"
#include "G4VStatMFEnsemble.hh"
#include "G4StatMFMicroPartition.hh"
#include "G4StatMFMicroManager.hh"
#include "G4StatMFParameters.hh"
#include "G4StatMFChannel.hh"

#include "G4Fragment.hh"
#include "G4VStatMFMacroCluster.hh"
#include "G4FunctionSolver.hh"

class G4Pow;

class G4StatMFMicroCanonical : public G4VStatMFEnsemble {

public:

  G4StatMFMicroCanonical();

  ~G4StatMFMicroCanonical() override;

  // Initialise for a given G4Fragment
  void Initialise(const G4Fragment& theFragment) override;

  // Choice of the channel
  G4StatMFChannel* ChooseAandZ(const G4Fragment &theFragment) override;

  G4double Function(G4double T)
  { return (fExEnergy + pFreeInternalE0 - CalcFreeInternalEnergy(T)); }

  // copy constructor
  G4StatMFMicroCanonical(const G4StatMFMicroCanonical& right) = delete;
  G4StatMFMicroCanonical& operator=(const G4StatMFMicroCanonical& right) = delete;
  G4bool operator==(const G4StatMFMicroCanonical& right) const = delete;
  G4bool operator!=(const G4StatMFMicroCanonical& right) const = delete;

private:

  // Calculate Entropy of Compound Nucleus
  G4double CalcEntropyOfCompoundNucleus(G4double& T);

  G4double CalcFreeInternalEnergy(G4double T);

  G4int Z{0};
  G4int A{0};
	
  // Statistical weight of compound nucleus
  G4double fWCompoundNucleus{0.0};
  G4double fExEnergy{0.0};

  G4double A13{0.0};
  G4double fInvLevelDensity{1.0};
  G4double fSymmetryTerm{1.0};
  G4double fCoulombTerm{0.0};
  
  G4Pow* g4calc;
  G4FunctionSolver<G4StatMFMicroCanonical>* fSolver;
  // This is a vector of partitions provided different multiplicities
  std::vector<G4StatMFMicroManager*> fPartitionManagerVector;

};

#endif
