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
/// @file G4VMPIseedGenerator.hh
/// @brief A base class for random number seed distribution

#ifndef G4VMPI_SEED_GENERATOR_H
#define G4VMPI_SEED_GENERATOR_H

#include "globals.hh"
#include <vector>

class G4VMPIseedGenerator {
public:
  G4VMPIseedGenerator();
  virtual ~G4VMPIseedGenerator();

  // set/get methods
  void SetMasterSeed(G4long aseed); // trigger GenerateSeeds()
  G4long GetMasterSeed() const;

  const std::vector<G4long>& GetSeedList() const;

protected:
  G4long master_seed_;
  std::vector<G4long> seed_list_;

  // generate seeds for MPI nodes
  virtual void GenerateSeeds() = 0;
};

// ====================================================================
inline G4long G4VMPIseedGenerator::GetMasterSeed() const
{
  return master_seed_;
}

inline const std::vector<G4long>& G4VMPIseedGenerator::GetSeedList() const
{
  return seed_list_;
}

#endif
