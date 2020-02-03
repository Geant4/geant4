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
/// @file G4MPIrandomSeedGenerator.cc
/// @brief An implementation of random number seed distribution

#include "mpi.h"
#include "Randomize.hh"
#include "G4MPImanager.hh"
#include "G4MPIrandomSeedGenerator.hh"

// --------------------------------------------------------------------------
G4MPIrandomSeedGenerator::G4MPIrandomSeedGenerator()
  : G4VMPIseedGenerator()
{
  GenerateSeeds();
}

// --------------------------------------------------------------------------
G4MPIrandomSeedGenerator::~G4MPIrandomSeedGenerator()
{
}

// --------------------------------------------------------------------------
G4bool G4MPIrandomSeedGenerator::CheckDoubleCount()
{
  G4MPImanager* g4mpi = G4MPImanager::GetManager();

  G4int nsize = g4mpi-> GetActiveSize();

  for ( G4int i = 0; i < nsize; i++ ) {
    for ( G4int j = 0; j < nsize; j++ ) {
      if( (i != j) && (seed_list_[i] == seed_list_[j]) ) {
        G4double x = G4UniformRand();
        seed_list_[j] = G4long(x*LONG_MAX);
        return false;
      }
    }
  }

  return true;
}

// --------------------------------------------------------------------------
void G4MPIrandomSeedGenerator::GenerateSeeds()
{
  G4MPImanager* g4mpi = G4MPImanager::GetManager();

  G4int nsize = g4mpi-> GetActiveSize();

  seed_list_.clear();

  for ( G4int i = 0; i < nsize; i++ ) {
    G4double x = G4UniformRand();
    G4int seed = G4long(x*LONG_MAX);
    seed_list_.push_back(seed);
  }

  while(! CheckDoubleCount()) {
    continue;
  }
}
