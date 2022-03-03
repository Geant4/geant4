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
#ifndef G4DNAGillespieDirectMethod_hh
#define G4DNAGillespieDirectMethod_hh 1
#include "globals.hh"
#include "G4DNAMesh.hh"
#include "G4Track.hh"
#include "G4Molecule.hh"
#include "G4MolecularConfiguration.hh"
#include "G4MoleculeTable.hh"
#include "G4ITTrackHolder.hh"
#include "G4DNAEventSet.hh"

class G4DNAMolecularReactionTable;
class G4DNAMolecularReactionData;
class G4DNAScavengerMaterial;
class G4DNAGillespieDirectMethod
{
 public:
  G4DNAGillespieDirectMethod();
  ~G4DNAGillespieDirectMethod();
  using MolType      = const G4MolecularConfiguration*;
  using Key          = unsigned int;
  using Index        = G4Voxel::Index;
  using JumpingData  = std::pair<MolType, Index>;
  using ReactionData = const G4DNAMolecularReactionData;
  using EventIt      = G4DNAEventSet::EventSet::iterator;

  G4double PropensityFunction(const Index& index, ReactionData* data);
  G4double PropensityFunction(const Index& index, MolType moleType);
  inline void SetVoxelMesh(G4DNAMesh& mesh) { fpMesh = &mesh; }
  void SetTimeStep(const G4double& stepTime);
  G4double Reaction(const Index& index);
  G4double DiffusiveJumping(const Index& index);
  G4double ComputeNumberInNode(const Index& index, MolType type);
  G4double VolumeOfNode(const Index& index);
  void Initialize();
  void CreateEvent(unsigned int key);
  void SetEventSet(G4DNAEventSet*);

 private:
  G4DNAMolecularReactionTable* fMolecularReactions;
  G4DNAMesh* fpMesh;
  G4double fTimeStep;
  G4DNAEventSet* fpEventSet;
  G4double fVerbose;
  std::map<G4double /*Propensity*/, ReactionData*> fReactionDataMap;
  std::map<G4double /*Propensity*/, JumpingData> fJumpingDataMap;
  G4bool FindScavenging(const Index& index, MolType, G4double&);
  G4DNAScavengerMaterial* fpScavengerMaterial;
};
#endif
