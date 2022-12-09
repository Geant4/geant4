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
class G4DNAMolecularReactionData;
class G4DNAScavengerMaterial;
class G4MolecularConfiguration;
class G4DNAGillespieDirectMethod
{
 public:
  G4DNAGillespieDirectMethod();
  ~G4DNAGillespieDirectMethod();
  using MolType      = const G4MolecularConfiguration*;
  using Index        = G4VDNAMesh::Index;
  using Voxel        = G4DNAMesh::Voxel;
  using JumpingData  = std::pair<MolType, Index>;
  using ReactionData = const G4DNAMolecularReactionData;
  using EventIt      = G4DNAEventSet::EventSet::iterator;

  G4double PropensityFunction(const Voxel& voxel, ReactionData* data);
  G4double PropensityFunction(const Voxel& voxel, MolType moleType);
  inline void SetVoxelMesh(G4DNAMesh& mesh) { fpMesh = &mesh; }
  void SetTimeStep(const G4double& stepTime);
  void Initialize();
  void CreateEvent(const Index& index);
  void SetEventSet(G4DNAEventSet*);

 private:
  G4double Reaction(const Voxel& voxel);
  G4double DiffusiveJumping(const Voxel& voxel);
  G4double ComputeNumberInNode(const Voxel& voxel, MolType type);
  G4double VolumeOfNode(const Voxel& voxel);
  G4DNAMolecularReactionTable* fMolecularReactions;
  G4DNAMesh* fpMesh = nullptr;
  G4double fTimeStep = DBL_MAX;
  G4DNAEventSet* fpEventSet = nullptr;
  G4double fVerbose         = 0;
  std::map<G4double /*Propensity*/, ReactionData*> fReactionDataMap;
  std::map<G4double /*Propensity*/, JumpingData> fJumpingDataMap;
  G4bool FindScavenging(const Voxel& voxel, MolType, G4double&);
  G4DNAScavengerMaterial* fpScavengerMaterial = nullptr;
};
#endif
