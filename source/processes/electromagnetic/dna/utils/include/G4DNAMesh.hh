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
#ifndef G4DNAMesh_hh
#define G4DNAMesh_hh 1

#include "globals.hh"
#include "G4DNABoundingBox.hh"
#include <vector>
#include <array>
#include <cmath>
#include <unordered_map>
#include <memory>
#include <sstream>
#include "G4MolecularConfiguration.hh"
#include "G4VDNAMesh.hh"

class G4DNAMesh : public G4VDNAMesh
{
 public:
  using Box         = G4DNABoundingBox;
  using MolType     = const G4MolecularConfiguration*;
  using Data        = std::map<MolType, size_t>;
  using Voxel       = std::tuple<Index, Box, Data>;
  using IndexMap    = std::unordered_map<Index, G4int, G4VDNAMesh::hashFunc>;
  using VoxelVector = std::vector<Voxel>;
  G4DNAMesh(const G4DNABoundingBox&, G4int);
  ~G4DNAMesh() override;
  Index GetIndex(const G4ThreeVector& position) const;
  Voxel& GetVoxel(const Index& index);  // GetorCreateVoxel
  size_t size() { return fVoxelVector.size(); };
  Index ConvertIndex(const Index& index, const G4int&) const;
  std::vector<Index> FindNeighboringVoxels(const Index& index) const;
  void Reset();
  Data& GetVoxelMapList(const Index& index);
  auto end() { return fVoxelVector.end(); }
  auto begin() { return fVoxelVector.begin(); }
  VoxelVector::const_iterator const_end() const { return fVoxelVector.end(); }
  VoxelVector::const_iterator const_begin() const
  {
    return fVoxelVector.begin();
  }
  void PrintMesh();
  void PrintVoxel(const Index& index);
  const G4DNABoundingBox& GetBoundingBox() const;
  G4DNABoundingBox GetBoundingBox(const Index& index);
  G4int GetNumberOfType(MolType type) const;
  void InitializeVoxel(const Index& key, Data&& mapList);
  G4double GetResolution() const;

 private:
  IndexMap fIndexMap;
  VoxelVector fVoxelVector;
  const G4DNABoundingBox* fpBoundingMesh;
  G4double fResolution;
};
#endif
