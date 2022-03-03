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
#include <map>
#include <memory>
#include <sstream>
#include "G4MolecularConfiguration.hh"
#include "G4Track.hh"

class G4Voxel
{
 public:
  struct Index
  {
    Index()
      : x(0)
      , y(0)
      , z(0)
    {}
    Index(int x0, int y0, int z0)
      : x(x0)
      , y(y0)
      , z(z0)
    {}

    friend std::ostream& operator<<(std::ostream& stream, const Index& rhs)
    {
      stream << "(" << rhs.x << ", " << rhs.y << ", " << rhs.z << ")";
      return stream;
    }
    G4bool operator==(const Index& rhs) const
    {
      return x == rhs.x && y == rhs.y && z == rhs.z;
    }
    G4bool operator!=(const Index& rhs) const
    {
      return x != rhs.x || y != rhs.y || z != rhs.z;
    }
    Index operator+(const Index& rhs) const
    {
      return { x + rhs.x, y + rhs.y, z + rhs.z };
    }
    int x, y, z;
  };
  using MolType = const G4MolecularConfiguration*;
  using MapList = std::map<MolType, size_t>;

  G4Voxel(MapList&& list, Index& index, G4DNABoundingBox&& box)
    : fMapList(std::move(list))
    , fIndex(index)
    , fBox(std::move(box))
  {}

  ~G4Voxel() = default;

  const Index& GetIndex() const { return fIndex; }

  MapList& GetMapList() { return fMapList; }

  void SetMapList(MapList&& mapList) { fMapList = std::move(mapList); }

  G4double GetVolume() const
  {
    auto xlo = fBox.Getxlo();
    auto ylo = fBox.Getylo();
    auto zlo = fBox.Getzlo();

    auto xhi = fBox.Getxhi();
    auto yhi = fBox.Getyhi();
    auto zhi = fBox.Getzhi();

    return (xhi - xlo) * (yhi - ylo) * (zhi - zlo);
  }

 private:
  MapList fMapList;
  Index fIndex;
  G4DNABoundingBox fBox;
};

class G4DNAMesh
{
 public:
  using Index    = G4Voxel::Index;
  using Key      = unsigned int;
  using VoxelMap = std::map<Key, G4Voxel*>;
  G4DNAMesh(const G4DNABoundingBox&, G4int);
  ~G4DNAMesh();

  Key GetKey(const G4ThreeVector& pos) const;
  Index GetIndex(Key key) const;
  Index GetIndex(const G4ThreeVector& position) const;

  G4Voxel* GetVoxel(Key key);
  [[maybe_unused]] G4Voxel* GetVoxel(const Index& index);
  size_t size() { return fMesh.size(); }
  Index GetIndex(const Index& index, int) const;
  std::vector<Index> FindVoxelNeighbors(const Index& index) const;
  std::vector<Index> FindNeighboringVoxels(const Index& index) const;
  void Reset();

  G4Voxel::MapList& GetVoxelMapList(Key key);
  G4Voxel::MapList& GetVoxelMapList(const Index& index);

  Key GetKey(const Index& index) const;
  VoxelMap::iterator end() { return fMesh.end(); }
  VoxelMap::iterator begin() { return fMesh.begin(); }
  VoxelMap::const_iterator end() const { return fMesh.end(); }
  VoxelMap::const_iterator begin() const { return fMesh.begin(); }

  void PrintMesh();
  void PrintVoxel(const Index& index);
  const G4DNABoundingBox& GetBoundingBox() const;
  G4DNABoundingBox GetBoundingBox(const Index& index);
  G4int GetNumberOfType(G4Voxel::MolType type) const;
  void SetVoxelMapList(const Key& key, G4Voxel::MapList&& mapList);

  G4double GetResolution() const;

 private:
  VoxelMap fMesh;
  const G4DNABoundingBox* fpBoundingMesh;
  G4double fResolution;
};
#endif
