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
#include "G4DNAMesh.hh"
#include <algorithm>
#include <cassert>
#include <ostream>

G4DNAMesh::G4DNAMesh(const G4DNABoundingBox& boundingBox, G4int pixel)
  : fpBoundingMesh(&boundingBox)
  , fResolution((2 * boundingBox.halfSideLengthInY() / pixel))
{}

G4DNAMesh::~G4DNAMesh() { Reset(); }

G4Voxel::MapList& G4DNAMesh::GetVoxelMapList(Key key)
{
  auto iter = fMesh.find(key);
  if(iter == fMesh.end())
  {
    G4Voxel::MapList maplist;
    SetVoxelMapList(key, std::move(maplist));
    return GetVoxelMapList(key);
  }
  else
  {
    return iter->second->GetMapList();
  }
}

G4Voxel::MapList& G4DNAMesh::GetVoxelMapList(const Index& index)
{
  auto key = GetKey(index);
  return GetVoxelMapList(key);
}

G4DNAMesh::Key G4DNAMesh::GetKey(const Index& index) const
{
  auto xmax = (unsigned int) (std::floor(
    (fpBoundingMesh->Getxhi() - fpBoundingMesh->Getxlo()) / fResolution));
  auto ymax = (unsigned int) (std::floor(
    (fpBoundingMesh->Getyhi() - fpBoundingMesh->Getylo()) / fResolution));
  return index.z * ymax * xmax + index.y * xmax + index.x;
}

void G4DNAMesh::PrintMesh()
{
  G4cout << "*********PrintMesh::Size : " << fMesh.size() << G4endl;
  auto iter = fMesh.begin();
  for(; iter != fMesh.end(); iter++)
  {
    auto index = iter->second->GetIndex();
    PrintVoxel(index);
  }
  G4cout << G4endl;
}
G4int G4DNAMesh::GetNumberOfType(G4Voxel::MolType type) const
{
  G4int output = 0;
  auto iter    = fMesh.begin();
  for(; iter != fMesh.end(); iter++)
  {
    auto node = dynamic_cast<G4Voxel*>(iter->second);
    if(node == nullptr)
    {
      continue;
    }
    auto it = node->GetMapList().find(type);
    if(it != node->GetMapList().end())
    {
      output += it->second;
    }
  }
  return output;
}

void G4DNAMesh::PrintVoxel(const Index& index)
{
  G4cout << "*********PrintVoxel::";
  G4cout << "key: " << GetKey(index) << " index : " << index
         << " number of type : " << this->GetVoxelMapList(index).size()
         << G4endl;

  for(const auto& it : this->GetVoxelMapList(index))
  {
    G4cout << "_____________" << it.first->GetName() << " : " << it.second
           << G4endl;
  }
  G4cout << G4endl;
}

void G4DNAMesh::SetVoxelMapList(const Key& key, G4Voxel::MapList&& mapList)
{
  auto index  = GetIndex(key);
  auto pVoxel = fMesh[key];
  if(nullptr == pVoxel)
  {
    pVoxel     = new G4Voxel(std::move(mapList), index, GetBoundingBox(index));
    fMesh[key] = pVoxel;
  }
  else
  {
    assert(pVoxel->GetMapList().empty());  // check if map list is empty
    pVoxel->SetMapList(std::move(mapList));
  }
}

G4Voxel::Index G4DNAMesh::GetIndex(const G4ThreeVector& position) const
{
  int dx = std::floor((position.x() - fpBoundingMesh->Getxlo()) / fResolution);
  int dy = std::floor((position.y() - fpBoundingMesh->Getylo()) / fResolution);
  int dz = std::floor((position.z() - fpBoundingMesh->Getzlo()) / fResolution);
  assert(dx >= 0 && dy >= 0 && dz >= 0);
  return G4Voxel::Index{ dx, dy, dz };
}
G4Voxel::Index G4DNAMesh::GetIndex(const Index& index, int pixels) const
{
  int xmax =
    std::floor((fpBoundingMesh->Getxhi() - fpBoundingMesh->Getxlo()) / fResolution);
  int ymax =
    std::floor((fpBoundingMesh->Getyhi() - fpBoundingMesh->Getylo()) / fResolution);
  int zmax =
    std::floor((fpBoundingMesh->Getzhi() - fpBoundingMesh->Getzlo()) / fResolution);
  int dx = (int) (index.x * pixels / xmax);
  int dy = (int) (index.y * pixels / ymax);
  int dz = (int) (index.z * pixels / zmax);
  assert(dx >= 0 && dy >= 0 && dz >= 0);
  return Index{ dx, dy, dz };
}

G4DNABoundingBox G4DNAMesh::GetBoundingBox(const Index& index)
{
  auto xlo = fpBoundingMesh->Getxlo() + index.x * fResolution;
  auto ylo = fpBoundingMesh->Getylo() + index.y * fResolution;
  auto zlo = fpBoundingMesh->Getzlo() + index.z * fResolution;

  auto xhi = fpBoundingMesh->Getxlo() + (index.x + 1) * fResolution;
  auto yhi = fpBoundingMesh->Getylo() + (index.y + 1) * fResolution;
  auto zhi = fpBoundingMesh->Getzlo() + (index.z + 1) * fResolution;
  return G4DNABoundingBox({ xhi, xlo, yhi, ylo, zhi, zlo });
}

G4Voxel::Index G4DNAMesh::GetIndex(Key key) const
{
  G4int xmax =
    std::floor((fpBoundingMesh->Getxhi() - fpBoundingMesh->Getxlo()) / fResolution);
  G4int ymax =
    std::floor((fpBoundingMesh->Getyhi() - fpBoundingMesh->Getylo()) / fResolution);
  G4int id = key;
  G4int x_ = id % xmax;
  id /= xmax;
  G4int y_ = id % ymax;
  id /= ymax;
  G4int z_ = id;

  if(xmax != ymax)
  {
    G4cout << xmax << " " << ymax << " " << key << G4endl;
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "xmax != ymax";
    G4Exception("G4DNAMesh::GetIndex", "G4DNAMesh006", FatalErrorInArgument,
                exceptionDescription);
  }

  if(x_ < 0 || y_ < 0 || z_ < 0)
  {
    G4cout << xmax << " " << ymax << " " << key << G4endl;
    G4cout << x_ << " " << y_ << " " << z_ << G4endl;
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "x_ < 0 || y_ < 0 || z_ < 0";
    G4Exception("G4DNAMesh::GetIndex", "G4DNAMesh005", FatalErrorInArgument,
                exceptionDescription);
  }
  return Index{ x_, y_, z_ };
}

void G4DNAMesh::Reset()
{
  if(fMesh.empty())
  {
    return;
  }
  for(auto iter : fMesh)  // should use smart ptr
  {
    delete iter.second;
  }
  fMesh.clear();
}

const G4DNABoundingBox& G4DNAMesh::GetBoundingBox() const
{
  return *fpBoundingMesh;
}

G4Voxel* G4DNAMesh::GetVoxel(Key key)
{
  auto it = fMesh.find(key);
  if(it != fMesh.end())
  {
    return it->second;
  }
  return nullptr;
}

[[maybe_unused]] G4Voxel* G4DNAMesh::GetVoxel(const Index& index)
{
  return GetVoxel(GetKey(index));
}

std::vector<G4Voxel::Index>  // array is better ?
G4DNAMesh::FindVoxelNeighbors(const Index& index) const
{
  std::vector<Index> neighbors;

  auto xMax = (int) (std::floor(
    (fpBoundingMesh->Getxhi() - fpBoundingMesh->Getxlo()) / fResolution));
  auto yMax = (int) (std::floor(
    (fpBoundingMesh->Getyhi() - fpBoundingMesh->Getylo()) / fResolution));
  auto zMax = (int) (std::floor(
    (fpBoundingMesh->Getzhi() - fpBoundingMesh->Getzlo()) / fResolution));

  auto xmin = (index.x - 1) < 0 ? 0 : (index.x - 1);
  auto ymin = (index.y - 1) < 0 ? 0 : (index.y - 1);
  auto zmin = (index.z - 1) < 0 ? 0 : (index.z - 1);

  auto xmax = (index.x + 1) > xMax ? xMax : (index.x + 1);
  auto ymax = (index.y + 1) > yMax ? yMax : (index.y + 1);
  auto zmax = (index.z + 1) > zMax ? zMax : (index.z + 1);
  for(int ix = xmin; ix <= xmax; ix++)
  {
    for(int iy = ymin; iy <= ymax; iy++)
    {
      for(int iz = zmin; iz <= zmax; iz++)
      {
        auto key = iz * yMax * xMax + iy * xMax + ix;
        if(GetIndex(key) != index)
        {  // deleting the middle element
          neighbors.push_back(GetIndex(key));
        }
      }
    }
  }
  if(neighbors.empty())
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "neighbors.empty()";
    G4Exception("G4DNAMesh::FindVoxelNeighbors", "G4DNAMesh001",
                FatalErrorInArgument, exceptionDescription);
  }

  return neighbors;
}

std::vector<G4Voxel::Index>  // array is better ?
G4DNAMesh::FindNeighboringVoxels(const Index& index) const
{
  std::vector<Index> neighbors;
  // auto key = GetKey(index);
  auto xMax = (int) (std::floor(
    (fpBoundingMesh->Getxhi() - fpBoundingMesh->Getxlo()) / fResolution));
  auto yMax = (int) (std::floor(
    (fpBoundingMesh->Getyhi() - fpBoundingMesh->Getylo()) / fResolution));
  auto zMax = (int) (std::floor(
    (fpBoundingMesh->Getzhi() - fpBoundingMesh->Getzlo()) / fResolution));

  if(index.x - 1 >= 0)
  {
    neighbors.emplace_back(Index(index.x - 1, index.y, index.z));
  }
  if(index.y - 1 >= 0)
  {
    neighbors.emplace_back(Index(index.x, index.y - 1, index.z));
  }
  if(index.z - 1 >= 0)
  {
    neighbors.emplace_back(Index(index.x, index.y, index.z - 1));
  }
  if(index.x + 1 < xMax)
  {
    neighbors.emplace_back(Index(index.x + 1, index.y, index.z));
  }
  if(index.y + 1 < yMax)
  {
    neighbors.emplace_back(Index(index.x, index.y + 1, index.z));
  }
  if(index.z + 1 < zMax)
  {
    neighbors.emplace_back(Index(index.x, index.y, index.z + 1));
  }

#ifdef DEBUG
  G4cout << "Neighbors of : " << index << G4endl;
  for(const auto& it : neighbors)
  {
    G4cout << it << G4endl;
  }
#endif

  if(neighbors.size() > 6)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "neighbors.size() > 6";
    G4Exception("G4DNAMesh::FindVoxelNeighbors", "G4DNAMesh002",
                FatalErrorInArgument, exceptionDescription);
  }
  return neighbors;
}

G4double G4DNAMesh::GetResolution() const { return fResolution; }

G4DNAMesh::Key G4DNAMesh::GetKey(const G4ThreeVector& position) const
{
  if(!fpBoundingMesh->contains(position))
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "the position: " << position
                         << " is not in the box";
    G4Exception("G4DNAMesh::GetKey", "G4DNAMesh010", FatalErrorInArgument,
                exceptionDescription);
  }

  auto dx = std::floor((position.x() - fpBoundingMesh->Getxlo()) / fResolution);
  auto dy = std::floor((position.y() - fpBoundingMesh->Getylo()) / fResolution);
  auto dz = std::floor((position.z() - fpBoundingMesh->Getzlo()) / fResolution);
  auto xmax =
    std::floor((fpBoundingMesh->Getxhi() - fpBoundingMesh->Getxlo()) / fResolution);
  auto ymax =
    std::floor((fpBoundingMesh->Getyhi() - fpBoundingMesh->Getylo()) / fResolution);
  return dz * ymax * xmax + dy * xmax + dx;
}
