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
#include <ostream>
#include "G4ITTrackHolder.hh"

std::ostream& operator<<(std::ostream& stream, const G4VDNAMesh::Index& rhs)
{
  stream << "{" << rhs.x << ", " << rhs.y << ", " << rhs.z << "}";
  return stream;
}

G4DNAMesh::Voxel& G4DNAMesh::GetVoxel(const Index& key)
{
  auto iter = fIndexMap.find(key);
  if(iter == fIndexMap.end())
  {
    auto box = GetBoundingBox(key);
    Data mapList;
    G4DNAMesh::Voxel& voxel =
      fVoxelVector.emplace_back(std::make_tuple(key, box, std::move(mapList)));
    fIndexMap[key] = G4int(fVoxelVector.size() - 1);
    return voxel;
  }
  else
  {
    auto index = fIndexMap[key];
    return fVoxelVector[index];
  }
}

G4DNAMesh::G4DNAMesh(const G4DNABoundingBox& boundingBox, G4int pixel)
  : fpBoundingMesh(&boundingBox)
  , fResolution((2 * boundingBox.halfSideLengthInY() / pixel))
{}

G4DNAMesh::~G4DNAMesh() { Reset(); }

G4DNAMesh::Data& G4DNAMesh::GetVoxelMapList(const Index& key)
{
  auto& pVoxel = GetVoxel(key);
  return std::get<2>(pVoxel);
}

void G4DNAMesh::PrintMesh()
{
  G4cout << "*********PrintMesh::Size : " << fVoxelVector.size() << G4endl;
  for(const auto& iter : fVoxelVector)
  {
    auto data = std::get<2>(iter);
    G4cout << "Index : " << std::get<0>(iter)
           << " number of type : " << std::get<2>(iter).size() << G4endl;
    for(const auto& it : data)
    {
      G4cout << "_____________" << it.first->GetName() << " : " << it.second
             << G4endl;
    }
    G4cout << G4endl;
  }
  G4cout << G4endl;
}

G4int G4DNAMesh::GetNumberOfType(G4DNAMesh::MolType type) const
{
  G4int output = 0;

  for(const auto& iter : fVoxelVector)
  {
    auto data = std::get<2>(iter);
    auto it   = data.find(type);
    if(it != data.end())
    {
      output += it->second;
    }
  }
  return output;
}

void G4DNAMesh::PrintVoxel(const Index& index)
{
  G4cout << "*********PrintVoxel::";
  G4cout << " index : " << index
         << " number of type : " << this->GetVoxelMapList(index).size()
         << G4endl;

  for(const auto& it : this->GetVoxelMapList(index))
  {
    G4cout << "_____________" << it.first->GetName() << " : " << it.second
           << G4endl;
  }
  G4cout << G4endl;
}

void G4DNAMesh::InitializeVoxel(const Index& index, Data&& mapList)
{
  auto& pVoxel        = GetVoxel(index);
  std::get<2>(pVoxel) = std::move(mapList);
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

void G4DNAMesh::Reset()
{
  fIndexMap.clear();
  fVoxelVector.clear();
}

const G4DNABoundingBox& G4DNAMesh::GetBoundingBox() const
{
  return *fpBoundingMesh;
}

std::vector<G4DNAMesh::Index>  // array is better ?
G4DNAMesh::FindNeighboringVoxels(const Index& index) const
{
  std::vector<Index> neighbors;
  neighbors.reserve(6);
  auto xMax = (G4int) (std::floor(
    (fpBoundingMesh->Getxhi() - fpBoundingMesh->Getxlo()) / fResolution));
  auto yMax = (G4int) (std::floor(
    (fpBoundingMesh->Getyhi() - fpBoundingMesh->Getylo()) / fResolution));
  auto zMax = (G4int) (std::floor(
    (fpBoundingMesh->Getzhi() - fpBoundingMesh->Getzlo()) / fResolution));

  if(index.x - 1 >= 0)
  {
    neighbors.push_back(Index(index.x - 1, index.y, index.z));
  }
  if(index.y - 1 >= 0)
  {
    neighbors.push_back(Index(index.x, index.y - 1, index.z));
  }
  if(index.z - 1 >= 0)
  {
    neighbors.push_back(Index(index.x, index.y, index.z - 1));
  }
  if(index.x + 1 < xMax)
  {
    neighbors.push_back(Index(index.x + 1, index.y, index.z));
  }
  if(index.y + 1 < yMax)
  {
    neighbors.push_back(Index(index.x, index.y + 1, index.z));
  }
  if(index.z + 1 < zMax)
  {
    neighbors.push_back(Index(index.x, index.y, index.z + 1));
  }

  return neighbors;
}

G4double G4DNAMesh::GetResolution() const { return fResolution; }

G4DNAMesh::Index G4DNAMesh::GetIndex(const G4ThreeVector& position) const
{
  if(!fpBoundingMesh->contains(position))
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "the position: " << position
                         << " is not in the box : " << *fpBoundingMesh;
    G4Exception("G4DNAMesh::GetKey", "G4DNAMesh010", FatalErrorInArgument,
                exceptionDescription);
  }

  G4int dx =
    std::floor((position.x() - fpBoundingMesh->Getxlo()) / fResolution);
  G4int dy =
    std::floor((position.y() - fpBoundingMesh->Getylo()) / fResolution);
  G4int dz =
    std::floor((position.z() - fpBoundingMesh->Getzlo()) / fResolution);
  if(dx < 0 || dy < 0 || dz < 0)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "the old index: " << position
                         << "  to new index : " << Index(dx, dx, dx);
    G4Exception("G4DNAMesh::CheckIndex", "G4DNAMesh015", FatalErrorInArgument,
                exceptionDescription);
  }
  return Index{ dx, dy, dz };
}

G4VDNAMesh::Index G4DNAMesh::ConvertIndex(const Index& index,
                                          const G4int& pixels) const
{
  G4int xmax = std::floor(
    (fpBoundingMesh->Getxhi() - fpBoundingMesh->Getxlo()) / fResolution);
  G4int ymax = std::floor(
    (fpBoundingMesh->Getyhi() - fpBoundingMesh->Getylo()) / fResolution);
  G4int zmax = std::floor(
    (fpBoundingMesh->Getzhi() - fpBoundingMesh->Getzlo()) / fResolution);
  auto dx = (G4int) (index.x * pixels / xmax);
  auto dy = (G4int) (index.y * pixels / ymax);
  auto dz = (G4int) (index.z * pixels / zmax);
  if(dx < 0 || dy < 0 || dz < 0)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "the old index: " << index
                         << "  to new index : " << Index(dx, dx, dx);
    G4Exception("G4DNAMesh::CheckIndex", "G4DNAMesh013", FatalErrorInArgument,
                exceptionDescription);
  }
  return Index{ dx, dy, dz };
}
