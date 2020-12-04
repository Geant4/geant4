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
// G4PhysicsOrderedFreeVector class implementation
//
// Author: Juliet Armstrong (TRIUMF), 13 August 1996
// Revisions:
// - 11.11.2000, H.Kurashige: use STL vector for dataVector and binVector
// - 19.06.2009, V.Ivanchenko: removed hidden bin
// --------------------------------------------------------------------

#include "G4PhysicsOrderedFreeVector.hh"

// --------------------------------------------------------------------
G4PhysicsOrderedFreeVector::G4PhysicsOrderedFreeVector()
  : G4PhysicsVector()
{
  type = T_G4PhysicsOrderedFreeVector;
}

// --------------------------------------------------------------------
G4PhysicsOrderedFreeVector::G4PhysicsOrderedFreeVector(G4double* Energies,
                                                       G4double* Values,
                                                       std::size_t VectorLength)
  : G4PhysicsVector()
{
  type = T_G4PhysicsOrderedFreeVector;

  dataVector.reserve(VectorLength);
  binVector.reserve(VectorLength);

  for(std::size_t i = 0; i < VectorLength; ++i)
  {
    InsertValues(Energies[i], Values[i]);
  }
}

// --------------------------------------------------------------------
G4PhysicsOrderedFreeVector::G4PhysicsOrderedFreeVector(
  const std::vector<G4double>& Energies, const std::vector<G4double>& Values)
  : G4PhysicsVector()
{
  if(Energies.size() != Values.size())
  {
    G4ExceptionDescription ed;
    ed << "The sizes of the two std::vector arguments must be the same";
    G4Exception("G4PhysicsOrderedFreeVector::G4PhysicsOrderedFreeVector()",
                "glob04", FatalException, ed);
  }

  type = T_G4PhysicsOrderedFreeVector;

  dataVector.reserve(Energies.size());
  binVector.reserve(Energies.size());

  for(std::size_t i = 0; i < Energies.size(); ++i)
  {
    InsertValues(Energies[i], Values[i]);
  }
}

// --------------------------------------------------------------------
G4PhysicsOrderedFreeVector::~G4PhysicsOrderedFreeVector() {}

// --------------------------------------------------------------------
void G4PhysicsOrderedFreeVector::InsertValues(G4double energy, G4double value)
{
  auto binLoc = std::lower_bound(binVector.cbegin(), binVector.cend(), energy);

  std::size_t binIdx = binLoc - binVector.cbegin();  // Iterator difference!

  auto dataLoc = dataVector.cbegin() + binIdx;

  binVector.insert(binLoc, energy);
  dataVector.insert(dataLoc, value);

  ++numberOfNodes;
  edgeMin = binVector.front();
  edgeMax = binVector.back();
}

// --------------------------------------------------------------------
G4double G4PhysicsOrderedFreeVector::GetEnergy(G4double aValue)
{
  G4double e;
  if(aValue <= GetMinValue())
  {
    e = edgeMin;
  }
  else if(aValue >= GetMaxValue())
  {
    e = edgeMax;
  }
  else
  {
    std::size_t closestBin = FindValueBinLocation(aValue);
    e                      = LinearInterpolationOfEnergy(aValue, closestBin);
  }
  return e;
}

// --------------------------------------------------------------------
std::size_t G4PhysicsOrderedFreeVector::FindValueBinLocation(G4double aValue)
{
  std::size_t bin =
    std::lower_bound(dataVector.cbegin(), dataVector.cend(), aValue) -
    dataVector.cbegin() - 1;
  bin = std::min(bin, numberOfNodes - 2);
  return bin;
}

// --------------------------------------------------------------------
G4double G4PhysicsOrderedFreeVector::LinearInterpolationOfEnergy(
  G4double aValue, std::size_t bin)
{
  G4double res = binVector[bin];
  G4double del = dataVector[bin + 1] - dataVector[bin];
  if(del > 0.0)
  {
    res += (aValue - dataVector[bin]) * (binVector[bin + 1] - res) / del;
  }
  return res;
}
