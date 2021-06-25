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
// G4PhysicsFreeVector class implementation
//
// Authors:
// - 02 Dec. 1995, G.Cosmo: Structure created based on object model
// - 06 Jun. 1996, K.Amako: Implemented the 1st version
// Revisions:
// - 11 Nov. 2000, H.Kurashige: Use STL vector for dataVector and binVector
// --------------------------------------------------------------------

#include "G4PhysicsFreeVector.hh"

// --------------------------------------------------------------------
G4PhysicsFreeVector::G4PhysicsFreeVector(std::size_t length, G4bool spline)
  : G4PhysicsVector(spline)
{
  type          = T_G4PhysicsFreeVector;
  numberOfNodes = length;

  if(0 < length) {
    dataVector.reserve(numberOfNodes);
    binVector.reserve(numberOfNodes);

    for(std::size_t i = 0; i < numberOfNodes; ++i)
    {
      binVector.push_back(0.0);
      dataVector.push_back(0.0);
    }
  }
}

// --------------------------------------------------------------------
G4PhysicsFreeVector::G4PhysicsFreeVector(std::size_t length, G4double emin,
                                         G4double emax, G4bool spline)
  : G4PhysicsFreeVector(length, spline)
{
  edgeMin = emin;
  edgeMax = emax;
}

// --------------------------------------------------------------------
G4PhysicsFreeVector::G4PhysicsFreeVector(const std::vector<G4double>& energies,
                                         const std::vector<G4double>& values,
                                         G4bool spline)
  : G4PhysicsVector(spline)
{
  type          = T_G4PhysicsFreeVector;
  numberOfNodes = energies.size();

  if(numberOfNodes != values.size())
  {
    G4ExceptionDescription ed;
    ed << "The size of energy vector " << numberOfNodes << " != " << values.size();
    G4Exception("G4PhysicsFreeVector constructor: ","glob04", FatalException, ed);
  }

  if(0 < numberOfNodes)
  {
    binVector.reserve(numberOfNodes);
    dataVector.reserve(numberOfNodes);

    for(std::size_t i = 0; i < numberOfNodes; ++i)
    {
      binVector.push_back(energies[i]);
      dataVector.push_back(values[i]);
    }
    edgeMin = binVector[0];
    edgeMax = binVector[numberOfNodes - 1];
  }
}

// --------------------------------------------------------------------
G4PhysicsFreeVector::G4PhysicsFreeVector(const G4double* energies,
                                         const G4double* values,
                                         std::size_t length,
                                         G4bool spline)
  : G4PhysicsVector(spline)
{
  type          = T_G4PhysicsFreeVector;
  numberOfNodes = length;

  if(0 < numberOfNodes) 
  {
    binVector.reserve(numberOfNodes);
    dataVector.reserve(numberOfNodes);

    for(std::size_t i = 0; i < numberOfNodes; ++i)
    {
      binVector.push_back(energies[i]);
      dataVector.push_back(values[i]);
    }
    edgeMin = binVector[0];
    edgeMax = binVector[numberOfNodes - 1];
  }
}

// --------------------------------------------------------------------
G4PhysicsFreeVector::~G4PhysicsFreeVector() 
{}

// --------------------------------------------------------------------
void G4PhysicsFreeVector::PutValues(std::size_t index, G4double e,
                                    G4double value)
{
  if(index >= numberOfNodes)
  {
    PrintPutValueError(index, binVector[numberOfNodes - 1],  e);
  }
  binVector[index]  = e;
  dataVector[index] = value;
  if(index == 0)
  {
    edgeMin = edgeMax = e;
  }
  else if(numberOfNodes == index + 1)
  {
    edgeMax = e;
  }
}

// --------------------------------------------------------------------
void G4PhysicsFreeVector::InsertValues(G4double energy, G4double value)
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
G4double G4PhysicsFreeVector::GetEnergy(G4double aValue)
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
std::size_t G4PhysicsFreeVector::FindValueBinLocation(G4double aValue)
{
  std::size_t bin =
    std::lower_bound(dataVector.cbegin(), dataVector.cend(), aValue) -
    dataVector.cbegin() - 1;
  bin = std::min(bin, numberOfNodes - 2);
  return bin;
}

// --------------------------------------------------------------------
G4double G4PhysicsFreeVector::LinearInterpolationOfEnergy(G4double aValue, 
                                                          std::size_t bin)
{
  G4double res = binVector[bin];
  G4double del = dataVector[bin + 1] - dataVector[bin];
  if(del > 0.0)
  {
    res += (aValue - dataVector[bin]) * (binVector[bin + 1] - res) / del;
  }
  return res;
}
