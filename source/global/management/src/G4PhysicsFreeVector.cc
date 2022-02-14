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
// - 25 Aug. 2021, V.Ivanchenko updated for Geant4 11.0 
// --------------------------------------------------------------------

#include "G4PhysicsFreeVector.hh"

// --------------------------------------------------------------------
G4PhysicsFreeVector::G4PhysicsFreeVector(G4bool spline)
  : G4PhysicsVector(spline)
{}

// --------------------------------------------------------------------
G4PhysicsFreeVector::G4PhysicsFreeVector(G4int length)
  : G4PhysicsFreeVector((std::size_t)length, false)
{}

// --------------------------------------------------------------------
G4PhysicsFreeVector::G4PhysicsFreeVector(std::size_t length, G4bool spline)
  : G4PhysicsVector(spline)
{
  numberOfNodes = length;

  if(0 < length) {
    binVector.resize(numberOfNodes, 0.0);
    dataVector.resize(numberOfNodes, 0.0);
  }
  Initialise();
}

// --------------------------------------------------------------------
G4PhysicsFreeVector::G4PhysicsFreeVector(std::size_t length, G4double,
                                         G4double, G4bool spline)
  : G4PhysicsFreeVector(length, spline)
{}

// --------------------------------------------------------------------
G4PhysicsFreeVector::G4PhysicsFreeVector(const std::vector<G4double>& energies,
                                         const std::vector<G4double>& values,
                                         G4bool spline)
  : G4PhysicsVector(spline)
{
  numberOfNodes = energies.size();

  if(numberOfNodes != values.size())
  {
    G4ExceptionDescription ed;
    ed << "The size of energy vector " << numberOfNodes << " != " << values.size();
    G4Exception("G4PhysicsFreeVector constructor: ","glob04", FatalException, ed);
  }

  binVector = energies;
  dataVector = values;
  Initialise();
}

// --------------------------------------------------------------------
G4PhysicsFreeVector::G4PhysicsFreeVector(const G4double* energies,
                                         const G4double* values,
                                         std::size_t length,
                                         G4bool spline)
  : G4PhysicsVector(spline)
{
  numberOfNodes = length;

  if(0 < numberOfNodes) 
  {
    binVector.resize(numberOfNodes);
    dataVector.resize(numberOfNodes);

    for(std::size_t i = 0; i < numberOfNodes; ++i)
    {
      binVector[i] = energies[i];
      dataVector[i] = values[i];
    }
  }
  Initialise();
}

// --------------------------------------------------------------------
void G4PhysicsFreeVector::PutValues(const std::size_t index, 
                                    const G4double e,
                                    const G4double value)
{
  if(index >= numberOfNodes)
  {
    PrintPutValueError(index, value, "G4PhysicsFreeVector::PutValues ");
    return;
  }
  binVector[index]  = e;
  dataVector[index] = value;
  if(index == 0)
  {
    edgeMin = e;
  }
  else if(numberOfNodes == index + 1)
  {
    edgeMax = e;
  }
}

// --------------------------------------------------------------------
void G4PhysicsFreeVector::InsertValues(const G4double energy, 
                                       const G4double value)
{
  auto binLoc = std::lower_bound(binVector.cbegin(), binVector.cend(), energy);
  auto dataLoc = dataVector.cbegin();
  dataLoc += binLoc - binVector.cbegin(); 

  binVector.insert(binLoc, energy);
  dataVector.insert(dataLoc, value);

  ++numberOfNodes;
  Initialise();
}

// --------------------------------------------------------------------
