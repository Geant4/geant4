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
// G4PhysicsOrderedFreeVector
//
// Class description:
//
// A physics ordered free vector inherits from G4PhysicsVector which
// has values of energy-loss, cross-section, and other physics values
// of a particle in matter in a given range of energy, momentum, etc.).
// In addition, the ordered free vector provides a method for
// the user to insert energy/value pairs in sequence. Methods to
// retrieve the Max and Min energies and values from the vector are
// also provided.

// Author: Juliet Armstrong (TRIUMF), 13 August 1996
// Revisions:
// - 11.11.2000, H.Kurashige: use STL vector for dataVector and binVector
// --------------------------------------------------------------------
#ifndef G4PhysicsOrderedFreeVector_hh
#define G4PhysicsOrderedFreeVector_hh 1

#include "G4PhysicsVector.hh"

class G4PhysicsOrderedFreeVector : public G4PhysicsVector
{
 public:
  G4PhysicsOrderedFreeVector();
  // The vector will be filled from extern file using Retrieve()
  // or InsertValues() methods

  G4PhysicsOrderedFreeVector(const std::vector<G4double>& Energies,
                             const std::vector<G4double>& Values);
  G4PhysicsOrderedFreeVector(G4double* Energies, G4double* Values,
                             std::size_t VectorLength);
  // The vector is filled in this constructor.
  // 'Energies' and 'Values' need to have the same vector length
  // 'Energies' assumed to be ordered

  virtual ~G4PhysicsOrderedFreeVector();

  void InsertValues(G4double energy, G4double value);

  G4double GetEnergy(G4double aValue);

  inline G4double GetMaxValue();

  inline G4double GetMinValue();

  inline G4double GetMaxLowEdgeEnergy();

  inline G4double GetMinLowEdgeEnergy();

 private:
  std::size_t FindValueBinLocation(G4double aValue);

  G4double LinearInterpolationOfEnergy(G4double aValue, std::size_t locBin);
};

// -----------------------------
// Inline methods implementation
// -----------------------------

inline G4double G4PhysicsOrderedFreeVector::GetMaxValue()
{
  return dataVector.back();
}

inline G4double G4PhysicsOrderedFreeVector::GetMinValue()
{
  return dataVector.front();
}

inline G4double G4PhysicsOrderedFreeVector::GetMaxLowEdgeEnergy()
{
  return binVector.back();
}

inline G4double G4PhysicsOrderedFreeVector::GetMinLowEdgeEnergy()
{
  return binVector.front();
}

#endif
