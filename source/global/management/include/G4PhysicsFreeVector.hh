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
// G4PhysicsFreeVector
//
// Class description:
//
// A physics vector which has values of energy-loss, cross-section,
// and other physics values of a particle in matter in a given
// range of the energy, momentum, etc. The scale of energy/momentum
// bins is in free, i.e. it is NOT need to be linear or log. Only
// restriction is that bin values alway have to increase from
// a lower bin to a higher bin. This is necessary for the binary
// search to work correctly.

// Authors:
// - 02 Dec. 1995, G.Cosmo: Structure created based on object model
// - 06 Jun. 1996, K.Amako: Implemented the 1st version
// Revisions:
// - 11 Nov. 2000, H.Kurashige: Use STL vector for dataVector and binVector
// - 04 Feb. 2021, V.Ivanchenko moved implementation of all free vectors 
//                 to this class
// --------------------------------------------------------------------
#ifndef G4PhysicsFreeVector_hh
#define G4PhysicsFreeVector_hh 1

#include "G4PhysicsVector.hh"
#include "globals.hh"
#include <vector>

class G4PhysicsFreeVector : public G4PhysicsVector
{
public:

  // Constructors of a free vector without filling of data.
  // The vector will be filled using PutValues(..), Retrieve(..), or 
  // InsertValues(..) methods. 
  // If length > 0 energy and data vectors are initialized with zeros
  explicit G4PhysicsFreeVector(G4bool spline = false);
  explicit G4PhysicsFreeVector(G4int length);
  explicit G4PhysicsFreeVector(std::size_t length, G4bool spline = false);

  // Obsolete constructor - emin and emax are not used 
  explicit G4PhysicsFreeVector(std::size_t length, G4double emin,
                               G4double emax, G4bool spline = false);

  // The vector is filled in these constructor;
  // 'energies' and 'values' need to have the same vector length;
  // 'energies' assumed to increase, it is allowed to have consequtive
  // equal energies
  explicit G4PhysicsFreeVector(const std::vector<G4double>& energies,
                               const std::vector<G4double>& values,
                               G4bool spline = false);
  explicit G4PhysicsFreeVector(const G4double* energies, const G4double* values,
                               std::size_t length, G4bool spline = false);

  ~G4PhysicsFreeVector() override = default;

  // Filling of the vector with the check on index and energy
  void PutValues(const std::size_t index, 
                 const G4double energy, const G4double value);

  // Insert extra pair of energy and value
  // If energy coincide with previously added energy then
  // this new pair is added after 
  void InsertValues(const G4double energy, const G4double value);
 
  // Obsolete method
  inline void PutValue(const std::size_t index, 
                       const G4double e, const G4double value);
};

inline void G4PhysicsFreeVector::PutValue(const std::size_t index, 
                                          const G4double e,
                                          const G4double value)
{
  PutValues(index, e, value);
}

#endif
