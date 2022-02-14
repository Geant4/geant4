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
// G4Physics2DVector
//
// Class description:
//
// A 2-dimentional vector with linear interpolation.

// Author: Vladimir Ivanchenko, 25.09.2011
// --------------------------------------------------------------------
#ifndef G4Physics2DVector_hh
#define G4Physics2DVector_hh 1

#include <fstream>
#include <iostream>
#include <vector>

#include "G4PhysicsVectorType.hh"
#include "G4ios.hh"
#include "globals.hh"

using G4PV2DDataVector = std::vector<G4double>;

class G4Physics2DVector
{
 public:
  G4Physics2DVector();
  // Vector will be filled via Retrieve method

  explicit G4Physics2DVector(std::size_t nx, std::size_t ny);
  // Vector will be filled via Put methods

  G4Physics2DVector(const G4Physics2DVector&);
  G4Physics2DVector& operator=(const G4Physics2DVector&);
  // Copy constructor and assignment operator

  G4bool operator==(const G4Physics2DVector& right) const = delete;
  G4bool operator!=(const G4Physics2DVector& right) const = delete;

  ~G4Physics2DVector();
  // Destructor

  G4double Value(G4double x, G4double y, std::size_t& lastidx,
                 std::size_t& lastidy) const;
  G4double Value(G4double x, G4double y) const;
  // Main method to interpolate 2D vector
  // Consumer class should provide initial values of lastidx and lastidy

  inline void PutX(std::size_t idx, G4double value);
  inline void PutY(std::size_t idy, G4double value);
  inline void PutValue(std::size_t idx, std::size_t idy, G4double value);
  void PutVectors(const std::vector<G4double>& vecX,
                  const std::vector<G4double>& vecY);
  // Methods to fill vector
  // Take note that the 'index' starts from '0'

  void ScaleVector(G4double factor);
  // Scale all values of the vector by factor.
  // This method may be applied for example after Retrieve a vector
  // from an external file to convert values into Geant4 units

  G4double FindLinearX(G4double rand, G4double y, std::size_t& lastidy) const;
  inline G4double FindLinearX(G4double rand, G4double y) const;
  // Find Y using linear interpolation for Y-vector filled by cumulative
  // probability function value of rand should be between 0 and 1

  inline G4double GetX(std::size_t index) const;
  inline G4double GetY(std::size_t index) const;
  inline G4double GetValue(std::size_t idx, std::size_t idy) const;
  // Returns simply the values of the vector by index
  // of the energy vector. The boundary check will not be done

  inline std::size_t FindBinLocationX(const G4double x,
                                      const std::size_t lastidx) const;
  inline std::size_t FindBinLocationY(const G4double y,
                                      const std::size_t lastidy) const;
  // Find the bin# in which theEnergy belongs. Starting from 0

  inline std::size_t GetLengthX() const;
  inline std::size_t GetLengthY() const;
  // Get the lengths of the vector

  inline G4PhysicsVectorType GetType() const;
  // Get physics vector type

  inline void SetBicubicInterpolation(G4bool);
  // Activate/deactivate bicubic interpolation

  void Store(std::ofstream& fOut) const;
  G4bool Retrieve(std::ifstream& fIn);
  // To store/retrieve persistent data to/from file streams

  inline void SetVerboseLevel(G4int value);

 protected:
  void PrepareVectors();

  void ClearVectors();

  void CopyData(const G4Physics2DVector& vec);

  G4double BicubicInterpolation(const G4double x, const G4double y,
                                const std::size_t idx,
                                const std::size_t idy) const;
  // Bicubic interpolation of 2D vector

  inline std::size_t FindBin(const G4double z, const G4PV2DDataVector&,
                             const std::size_t idz,
                             const std::size_t idzmax) const;

 private:
  G4double InterpolateLinearX(G4PV2DDataVector& v, G4double rand) const;

  inline G4double DerivativeX(std::size_t idx, std::size_t idy,
                              G4double fac) const;
  inline G4double DerivativeY(std::size_t idx, std::size_t idy,
                              G4double fac) const;
  inline G4double DerivativeXY(std::size_t idx, std::size_t idy,
                               G4double fac) const;
  // computation of derivatives

  G4PhysicsVectorType type = T_G4PhysicsFreeVector;
  // The type of PhysicsVector (enumerator)

  std::size_t numberOfXNodes = 0;
  std::size_t numberOfYNodes = 0;

  G4PV2DDataVector xVector;
  G4PV2DDataVector yVector;
  std::vector<G4PV2DDataVector*> value;

  G4int verboseLevel = 0;
  G4bool useBicubic  = false;
};

#include "G4Physics2DVector.icc"

#endif
