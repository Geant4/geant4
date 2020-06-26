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
// G4PhysicsVector
//
// Class description:
//
// A physics vector which has values of energy-loss, cross-section,
// and other physics values of a particle in matter in a given
// range of energy, momentum, etc.
// This class serves as the base class for a vector having various
// energy scale, for example like 'log', 'linear', 'free', etc.

// Authors:
// - 02 Dec. 1995, G.Cosmo: Structure created based on object model
// - 03 Mar. 1996, K.Amako: Implemented the 1st version
// Revisions:
// - 11 Nov. 2000, H.Kurashige: Use STL vector for dataVector and binVector
// - 02 Apr. 2008, A.Bagulya: Added SplineInterpolation() and SetSpline()
// - 19 Jun. 2009, V.Ivanchenko: Removed hidden bin
// - 15 Mar. 2019  M.Novak: added Value method with the known log-energy value
//                          that can avoid the log call in case of log-vectors
// --------------------------------------------------------------------
#ifndef G4PhysicsVector_hh
#define G4PhysicsVector_hh 1

#include <fstream>
#include <iostream>
#include <vector>

#include "G4Log.hh"
#include "G4PhysicsVectorType.hh"
#include "G4ios.hh"
#include "globals.hh"

using G4PVDataVector = std::vector<G4double>;

class G4PhysicsVector
{
 public:
  explicit G4PhysicsVector(G4bool spline = false);
  // Default constructor - vector will be filled via Retrieve() method

  G4PhysicsVector(const G4PhysicsVector&);
  G4PhysicsVector& operator=(const G4PhysicsVector&);
  // Copy constructor and assignment operator

  G4bool operator==(const G4PhysicsVector& right) const;
  G4bool operator!=(const G4PhysicsVector& right) const;
  // Equality operators

  virtual ~G4PhysicsVector();

  G4double Value(G4double theEnergy, std::size_t& lastidx) const;
  // Get the cross-section/energy-loss value corresponding to the
  // given energy. An appropriate interpolation is used to calculate
  // the value. Consumer code gets changed index and may reuse it
  // for the next call to save CPU for bin location.

  inline G4double LogVectorValue(const G4double theEnergy,
                                 const G4double theLogEnergy) const;
  // Same as the Value() method above but specialised for log-vector type.
  // Note, unlike the general Value() method above, this method will work
  // properly only in case of G4PhysicsLogVector-s.

  inline G4double Value(G4double theEnergy) const;
  // Get the cross-section/energy-loss value corresponding to the
  // given energy. An appropriate interpolation is used to calculate
  // the value. This method is kept for backward compatibility reason,
  // it should be used instead of the previous method if bin location
  // cannot be kept thread safe

  inline G4double GetValue(G4double theEnergy, G4bool& isOutRange) const;
  // Obsolete method to get value, 'isOutRange' is not used anymore.
  // This method is kept for the compatibility reason

  inline G4double operator[](const std::size_t index) const;
  // Returns the value for the specified index of the dataVector
  // The boundary check will not be done

  inline G4double operator()(const std::size_t index) const;
  // Returns the value for the specified index of the dataVector
  // The boundary check will not be done

  inline void PutValue(std::size_t index, G4double theValue);
  // Put 'theValue' into the dataVector specified by 'index'.
  // Take note that the 'index' starts from '0'.
  // To fill the vector, need to beforehand construct a vector
  // by the constructor with Emin, Emax, Nbin. 'theValue' should
  // be the cross-section/energy-loss value corresponding to the
  // energy of the index

  virtual void ScaleVector(G4double factorE, G4double factorV);
  // Scale all values of the vector and second derivatives
  // by factorV, energies by vectorE. This method may be applied
  // for example after retrieving a vector from an external file to
  // convert values into Geant4 units

  inline G4double Energy(std::size_t index) const;
  // Returns the value in the energy specified by 'index'
  // of the energy vector. The boundary check will not be done.
  // Use this function when compute cross-section or dEdx
  // before filling the vector by PutValue()

  inline G4double GetMaxEnergy() const;
  // Returns the energy of the last point of the vector

  G4double GetLowEdgeEnergy(std::size_t binNumber) const;
  // Obsolete method
  // Get the energy value at the low edge of the specified bin.
  // Take note that the 'binNumber' starts from '0'.
  // The boundary check will not be done

  inline std::size_t GetVectorLength() const;
  // Get the total length of the vector

  inline std::size_t FindBin(const G4double energy,
                             const std::size_t idx) const;
  // Find low edge index of a bin for given energy.
  // Min value 0, max value VectorLength-1.
  // idx is suggested bin number from user code

  inline std::size_t ComputeLogVectorBin(const G4double logenergy) const;
  // Computes the lower index the energy bin in case of log-vector i.e.
  // in case of vectors with equal bin widths on log-scale

  void FillSecondDerivatives();
  // Initialise second derivatives for Spline keeping
  // 3rd derivative continues - default algorithm.
  // Warning: this method should be called when the vector
  // is already filled

  void ComputeSecDerivatives();
  // Initialise second derivatives for Spline using algorithm
  // which garantee only 1st derivative continues.
  // Warning: this method should be called when the vector
  // is already filled

  void ComputeSecondDerivatives(G4double firstPointDerivative,
                                G4double endPointDerivative);
  // Initialise second derivatives for Spline using
  // user defined 1st derivatives at edge points.
  // Warning: this method should be called when the vector
  // is already filled

  G4double FindLinearEnergy(G4double rand) const;
  // Find energy using linear interpolation for vector
  // filled by cumulative probability function
  // value of rand should be between 0 and 1

  inline G4bool IsFilledVectorExist() const;
  // Is non-empty physics vector already exist?

  inline G4PhysicsVectorType GetType() const;
  // Get physics vector type

  inline void SetSpline(G4bool);
  // Activate/deactivate Spline interpolation

  G4bool Store(std::ofstream& fOut, G4bool ascii = false) const;
  virtual G4bool Retrieve(std::ifstream& fIn, G4bool ascii = false);
  // To store/retrieve persistent data to/from file streams.

  friend std::ostream& operator<<(std::ostream&, const G4PhysicsVector&);
  void DumpValues(G4double unitE = 1.0, G4double unitV = 1.0) const;
  // Print vector

  inline void SetVerboseLevel(G4int value);

 protected:
  void DeleteData();
  void CopyData(const G4PhysicsVector& vec);
  // Internal methods for allowing copy of objects

  void PrintPutValueError(std::size_t index);

  G4PhysicsVectorType type = T_G4PhysicsVector;
  // The type of PhysicsVector (enumerator)

  G4double edgeMin = 0.0;  // Energy of first point
  G4double edgeMax = 0.0;  // Energy of the last point

  G4double invdBin = 0.0;  // 1/Bin width - useful only for fixed binning
  G4double baseBin = 0.0;  // Set this in constructor for performance

  G4int verboseLevel        = 0;
  std::size_t numberOfNodes = 0;

  G4PVDataVector dataVector;     // Vector to keep the crossection/energyloss
  G4PVDataVector binVector;      // Vector to keep energy
  G4PVDataVector secDerivative;  // Vector to keep second derivatives

 private:
  G4bool SplinePossible();

  inline std::size_t FindBinLocation(const G4double theEnergy) const;
  // Find low edge index of a bin for given energy.
  // Min value 0, max value VectorLength-1

  inline G4double Interpolation(const std::size_t idx,
                                const G4double energy) const;

  G4bool useSpline = false;
};

#include "G4PhysicsVector.icc"

#endif
