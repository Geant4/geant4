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

class G4PhysicsVector
{
public:
  // Default constructor - vector will be filled via Retrieve() method
  // Free vector may be filled via InsertValue(..) method
  explicit G4PhysicsVector(G4bool spline = false);

  // Copy constructor and assignment operator
  G4PhysicsVector(const G4PhysicsVector&) = default;
  G4PhysicsVector& operator=(const G4PhysicsVector&) = default;

  // not used operators
  G4PhysicsVector(const G4PhysicsVector&&) = delete;
  G4PhysicsVector& operator=(const G4PhysicsVector&&) = delete;
  G4bool operator==(const G4PhysicsVector& right) const = delete;
  G4bool operator!=(const G4PhysicsVector& right) const = delete;

  virtual ~G4PhysicsVector() = default;

  // Get the cross-section/energy-loss value corresponding to the
  // given energy. An appropriate interpolation is used to calculate
  // the value. Consumer code gets changed index and may reuse it
  // for the next call to save CPU for bin location.
  inline G4double Value(const G4double energy, std::size_t& lastidx) const;

  // Get the cross-section/energy-loss value corresponding to the
  // given energy. An appropriate interpolation is used to calculate
  // the value. This method should be used if bin location cannot be 
  // kept in the user code.
  inline G4double Value(const G4double energy) const;

  // Obsolete method to get value, 'isOutRange' is not used anymore.
  // This method is kept for the compatibility reason
  inline G4double GetValue(const G4double energy, G4bool& isOutRange) const;

  // Same as the Value() method above but specialised for log-vector type.
  // Note, unlike the general Value() method above, this method will work
  // properly only for G4PhysicsLogVector.
  inline G4double LogVectorValue(const G4double energy,
                                 const G4double theLogEnergy) const;

  // Returns the value for the specified index of the dataVector
  // The boundary check will not be done
  inline G4double operator[](const std::size_t index) const;
  inline G4double operator()(const std::size_t index) const;

  // Put data into the vector at 'index' position.
  // Take note that the 'index' starts from '0'.
  // It is assumed that energies are already filled.
  inline void PutValue(const std::size_t index, const G4double value);

  // Returns the value in the energy specified by 'index'
  // of the energy vector. The boundary check will not be done.
  // Use this when compute cross-section, dEdx, or other value
  // before filling the vector by PutValue().
  inline G4double Energy(const std::size_t index) const;
  inline G4double GetLowEdgeEnergy(const std::size_t index) const;

  // Returns the energy of the first and the last point of the vector.
  inline G4double GetMinEnergy() const;
  inline G4double GetMaxEnergy() const;

  // Returns the data of the first and the last point of the vector.
  // If the vector is empty returns zeros.
  inline G4double GetMinValue() const;
  inline G4double GetMaxValue() const;

  // Get the total length of the vector
  inline std::size_t GetVectorLength() const;

  // Computes the lower index the energy bin in case of log-vector i.e.
  // in case of vectors with equal bin widths on log-scale
  // Note, that no check on the boundary is performed
  inline std::size_t ComputeLogVectorBin(const G4double logenergy) const;

  // Get physics vector type.
  inline G4PhysicsVectorType GetType() const;

  // True if using spline interpolation.
  inline G4bool GetSpline() const;

  // Define verbosity level.
  inline void SetVerboseLevel(G4int value);

  // Find energy using linear interpolation for vector
  // filled by cumulative probability function.
  // Assuming that vector is already filled.
  inline G4double FindLinearEnergy(const G4double rand) const;

  // Find low edge index of a bin for given energy.
  // Min value 0, max value idxmax.
  std::size_t FindBin(const G4double energy, std::size_t idx) const;

  // Scale all values of the vector by factorV, energies by vectorE.
  // AFter this method FillSecondDerivatives(...) should be called. 
  // This method may be applied for example after retrieving a vector 
  // from an external file to convert values into Geant4 units.
  void ScaleVector(const G4double factorE, const G4double factorV);

  // This method should be called when the vector is fully filled 
  // There are 3 types of second derivative computations:
  //    fSplineSimple -     2d derivative continues
  //    fSplineBase -       3d derivative continues (the default)
  //    fSplineFixedEdges - 3d derivatives continues, 1st and last 
  //                        derivatives are fixed  
  void FillSecondDerivatives(const G4SplineType = G4SplineType::Base,
                             const G4double dir1 = 0.0,
                             const G4double dir2 = 0.0);

  // This method can be applied if both energy and data values 
  // grow monotonically, for example, if in this vector a 
  // cumulative probability density function is stored. 
  G4double GetEnergy(const G4double value) const;

  // To store/retrieve persistent data to/from file streams.
  G4bool Store(std::ofstream& fOut, G4bool ascii = false) const;
  G4bool Retrieve(std::ifstream& fIn, G4bool ascii = false);

  // Print vector
  friend std::ostream& operator<<(std::ostream&, const G4PhysicsVector&);
  void DumpValues(G4double unitE = 1.0, G4double unitV = 1.0) const;

protected:

  // The default implements a free vector initialisation.
  virtual void Initialise();

  void PrintPutValueError(std::size_t index, G4double value, 
                          const G4String& text);

private:

  void ComputeSecDerivative0();
  void ComputeSecDerivative1();
  void ComputeSecDerivative2(const G4double firstPointDerivative,
                             const G4double endPointDerivative);
  // Internal methods for computing of spline coeffitients

  // Linear or spline interpolation.
  inline G4double Interpolation(const std::size_t idx,
                                const G4double energy) const;

  // Assuming (edgeMin <= energy <= edgeMax).
  inline std::size_t GetBin(const G4double energy) const;

protected:

  G4double edgeMin = 0.0;  // Energy of first point
  G4double edgeMax = 0.0;  // Energy of the last point

  G4double invdBin = 0.0;  // 1/Bin width for linear and log vectors
  G4double logemin = 0.0;  // used only for log vector

  G4int verboseLevel = 0;
  std::size_t idxmax = 0;
  std::size_t numberOfNodes = 0;

  G4PhysicsVectorType type = T_G4PhysicsFreeVector;
  // The type of PhysicsVector (enumerator)

  std::vector<G4double> binVector;      // energy
  std::vector<G4double> dataVector;     // crossection/energyloss
  std::vector<G4double> secDerivative;  // second derivatives

private:

  G4bool useSpline = false;
};

#include "G4PhysicsVector.icc"

#endif
