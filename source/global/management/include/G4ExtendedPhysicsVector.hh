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
// G4ExtendedPhysicsVector
//
// Class description:
//
// A data scructure which includes G4PhysicsVector and also data on
// partial x-sections
//
// Author: V.Ivanchenko 09.09.2025
//
// --------------------------------------------------------------------
#ifndef G4ExtendedPhysicsVector_hh
#define G4ExtendedPhysicsVector_hh 1

#include <vector>

#include "globals.hh"
#include "G4PhysicsVector.hh"

class G4ExtendedPhysicsVector
{
public:
 
  explicit G4ExtendedPhysicsVector(G4PhysicsVector*, G4int nxsec = 0);
  virtual ~G4ExtendedPhysicsVector();

  // Copy constructor and assignment operator
  G4ExtendedPhysicsVector(const G4ExtendedPhysicsVector&) = default;
  G4ExtendedPhysicsVector& operator=(const G4ExtendedPhysicsVector&) = default;

  // not used operators
  G4ExtendedPhysicsVector(const G4ExtendedPhysicsVector&&) = delete;
  G4ExtendedPhysicsVector& operator=(const G4ExtendedPhysicsVector&&) = delete;
  G4bool operator==(const G4ExtendedPhysicsVector& right) const = delete;

  // get G4PhysicsVector with total cross section
  inline const G4PhysicsVector* GetPhysicsVector() const;
  
  // Get the value from the total x-section using interpolation defined
  // in the G4PhysicsVector object. Consumer code gets changed
  // index and may reuse it for the next call to save CPU for bin location.
  inline G4double Value(const G4double energy, std::size_t& lastidx) const;
  
  // Get the cross-section/energy-loss value corresponding to the
  // given energy using log-log interpolation. Consumer code gets changed
  // index and may reuse it for the next call to save CPU for bin location.
  // Spline is not used in this method.
  G4double LogLogValue(const G4double energy, std::size_t& lastidx) const;

  // This method may be applied only once.
  // Force length of data using std::vector::resize() with the
  // the default value 0; partial cross section vector is resized
  // only if the number of partial x-sections is above zero.
  void SetDataLength(G4int dlength);
  
  // Filled partial cross sections by energy index (0 <= idx < numberOfNodes)
  // It is assumed that the array y does not have negative values.
  void PutPartialXSData(const std::size_t idx, const G4double* y);

  // Sample partial reaction channel for given energy, rand - uniform
  // random number, lastidx is the cache allowing to reduce CPU for energy
  // bin location. Data are stored as float but computations are in double.
  G4int SampleReactionChannel(const G4double energy, const G4double rand,
                              std::size_t& lastidx) const;
  
  // randomly select partial cross section using Log-Log interpolation
  G4int SampleReactionChannelLogLog(const G4double energy, const G4double rand,
                                    std::size_t& lastidx) const;
  
  // Print data
  void DumpValues(G4double unitE = 1.0, G4double unitV = 1.0) const;

private:

  G4int verboseLevel = 0;
  G4int nPartialXS = 0;
  std::size_t idxmax = 0;
  std::size_t numberOfNodes = 0;

  // The partial cumulative x-sections normalized to the sum    
  std::vector<std::vector<G4float>* >* dataPartialXS = nullptr;
  G4PhysicsVector* totalData = nullptr;
};

inline const G4PhysicsVector* G4ExtendedPhysicsVector::GetPhysicsVector() const
{
  return totalData;
}

inline G4double G4ExtendedPhysicsVector::Value(const G4double e, std::size_t& idx) const
{
  return totalData->Value(e, idx);
}

#endif
