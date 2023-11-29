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
// Helper namespace field_utils
//
// Description:
//
//  Simple methods to extract vectors from arrays in conventions of
//  the magnetic field integration.

// Author: Dmitry Sorokin, Google Summer of Code 2017
// Supervision: John Apostolakis, CERN
// --------------------------------------------------------------------
#ifndef G4FIELD_UTILS_HH
#define G4FIELD_UTILS_HH

#include "G4FieldTrack.hh"
#include "G4Types.hh"
#include "G4ThreeVector.hh"

namespace field_utils
{

  using State = G4double[G4FieldTrack::ncompSVEC];

  template <unsigned int N>
  using ShortState = G4double[N];
   
  enum class Value3D
  {
    Position = 0,
    Momentum = 3,
    Spin = 9
  };

  enum class Value1D
  {
    KineticEnergy = 6,
    LabTime = 7,
    ProperTime = 8
  };

  template <typename ArrayType>
  G4double getValue(const ArrayType& array, Value1D value);

  template <typename ArrayType>
  G4double getValue2(const ArrayType& array, Value1D value);

  template <typename ArrayType>
  G4double getValue(const ArrayType& array, Value3D value);

  template <typename ArrayType>
  G4double getValue2(const ArrayType& array, Value3D value);

  template <typename ArrayType>
  G4ThreeVector makeVector(const ArrayType& array, Value3D value);

  G4double absoluteError(
    const G4double y[],
    const G4double yerr[],
    G4double hstep);

  G4double relativeError2(
    const G4double y[],
    const G4double yerr[],
    G4double hstep,
    G4double errorTolerance);

  G4double relativeError(
    const G4double y[],
    const G4double yerr[],
    G4double hstep,
    G4double errorTolerance);

  template <typename SourceArray, typename TargetArray>
  void setValue(const SourceArray& src, Value1D value, TargetArray& trg);

  template <typename SourceArray, typename TargetArray, typename ...TargetArrays>
  void setValue(const SourceArray& src, Value1D value,
                TargetArray& trg, TargetArrays&... trgs);

  void copy(G4double dst[], const G4double src[],
            std::size_t size = G4FieldTrack::ncompSVEC);

  G4double inverseCurvatureRadius(G4double particleCharge,
                                  G4double momentum, G4double BField);

  template <typename T>
  T clamp(T value, T lo, T hi);

} // field_utils

#include "G4FieldUtils.icc"

#endif
