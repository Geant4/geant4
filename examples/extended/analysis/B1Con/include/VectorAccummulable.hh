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
//
/// \file VectorAccummulable.hh
/// \brief Definition of the VectorAccummulable class
///
/// This class defines an accummulable of a vector<T> type
/// that merges the vectors filled on workers in a single vector.

#ifndef VectorAccumulable_h
#define VectorAccumulable_h 1

#include "G4VAccumulable.hh"
#include "globals.hh"
#include <vector>

template <typename T>
class VectorAccumulable : public G4VAccumulable
{
  public:
    VectorAccumulable() = default;
    ~VectorAccumulable() override = default;

    void AddValue(T value);
    const std::vector<T>& GetVector() const;

    void Merge(const G4VAccumulable& other) override;
    void Reset() override;

  private:
    std::vector<T> fTVector;
};

// inline functions

template <typename T>
inline void VectorAccumulable<T>::AddValue(T value) { 
  fTVector.push_back(value); 
}

template <typename T>
inline const std::vector<T>& VectorAccumulable<T>::GetVector() const
{
  return fTVector;
}

template <typename T>
inline void VectorAccumulable<T>::Merge(const G4VAccumulable& other) { 
  for (const auto& value : static_cast<const VectorAccumulable<T>&>(other).fTVector )  {
    fTVector.push_back(value);
  }
}

template <typename T>
inline void VectorAccumulable<T>::Reset() { 
  fTVector.clear();
}

#endif
