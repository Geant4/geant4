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

// Class template for accumulable arrays handled by Geant4 analysis
//
// Author: Ivana Hrivnacova, IJCLab IN2P3/CNRS, 18/07/2024

#ifndef G4AccArray_h
#define G4AccArray_h 1

#include "G4VAccumulable.hh"
#include "G4MergeMode.hh"

#include "globals.hh"

#include <array>

template <class T,  std::size_t N>
class G4AccArray : public G4VAccumulable
{
  public:
    // Default constructor (1)
    // Constructs the container with the default init value
    G4AccArray(const G4String& name = "",
               const T& value = 0,
               G4MergeMode mergeMode = G4MergeMode::kAddition);

    // Constructor (2)
    // Constructs the container with the provided initializer list
    // and defaults for all other members
    template <typename... Args>
    G4AccArray(Args&&... args);

    // // Constructor (3)
    // // Constructs the container with the provided name and initializer list
    // // and defaults for the other members
    template <typename First, typename... Args>
    G4AccArray(const First& firstArg, Args&&... args);

    // Copy constructor
    G4AccArray(const G4AccArray& rhs) = default;
    // Move constructor
    G4AccArray(G4AccArray&& rhs) = default;

    // Destructor
    ~G4AccArray() override = default;

    // operators
    inline T& operator[](typename std::array<T,N>::size_type i) { return fArray[i]; }

    // Methods
    void Merge(const G4VAccumulable& other) final;
    void Reset() final;
    void Print(G4PrintOptions options = G4PrintOptions()) const final;
    void SetMergeMode(G4MergeMode value) final;

    // Get methods
    G4AccType GetType() const final { return G4AccType::kArray; }
    std::array<T,N>& GetArray();
    const std::array<T,N>& GetArray() const;

  private:
    // Data members
    std::array<T,N> fArray {};
    T fInitValue = 0;
    G4MergeFunction<T> fMergeFunction;
};

// inline functions

#include "G4AccArray.icc"

#endif
