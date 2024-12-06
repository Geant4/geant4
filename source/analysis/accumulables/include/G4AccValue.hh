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

// Class template for accumulable values handled by Geant4 analysis
//
// Author: Ivana Hrivnacova, IJCLab IN2P3/CNRS, 07/09/2015

#ifndef G4AccValue_h
#define G4AccValue_h 1

#include "G4VAccumulable.hh"
#include "G4MergeMode.hh"

#include "globals.hh"

template <typename T>
class G4AccValue : public G4VAccumulable
{
  public:
    G4AccValue(const G4String& name, T initValue,
                  G4MergeMode mergeMode = G4MergeMode::kAddition);
    G4AccValue(T initValue = 0,
                  G4MergeMode mergeMode = G4MergeMode::kAddition);
    G4AccValue(const G4AccValue& rhs);
    G4AccValue(G4AccValue&& rhs) noexcept;
    ~G4AccValue() override = default;

    // Operators
    G4AccValue<T>& operator= (const G4AccValue<T>& rhs);
    G4AccValue<T>& operator=(G4AccValue<T>&& rhs) noexcept;
    G4AccValue<T>& operator+=(const G4AccValue<T>& rhs);
    G4AccValue<T>& operator*=(const G4AccValue<T>& rhs);
    G4AccValue<T>  operator++(int); // postfix increment
    G4AccValue<T>& operator++();    // prefix increment

    G4AccValue<T>& operator= (const T& rhs);
    G4AccValue<T>& operator+=(const T& rhs);
    G4AccValue<T>& operator*=(const T& rhs);

    // Methods
    void Merge(const G4VAccumulable& other) final;
    void Reset() final;
    using G4VAccumulable::Print;
    void Print(G4PrintOptions options = G4PrintOptions()) const final;

    // using G4VAccumulable::SetMergeMode;
    void SetMergeMode(G4MergeMode value) final;

    G4AccType GetType() const final { return G4AccType::kValue; }

    // Get methods
    T  GetValue() const;

  private:
    // Data members
    T  fValue = 0;
    T  fInitValue = 0;
    G4MergeFunction<T> fMergeFunction;
 };

// inline functions

#include "G4AccValue.icc"

#endif
