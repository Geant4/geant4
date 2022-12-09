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

// Class template for accumulables handled by Geant4 analysis
//
// Author: Ivana Hrivnacova, 07/09/2015  (ivana@ipno.in2p3.fr)

#ifndef G4Accumulable_h
#define G4Accumulable_h 1

#include "G4VAccumulable.hh"
#include "G4MergeMode.hh"

#include "globals.hh"

template <typename T>
class G4Accumulable : public G4VAccumulable
{
  public:
    G4Accumulable(const G4String& name, T initValue,
                G4MergeMode mergeMode = G4MergeMode::kAddition);
    G4Accumulable(T initValue,
                G4MergeMode mergeMode = G4MergeMode::kAddition);
    G4Accumulable(const G4Accumulable& rhs);
    G4Accumulable(G4Accumulable&& rhs) noexcept;
    G4Accumulable() = delete;
    ~G4Accumulable() override = default;

    // Operators
    G4Accumulable<T>& operator= (const G4Accumulable<T>& rhs);
    G4Accumulable<T>& operator=(G4Accumulable<T>&& rhs) noexcept;
    G4Accumulable<T>& operator+=(const G4Accumulable<T>& rhs);
    G4Accumulable<T>& operator*=(const G4Accumulable<T>& rhs);
    G4Accumulable<T>  operator++(int); // postfix increment
    G4Accumulable<T>& operator++();    // prefix increment

    G4Accumulable<T>& operator= (const T& rhs);
    G4Accumulable<T>& operator+=(const T& rhs);
    G4Accumulable<T>& operator*=(const T& rhs);

    // Methods
    void Merge(const G4VAccumulable& other) final;
    void Reset() final;

    // Get methods
    T  GetValue() const;
    G4MergeMode GetMergeMode() const;

  private:
    // Data members
    T  fValue;
    T  fInitValue;
    G4MergeMode  fMergeMode;
    G4MergeFunction<T> fMergeFunction;
 };

// inline functions

#include "G4Accumulable.icc"

#endif
