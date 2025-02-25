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

// Author: Ivana Hrivnacova, IJCLab IN2P3/CNRS, 25/07/2024

#ifndef G4AccType_h
#define G4AccType_h 1

#include "globals.hh"

#include <bitset>

// Enumeration for definition available accummulables

enum class G4AccType {
  kValue,        // G4AccValue<T>
  kArray,        // G4AccArray<T>
  kMap,          // G4AccMap<T>
  kUnorderedMap, // G4AccUnorderedMap<T>
  kVector,       // G4AccumulableVector<T>
  kUser          // User type
};

// TODO: add G4String GetTypeName(G4AccType)

// Helper class for printing
//
class G4PrintOptions {
  public:
    enum Option : std::size_t {
        kName = 0,
        kType = 1,
        kId = 2
    };

    G4PrintOptions() {
        SetDefaults();
    }
    G4PrintOptions(std::initializer_list<Option> options) {
        SetDefaults();
        for (const Option& option : options) {
            fValue.set(option);
        }
    }
    ~G4PrintOptions() = default;

    void Set(Option option, bool value = true) { fValue.set(option, value); }
    bool Has(Option option) const { return fValue[option]; }

  private:
    // methods
    void SetDefaults() {
      fValue.set(kName);
      fValue.set(kType);
    }

    // data members
    static constexpr std::size_t kMax = G4PrintOptions::kId + 1;
    std::bitset<kMax> fValue{0};
};

namespace G4Accumulables
{

// Constant expressions
//
constexpr G4int kInvalidId { -1 };

// Verbose level
//
[[maybe_unused]] static G4int VerboseLevel {1};

}

#endif
