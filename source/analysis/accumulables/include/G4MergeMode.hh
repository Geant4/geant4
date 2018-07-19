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
// $Id$

// Author: Ivana Hrivnacova, 04/07/2012  (ivana@ipno.in2p3.fr)

#ifndef G4MergeMode_h
#define G4MergeMode_h 1

#include "globals.hh"

#include <functional>
#include <cmath>

// Enumeration for definition available binning schemes

enum class G4MergeMode {
  kAddition,        // "Or" if boolean type
  kMultiplication,  // "And" if boolean type
  kMaximum,         // "Or" if boolean type
  kMinimum          // "And" if boolean type
};  

// Generic merge function
template <typename T>
using G4MergeFunction = std::function<T(const T&, const T&)>;

// Utility functions
namespace  G4Accumulables {

G4MergeMode GetMergeMode(const G4String& mergeModeName);

template <typename T>
G4MergeFunction<T> GetMergeFunction(G4MergeMode mergeMode)
{
  switch ( mergeMode ) {
    case G4MergeMode::kAddition:
      // return std::bind([](const T& x, const T& y) { return x + y; });
      return [](const T& x, const T& y) { return x + y; };

    case G4MergeMode::kMultiplication:
      return [](const T& x, const T& y) { return x * y; };
      
    case G4MergeMode::kMaximum:
      // return std::bind([](const T& x, const T& y) { return std::max(x,y);});
      return [](const T& x, const T& y) { return std::max(x,y); };

    case G4MergeMode::kMinimum:
      // return std::bind([](const T& x, const T& y) { return std::min(x,y);});
      return [](const T& x, const T& y) { return std::min(x,y); };
  }
}

template <G4bool>
G4MergeFunction<G4bool> GetMergeFunction(G4MergeMode mergeMode)
{
  switch ( mergeMode ) {
    case G4MergeMode::kAddition:
    case G4MergeMode::kMaximum:
      // return std::bind([](const T& x, const T& y) { return x + y; });
      return [](const G4bool& x, const G4bool& y) { return x || y; };

    case G4MergeMode::kMultiplication:
    case G4MergeMode::kMinimum:
      return [](const G4bool& x, const G4bool& y) { return x && y; };
  }
}

}

#endif
