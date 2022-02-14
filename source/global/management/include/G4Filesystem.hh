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
#ifndef G4FILESYSTEM_HH
#define G4FILESYSTEM_HH
/**
 * @namespace G4fs
 * Namespace adapter for `filesystem` in `std::` or `std::experimental::`
 *
 * Developers should use this namespace to scope into the C++ filesystem
 * library, e.g.
 *
 * ```cpp
 * // Don't use std:: directly...
 * auto p = std::filesystem::path{ argv[0] };
 *
 * // ... use the alias
 * auto p = G4fs::path{ argv[0] };
 * ```
 *
 * This allows compatibility with C++ standard libraries that only provide
 * `filesystem` in the `std::experimental` namespace
 *
 */
#if __has_include(<filesystem>)
#  include <filesystem>
namespace G4fs = std::filesystem;
#elif __has_include(<experimental/filesystem>)
#  include <experimental/filesystem>
namespace G4fs = std::experimental::filesystem;
#endif

#endif /* G4FILESYSTEM_HH */