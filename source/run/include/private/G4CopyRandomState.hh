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
/// @file G4CopyRandomState.hh
/// @brief Helper function for copying random state files in G4run
/// @author Ben Morgan
/// @date 2023-06-05 

#ifndef G4CopyRandomState_hh
#define G4CopyRandomState_hh

#include "G4Exception.hh"
#include "G4Filesystem.hh"
#include "G4String.hh"
#include "G4Types.hh"

/// Return true if file `source` is successfully copied to `dest`
/// Convert any thrown exception to JustWarning G4Exception
inline G4bool G4CopyRandomState(const G4fs::path& source, const G4fs::path& dest,
                              const G4String& callsite)
{
  try {
    G4fs::copy_file(source, dest);
  }
  catch (G4fs::filesystem_error const& ex) {
    G4ExceptionDescription ed;
    ed << "Failed to copy " << ex.path1() << " to " << ex.path2() << " , error:\n"
       << "  code   : " << ex.code().value() << '\n'
       << "  message: " << ex.code().message() << '\n';
    G4Exception(callsite, "UnableToCopyRndmStateFile", JustWarning, ed);
    return false;
  }
  return true;
}

#endif  // G4CopyRandomState_hh
