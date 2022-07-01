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
// G4ProfilerMessenger
//
// Class description
//
// This class is the messenger of the class which maintain the profiling
// controls (located in global/management/include/G4Profiler.hh).
// In general, it forwards to the argument parser provided by timemory.

// Author: J.Madsen, 12 November 2020
// --------------------------------------------------------------------
#ifndef G4ProfilerMessenger_hh
#define G4ProfilerMessenger_hh 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4Profiler.hh"

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithABool;

// --------------------------------------------------------------------

class G4ProfilerMessenger : public G4UImessenger
{
 public:
  G4ProfilerMessenger();
  ~G4ProfilerMessenger() override;

  // no copy but move is fine
  G4ProfilerMessenger(const G4ProfilerMessenger&) = delete;
  G4ProfilerMessenger(G4ProfilerMessenger&&)      = default;

  // no copy but move is fine
  G4ProfilerMessenger& operator=(const G4ProfilerMessenger&) = delete;
  G4ProfilerMessenger& operator=(G4ProfilerMessenger&&) = default;

  void SetNewValue(G4UIcommand*, G4String) override;
  // G4String GetCurrentValue(G4UIcommand* command);

 private:
  using stringcmd_pair  = std::pair<G4UIcmdWithAString*, std::string>;
  using boolcmd_pair    = std::pair<G4UIcmdWithABool*, std::string>;
  using directory_array = std::array<G4UIdirectory*, G4ProfileType::TypeEnd>;
  using stringcmd_array = std::array<stringcmd_pair, G4ProfileType::TypeEnd>;
  using boolcmd_array   = std::array<boolcmd_pair, G4ProfileType::TypeEnd>;
  using boolcmd_vector  = std::vector<boolcmd_pair>;

  G4UIdirectory* profileDirectory;
  G4UIdirectory* profileOutputDirectory;
  directory_array profileTypeDirs;

  boolcmd_array profileEnableCmds;
  boolcmd_vector profileGeneralCmds;
  stringcmd_array profileCompCmds;
};

#endif
