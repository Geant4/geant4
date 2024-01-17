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

#ifndef PAR04PARALLELMESSENGER_H
#define PAR04PARALLELMESSENGER_H

#include <G4String.hh>       // for G4String
#include "G4UImessenger.hh"  // for G4UImessenger
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;
class G4UIcommand;
class G4UIdirectory;
class Par04ParallelFullWorld;

/**
 * @brief Parallel messenger.
 *
 * Provides UI commands to setup readout geometry in parallel world (prior to
 * initialization). Number of layers is taken from detector construction.
 *
 */

class Par04ParallelMessenger : public G4UImessenger
{
 public:
  Par04ParallelMessenger(Par04ParallelFullWorld*);
  ~Par04ParallelMessenger();

  /// Invokes appropriate methods based on the typed command
  virtual void SetNewValue(G4UIcommand*, G4String) final;
  /// Retrieves the current settings
  virtual G4String GetCurrentValue(G4UIcommand*) final;

 private:
  /// Parallel world to setup
  Par04ParallelFullWorld* fParallel = nullptr;
  /// Command to set the directory common to all messengers in this example
  /// /Par04
  G4UIdirectory* fExampleDir = nullptr;
  /// Command to set the directory for parallel settings /Par04/parallel
  G4UIdirectory* fParallelDir = nullptr;
  /// Command printing current settings
  G4UIcmdWithoutParameter* fPrintCmd;
  /// Command to set the number of slices
  G4UIcmdWithAnInteger* fNbSlicesCmd = nullptr;
  /// Command to set the number of rows
  G4UIcmdWithAnInteger* fNbRowsCmd = nullptr;
};

#endif
