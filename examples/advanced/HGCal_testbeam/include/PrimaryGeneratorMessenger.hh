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
#ifndef PRIMARYGENERATORMMESSENGER_HH
#define PRIMARYGENERATORMMESSENGER_HH

#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class PrimaryGeneratorAction;

/**
 * @brief Primary generator messenger
 *
 * Defines UI commands to set up the primary generator.
 *
 * By default particle gun is used as the primary generator. It can be
 * controlled with standard UI commands (/gun/).
 *
 * "/HGCalTestbeam/generator/momentumSpread <VALUE>" to change constant
 * particle energy to Gaussian distribution with sigma expressed in units of the
 * initial energy (e.g. value 0.05 means sigma of 0.05 * E).
 * By default it equals to 0 and constant energy value is used.
 * "/HGCalTestbeam/generator/beamSpread <none/Gaussian/flat>" to define type of
 * beam position spread. By default none is used.
 * "/HGCalTestbeam/generator/beamSpreadX <SIZE>" to define size of beam spread
 * along x axis. It is sigma of a Gaussian distribution, or half-width of a
 * flat distribution.
 * "/HGCalTestbeam/generator/beamSpreadY <SIZE>" to define size of beam spread
 * along y axis. It is sigma of a Gaussian distribution, or half-width of a
 * flat distribution.
 * "/HGCalTestbeam/generator/beamZ0 <POSITION>" to define beam position along z
 * axis. By default edge of the world volume is used.
 *
 * If installation was done with ROOT package (CMake was able to locate it),
 * an additional option of input read from the ROOT file is enabled.
 * It can be activated with "/HGCalTestbeam/generator/readInputFile true".
 * "/HGCalTestbeam/generator/pathInputFile <FILE>" sets the path to the input
 * file.
 * "/HGCalTestbeam/generator/startFromEvent <N>" allows to start simulation from
 * Nth event.
 * Please note that in current implementation input from file needs to be
 * executed in a non-multithreaded mode (or with 1 thread).
 *
 */

class PrimaryGeneratorMessenger : public G4UImessenger {
public:
  explicit PrimaryGeneratorMessenger(PrimaryGeneratorAction *aPrimaryGeneratorAction);
  ~PrimaryGeneratorMessenger();

public:
  void SetNewValue(G4UIcommand *command, G4String newValues);
  G4String GetCurrentValue(G4UIcommand *command);

private:
  /// Pointer to the primmary generator
  PrimaryGeneratorAction *fPrimaryGenerator;

private:
  /// Directory for UI commands
  G4UIdirectory *fDirectory;
#ifdef WITHROOT
  /// Command specyfing if primary event should be read from file
  G4UIcmdWithABool *fReadInputCmd;
  /// Command to set the path to the input file
  G4UIcmdWithAString *fPathInputCmd;
  /// Command to set ID of the first event to be read from the file
  G4UIcmdWithAnInteger *fStartFromEventCmd;
#endif
  /// Command to set the sigma of the Gaussian momentum spread
  G4UIcmdWithADouble *fMomentumSpreadCmd;
  /// Command to set the type of transverse beam position spread
  G4UIcmdWithAString *fBeamSpreadTypeCmd;
  /// Command to set the size of beam position spread along X axis
  G4UIcmdWithADoubleAndUnit *fBeamSpreadXCmd;
  /// Command to set the size of beam position spread along Y axis
  G4UIcmdWithADoubleAndUnit *fBeamSpreadYCmd;
  /// Command to set the initial beam position size along Z axis
  G4UIcmdWithADoubleAndUnit *fBeamZ0Cmd;
};

#endif /*PRIMARYGENERATORMMESSENGER_HH */