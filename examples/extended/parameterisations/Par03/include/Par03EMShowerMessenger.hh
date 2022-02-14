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
#ifndef PAR03EMSHOWERMESSENGER_HH
#define PAR03EMSHOWERMESSENGER_HH

class Par03EMShowerModel;
class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;

#include "G4UImessenger.hh"

/**
 * @brief Messenger for the example fast simulation model.
 *
 * Allows to set the parameters of the Par03EMShowerModel: parameters of the
 * distributions used in the parametrisation, number of created (same) energy
 * deposits, as well as the maximum depth of an EM shower.
 *
 */

class Par03EMShowerMessenger : public G4UImessenger
{
 public:
  Par03EMShowerMessenger(Par03EMShowerModel* aModel);
  ~Par03EMShowerMessenger();

 public:
  /// Invokes appropriate methods based on the typed command
  virtual void SetNewValue(G4UIcommand* aCommand, G4String aNewValues) final;
  /// Retrieves the current settings
  virtual G4String GetCurrentValue(G4UIcommand* aCommand) final;

 private:
  /// Model to setup
  Par03EMShowerModel* fModel;
  /// Command to set the up a directory for model settings /Par03/fastSim
  G4UIdirectory* fDirectory;
  /// Command printing current settings
  G4UIcmdWithoutParameter* fPrintCmd;
  /// Command to set the sigma parameter of the Gaussian distribution describing
  /// the transverse profile
  G4UIcmdWithADoubleAndUnit* fSigmaCmd;
  /// Command to set the alpha parameter of the Gamma distribution describing
  /// the longitudinal profile
  G4UIcmdWithADouble* fAlphaCmd;
  /// Command to set the beta parameter of the Gamma distribution describing the
  /// longitudinal profile
  G4UIcmdWithADouble* fBetaCmd;
  /// Command to set the number of (same energy) deposits to be created by
  /// fast simulation
  G4UIcmdWithAnInteger* fNbOfHitsCmd;
  /// Command to set the maximum shower depth
  G4UIcmdWithADouble* fLongMaxDepthCmd;
};

#endif /* PAR03EMSHOWERMESSENGER_HH */