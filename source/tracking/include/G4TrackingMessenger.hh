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
// G4TrackingMessenger
//
// Class description:
//
// This is a messenger class to interface and exchange information
// between tracking/stepping and UI.

// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Makoto  Asai   (e-mail: asai@slac.stanford.edu)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//---------------------------------------------------------------
#ifndef G4TrackingMessenger_hh
#define G4TrackingMessenger_hh 1

#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;
class G4TrackingManager;
class G4SteppingManager;

class G4TrackingMessenger : public G4UImessenger
{
 public:
  G4TrackingMessenger(G4TrackingManager* trMan);
  ~G4TrackingMessenger() override;
  void SetNewValue(G4UIcommand* command, G4String newValues) override;
  G4String GetCurrentValue(G4UIcommand* command) override;

 private:
  G4TrackingManager* trackingManager = nullptr;
  G4SteppingManager* steppingManager = nullptr;

  // commands

  G4UIdirectory* TrackingDirectory = nullptr;
  G4UIcmdWithoutParameter* AbortCmd = nullptr;
  G4UIcmdWithoutParameter* ResumeCmd = nullptr;
  G4UIcmdWithAnInteger* StoreTrajectoryCmd = nullptr;
  G4UIcmdWithAnInteger* VerboseCmd = nullptr;
};

#endif
