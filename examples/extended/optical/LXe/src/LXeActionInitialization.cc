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
//
/// \file LXeActionInitialization.cc
/// \brief Implementation of the LXeActionInitialization class

#include "LXeActionInitialization.hh"

#include "LXePrimaryGeneratorAction.hh"
#include "LXeDetectorConstruction.hh"

#include "LXeRunAction.hh"
#include "LXeEventAction.hh"
#include "LXeTrackingAction.hh"
#include "LXeSteppingAction.hh"
#include "LXeStackingAction.hh"
#include "LXeSteppingVerbose.hh"

#include "LXeRecorderBase.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

<<<<<<< HEAD
LXeActionInitialization::LXeActionInitialization(LXeRecorderBase* recorder)
 : G4VUserActionInitialization(), fRecorder(recorder)
=======
LXeActionInitialization::LXeActionInitialization(
   const LXeDetectorConstruction* det)
  : G4VUserActionInitialization(), fDetector(det)
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeActionInitialization::~LXeActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeActionInitialization::BuildForMaster() const
{
  SetUserAction(new LXeRunAction(fRecorder));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeActionInitialization::Build() const
{
  SetUserAction(new LXePrimaryGeneratorAction());

<<<<<<< HEAD
  SetUserAction(new LXeStackingAction());

  SetUserAction(new LXeRunAction(fRecorder));
  SetUserAction(new LXeEventAction(fRecorder));
  SetUserAction(new LXeTrackingAction(fRecorder));
  SetUserAction(new LXeSteppingAction(fRecorder));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
=======
  LXeEventAction* eventAction = new LXeEventAction(fDetector);
  SetUserAction(eventAction);
  SetUserAction(new LXeStackingAction(eventAction));
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c

G4VSteppingVerbose* LXeActionInitialization::InitializeSteppingVerbose() const
{
  return new LXeSteppingVerbose();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
