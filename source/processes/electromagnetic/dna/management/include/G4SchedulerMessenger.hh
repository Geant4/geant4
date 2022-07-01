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
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr)
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// Author: Mathieu Karamitros

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial.
// If you use Geant4-DNA chemistry and you publish papers about your software,
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508

#ifndef G4ITSTEPPINGMESSENGER_H
#define G4ITSTEPPINGMESSENGER_H

class G4Scheduler;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;
class G4UIcommand;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
#include "G4UImessenger.hh"
#include "globals.hh"

class G4SchedulerMessenger : public G4UImessenger
{
 public:
  explicit G4SchedulerMessenger(G4Scheduler* runMgr);
  ~G4SchedulerMessenger();
  void SetNewValue(G4UIcommand* command, G4String newValues);
  G4String GetCurrentValue(G4UIcommand* command);

 private:
  G4Scheduler* fScheduler;

 private:  // commands
  std::unique_ptr<G4UIdirectory> fITDirectory;
  std::unique_ptr<G4UIcmdWithADoubleAndUnit> fEndTime;
  std::unique_ptr<G4UIcmdWithADoubleAndUnit> fTimeTolerance;
  std::unique_ptr<G4UIcmdWithAnInteger> fVerboseCmd;
  std::unique_ptr<G4UIcmdWithAnInteger> fMaxStepNumber;
  std::unique_ptr<G4UIcmdWithoutParameter> fInitCmd;
  std::unique_ptr<G4UIcmdWithoutParameter> fProcessCmd;
  std::unique_ptr<G4UIcmdWithAnInteger> fMaxNULLTimeSteps;
  std::unique_ptr<G4UIcmdWithoutParameter> fWhyDoYouStop;
  std::unique_ptr<G4UIcmdWithABool> fUseDefaultTimeSteps;
  std::unique_ptr<G4UIcmdWithABool> fResetScavenger;
};

#endif  // G4ITSTEPPINGMESSENGER_H
