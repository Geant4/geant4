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
/// \file TrackingMessenger.cc
/// \brief Implementation of the TrackingMessenger class
//
// $Id: TrackingMessenger.cc 98257 2016-07-04 17:39:46Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingMessenger.hh"
#include "TrackingAction.hh"

#include "G4UIcmdWithABool.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingMessenger::TrackingMessenger(TrackingAction* trackA)
:G4UImessenger(),
 fTrackingAction(trackA),fTrackingCmd(0),fTimeWindowCmd(0)
{
  fTrackingCmd = new G4UIcmdWithABool("/rdecay01/fullChain",this);
  fTrackingCmd->SetGuidance("allow full decay chain");
  fTrackingCmd->SetParameterName("flag",true);
  fTrackingCmd->SetDefaultValue(true);

 fTimeWindowCmd = new G4UIcommand("/rdecay01/timeWindow",this);
 fTimeWindowCmd->SetGuidance("set time window to survey activity per nuclide");
 fTimeWindowCmd->SetGuidance("  t1, delta_t");
 //    
 G4UIparameter* t1Prm = new G4UIparameter("t1",'d',false);
 t1Prm->SetGuidance("time_1");
 t1Prm->SetParameterRange("t1>0.");
 fTimeWindowCmd->SetParameter(t1Prm);
 //
 G4UIparameter* unit1Prm = new G4UIparameter("unit1",'s',false);
 unit1Prm->SetGuidance("unit of t1");
 G4String unitList = G4UIcommand::UnitsList(G4UIcommand::CategoryOf("second"));
 unit1Prm->SetParameterCandidates(unitList);
 fTimeWindowCmd->SetParameter(unit1Prm);
 //
 G4UIparameter* dtPrm = new G4UIparameter("dt",'d',false);
 dtPrm->SetGuidance("delta_t");
 dtPrm->SetParameterRange("dt>0.");
 fTimeWindowCmd->SetParameter(dtPrm);
 //
 G4UIparameter* unit2Prm = new G4UIparameter("unit2",'s',false);
 unit2Prm->SetGuidance("unit of dt");
 unit2Prm->SetParameterCandidates(unitList);
 fTimeWindowCmd->SetParameter(unit2Prm);
 //
 fTimeWindowCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingMessenger::~TrackingMessenger()
{
  delete fTrackingCmd;
  delete fTimeWindowCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if (command == fTrackingCmd)
   { fTrackingAction->SetFullChain(fTrackingCmd->GetNewBoolValue(newValue));}

  if (command == fTimeWindowCmd)
   {
     G4double t1; G4double dt;
     G4String unt1, unt2;
     std::istringstream is(newValue);
     is >> t1 >> unt1 >> dt >> unt2;
     t1 *= G4UIcommand::ValueOf(unt1);
     dt *= G4UIcommand::ValueOf(unt2);
     fTrackingAction->SetTimeWindow (t1, dt);
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
