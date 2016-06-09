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
// $Id: G4BlineTracer.cc,v 1.3 2006-06-29 17:15:16 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// --------------------------------------------------------------------
//
// G4BlineTracer implementation
//
// --------------------------------------------------------------------
// Author: Laurent Desorgher (desorgher@phim.unibe.ch)
//         Created - 2003-10-06
// --------------------------------------------------------------------

#include "G4BlineTracer.hh"
#include "G4BlineTracerMessenger.hh"
#include "G4BlinePrimaryGeneratorAction.hh"
#include "G4BlineEventAction.hh"
#include "G4BlineSteppingAction.hh"
#include "G4BlineEquation.hh"

#include "G4RunManager.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4PropagatorInField.hh"
#include "G4CashKarpRKF45.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ChordFinder.hh"

//////////////////////////////////////////////////////////////////

G4BlineTracer::G4BlineTracer()
{
  fMessenger = new G4BlineTracerMessenger(this);
  fSteppingAction = new G4BlineSteppingAction(this) ;
  fEventAction = new G4BlineEventAction(this);
  fPrimaryGeneratorAction = new G4BlinePrimaryGeneratorAction();
  MaxTrackingStep =1000.*m;
  was_ResetChordFinders_already_called=false;
}

///////////////////////////////////////////////////////////////////////

G4BlineTracer::~G4BlineTracer()
{
  delete fMessenger;
  delete fSteppingAction;
  delete fEventAction; 
  delete fPrimaryGeneratorAction;
  for (size_t i=0; i< vecEquationOfMotion.size();i++)
  {
    if (vecEquationOfMotion[i]) delete vecEquationOfMotion[i];
    if (vecChordFinders[i]) delete vecChordFinders[i];
  }
}  

////////////////////////////////////////////////////////////////////

void G4BlineTracer::BeginOfRunAction(const G4Run*)
{
}  

///////////////////////////////////////////////////////////////////////

void G4BlineTracer::EndOfRunAction(const G4Run*)
{
}

////////////////////////////////////////////////////////////////

void G4BlineTracer::ComputeBlines(G4int n_of_lines)
{
  //the first time ResetChordFinders should be called
  //
  if (!was_ResetChordFinders_already_called)
  {
    ResetChordFinders();
    was_ResetChordFinders_already_called=true;
  }

  // Replace the user action by the ad-hoc actions for Blines
  
  G4RunManager* theRunManager =  G4RunManager::GetRunManager();
  G4UserRunAction* user_run_action =
    (G4UserRunAction*)theRunManager->GetUserRunAction();
  theRunManager->SetUserAction(this);

  G4UserSteppingAction* user_stepping_action =
    (G4UserSteppingAction*)theRunManager->GetUserSteppingAction();
  theRunManager->SetUserAction(fSteppingAction);
  
  G4VUserPrimaryGeneratorAction* fUserPrimaryAction =
    (G4VUserPrimaryGeneratorAction*)theRunManager->GetUserPrimaryGeneratorAction();
  if (fUserPrimaryAction) 
    fPrimaryGeneratorAction->SetUserPrimaryAction(fUserPrimaryAction);
  theRunManager->SetUserAction(fPrimaryGeneratorAction);

  G4UserEventAction* user_event_action =
    (G4UserEventAction*)theRunManager->GetUserEventAction();
  theRunManager->SetUserAction(fEventAction);
 
  G4UserTrackingAction* user_tracking_action =
    (G4UserTrackingAction*)theRunManager->GetUserTrackingAction();
  G4UserTrackingAction* aNullTrackingAction = 0;
  theRunManager->SetUserAction(aNullTrackingAction);

  G4UserStackingAction* user_stacking_action = 
    (G4UserStackingAction*)theRunManager->GetUserStackingAction();
  G4UserStackingAction* aNullStackingAction = 0;
  theRunManager->SetUserAction(aNullStackingAction);

  // replace the user defined chordfinder by the element of vecChordFinders  
  
  std::vector<G4ChordFinder*> user_chord_finders;
  std::vector<G4double> user_largest_acceptable_step;
  for (size_t i=0;i<vecChordFinders.size();i++)
  {
    user_largest_acceptable_step.push_back(-1.);
    if (vecChordFinders[i])
    {
      user_chord_finders.push_back(vecFieldManagers[i]->GetChordFinder());
      vecChordFinders[i]->SetDeltaChord(user_chord_finders[i]->GetDeltaChord());
      vecFieldManagers[i]->SetChordFinder(vecChordFinders[i]);
    }
    else user_chord_finders.push_back(0); 
  }

  // I have tried to use the smooth line filter ability but I could not obtain 
  // a smooth trajectory in the G4TrajectoryContainer after an event
  // Another solution for obtaining a smooth trajectory is to limit
  // the LargestAcceptableStep in the G4PropagatorInField object.
  // This is the solution I used. 

  // Old solution:
  // G4TransportationManager::GetTransportationManager()
  //     ->GetPropagatorInField()->SetTrajectoryFilter(fTrajectoryFilter);

  // New solution:
  // set the largest_acceptable_step to max_step:length

  G4TransportationManager* tmanager =
    G4TransportationManager::GetTransportationManager();
  G4double previous_largest_acceptable_step =
    tmanager->GetPropagatorInField()->GetLargestAcceptableStep();

  tmanager->GetPropagatorInField()
          ->SetLargestAcceptableStep(MaxTrackingStep);

  // Start the integration of n_of_lines different magnetic field lines

  for (G4int i=0; i<n_of_lines;i++)
  {
    // for each magnetic field line we integrate once backward and once
    // forward from the same starting point

    // backward integration

    for (size_t i=0; i< vecEquationOfMotion.size();i++)
    {
      if (vecEquationOfMotion[i]) 
        vecEquationOfMotion[i]->SetBackwardDirectionOfIntegration(true);
    }
    theRunManager->BeamOn(1);

    // forward integration

    for (size_t i=0; i < vecEquationOfMotion.size();i++)
    {
      if (vecEquationOfMotion[i]) 
        vecEquationOfMotion[i]->SetBackwardDirectionOfIntegration(false);
    }
    theRunManager->BeamOn(1);
  }

  // Remove trajectory filter to PropagatorInField
  // It was for old solution when using smooth trajectory filter

  // tmanager->GetPropagatorInField()->SetTrajectoryFilter(0);

  // back to User defined actions and other parameters
  // -------------------------------------------------

  tmanager->GetPropagatorInField()
          ->SetLargestAcceptableStep(previous_largest_acceptable_step);

  // return to User actions

  theRunManager->SetUserAction(user_run_action);
  theRunManager->SetUserAction(user_event_action);
  theRunManager->SetUserAction(fUserPrimaryAction);
  theRunManager->SetUserAction(user_stepping_action);
  theRunManager->SetUserAction(user_tracking_action);
  theRunManager->SetUserAction(user_stacking_action);

  // set user defined chord finders and largest acceptable step

  for (size_t i=0;i<vecFieldManagers.size();i++)
  {
    if (user_chord_finders[i])
      vecFieldManagers[i]->SetChordFinder(user_chord_finders[i]);
  } 
}

////////////////////////////////////////////////////////////////

/*
G4bool G4BlineTracer::CheckMagneticFields()
{
  // Check FieldManagers

  G4TransportationManager* tmanager =
    G4TransportationManager::GetTransportationManager();

  if (vecFieldManagers[0] != tmanager->GetFieldManager())
    return false;
  if (vecMagneticFields[0] != tmanager->GetFieldManager()->GetDetectorField())
    return false;
  G4LogicalVolumeStore* theVolumeStore = G4LogicalVolumeStore::GetInstance();  
   
  std::vector<G4FieldManagers*> LogicalVolumeFields;
  size_t j=0;
  for (size_t i=0; i<theVolumeStore.size();i++)
  {
    if (theVolumeStore[i]->GetFieldManager())
    {
      j++;
      if (j >= vecFieldManagers.size()) return false;
      if (vecFieldManagers[j] != theVolumeStore[i]->GetFieldManager())
        return false;
      if (vecMagneticFields[j] !=
          theVolumeStore[i]->GetFieldManager()->GetDetectorField())
        return false;
    }
  }
  if (j<vecFieldManagers.size()) return false;

 return true;
}
*/

////////////////////////////////////////////////////////////////

void G4BlineTracer::ResetChordFinders()
{
  for (size_t i=0; i<vecEquationOfMotion.size();i++)
  {
    delete vecEquationOfMotion[i];
    delete vecChordFinders[i];
  } 

  vecChordFinders.clear();
  vecFieldManagers.clear();
  vecMagneticFields.clear();
  vecEquationOfMotion.clear();

  // global field

  vecChordFinders.push_back(0);
  vecMagneticFields.push_back(0);
  vecEquationOfMotion.push_back(0);
  vecFieldManagers.push_back(G4TransportationManager::GetTransportationManager()
                             ->GetFieldManager());
  if (vecFieldManagers[0])
  {
    vecMagneticFields[0] =
      (G4MagneticField*) vecFieldManagers[0]->GetDetectorField();
    if (vecMagneticFields[0])
    {
      vecEquationOfMotion[0] = new G4BlineEquation(vecMagneticFields[0]);
      G4CashKarpRKF45* pStepper = new G4CashKarpRKF45(vecEquationOfMotion[0]);
      G4MagInt_Driver* pIntgrDriver =
        new G4MagInt_Driver(0.01*mm,pStepper,pStepper->GetNumberOfVariables());
      vecChordFinders[0] = new G4ChordFinder(pIntgrDriver);
    }
  } 

  // local fields   

  G4LogicalVolumeStore* theVolumeStore = G4LogicalVolumeStore::GetInstance();

  size_t j=0;
  for (size_t i=0; i<theVolumeStore->size();i++)
  {
    if ((*theVolumeStore)[i]->GetFieldManager())
    {
      j++;
      vecFieldManagers.push_back(((*theVolumeStore)[i])->GetFieldManager());
      vecMagneticFields.push_back((G4MagneticField*)
                                   vecFieldManagers[j]->GetDetectorField());
      vecEquationOfMotion.push_back(0);
      vecChordFinders.push_back(0);
      if (vecMagneticFields[j])
      {
        vecEquationOfMotion[j]= new G4BlineEquation(vecMagneticFields[j]);
        G4CashKarpRKF45* pStepper = new G4CashKarpRKF45(vecEquationOfMotion[j]);
        G4MagInt_Driver* pIntgrDriver =
          new G4MagInt_Driver(.01*mm,pStepper,pStepper->GetNumberOfVariables());
        vecChordFinders[j] = new G4ChordFinder(pIntgrDriver);
      } 
    }
  }             
}
