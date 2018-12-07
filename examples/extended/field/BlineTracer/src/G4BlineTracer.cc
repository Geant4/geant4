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
/// \file field/BlineTracer/src/G4BlineTracer.cc
/// \brief Implementation of the G4BlineTracer class
//
//
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
#include "G4MagIntegratorDriver.hh"
#include "G4CashKarpRKF45.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ChordFinder.hh"
#include "G4SystemOfUnits.hh"

//////////////////////////////////////////////////////////////////

G4BlineTracer::G4BlineTracer()
{
  fMessenger = new G4BlineTracerMessenger(this);
  fSteppingAction = new G4BlineSteppingAction(this) ;
  fEventAction = new G4BlineEventAction(this);
  fPrimaryGeneratorAction = new G4BlinePrimaryGeneratorAction();
  fMaxTrackingStep =1000.*m;
  fWas_ResetChordFinders_already_called=false;
}

///////////////////////////////////////////////////////////////////////

G4BlineTracer::~G4BlineTracer()
{
  delete fMessenger;
  delete fSteppingAction;
  delete fEventAction; 
  delete fPrimaryGeneratorAction;
  for (size_t i=0; i< fVecEquationOfMotion.size();i++)
  {
    if (fVecEquationOfMotion[i]) delete fVecEquationOfMotion[i];
    if (fVecChordFinders[i]) delete fVecChordFinders[i];
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
  if (!fWas_ResetChordFinders_already_called)
  {
    ResetChordFinders();
    fWas_ResetChordFinders_already_called=true;
  }

  // Replace the user action by the ad-hoc actions for Blines
  
  G4RunManager* theRunManager =  G4RunManager::GetRunManager();
  G4UserRunAction* user_run_action =
    (G4UserRunAction*)theRunManager->GetUserRunAction();
  theRunManager->SetUserAction(this);

  G4UserSteppingAction* user_stepping_action =
    (G4UserSteppingAction*)theRunManager->GetUserSteppingAction();
  theRunManager->SetUserAction(fSteppingAction);
  
  G4VUserPrimaryGeneratorAction* userPrimaryAction =
    (G4VUserPrimaryGeneratorAction*)theRunManager->GetUserPrimaryGeneratorAction();
  if (userPrimaryAction) 
    fPrimaryGeneratorAction->SetUserPrimaryAction(userPrimaryAction);
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

  // replace the user defined chordfinder by the element of fVecChordFinders  
  
  std::vector<G4ChordFinder*> user_chord_finders;
  std::vector<G4double> user_largest_acceptable_step;
  for (size_t i=0;i<fVecChordFinders.size();i++)
  {
    user_largest_acceptable_step.push_back(-1.);
    if (fVecChordFinders[i])
    {
      user_chord_finders.push_back(fVecFieldManagers[i]->GetChordFinder());
      fVecChordFinders[i]->SetDeltaChord(user_chord_finders[i]->GetDeltaChord());
      fVecFieldManagers[i]->SetChordFinder(fVecChordFinders[i]);
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
          ->SetLargestAcceptableStep(fMaxTrackingStep);

  // Start the integration of n_of_lines different magnetic field lines

  for (G4int il=0; il<n_of_lines;il++)
  {
    // for each magnetic field line we integrate once backward and once
    // forward from the same starting point

    // backward integration

    for (size_t i=0; i< fVecEquationOfMotion.size();i++)
    {
      if (fVecEquationOfMotion[i]) 
        fVecEquationOfMotion[i]->SetBackwardDirectionOfIntegration(true);
    }
    theRunManager->BeamOn(1);

    // forward integration

    for (size_t i=0; i < fVecEquationOfMotion.size();i++)
    {
      if (fVecEquationOfMotion[i]) 
        fVecEquationOfMotion[i]->SetBackwardDirectionOfIntegration(false);
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
  theRunManager->SetUserAction(userPrimaryAction);
  theRunManager->SetUserAction(user_stepping_action);
  theRunManager->SetUserAction(user_tracking_action);
  theRunManager->SetUserAction(user_stacking_action);

  // set user defined chord finders and largest acceptable step

  for (size_t i=0;i<fVecFieldManagers.size();i++)
  {
    if (user_chord_finders[i])
      fVecFieldManagers[i]->SetChordFinder(user_chord_finders[i]);
  } 
}

////////////////////////////////////////////////////////////////

/*
G4bool G4BlineTracer::CheckMagneticFields()
{
  // Check FieldManagers

  G4TransportationManager* tmanager =
    G4TransportationManager::GetTransportationManager();

  if (fVecFieldManagers[0] != tmanager->GetFieldManager())
    return false;
  if (fVecMagneticFields[0] != tmanager->GetFieldManager()->GetDetectorField())
    return false;
  G4LogicalVolumeStore* theVolumeStore = G4LogicalVolumeStore::GetInstance();  
   
  std::vector<G4FieldManagers*> LogicalVolumeFields;
  size_t j=0;
  for (size_t i=0; i<theVolumeStore.size();i++)
  {
    if (theVolumeStore[i]->GetFieldManager())
    {
      j++;
      if (j >= fVecFieldManagers.size()) return false;
      if (fVecFieldManagers[j] != theVolumeStore[i]->GetFieldManager())
        return false;
      if (fVecMagneticFields[j] !=
          theVolumeStore[i]->GetFieldManager()->GetDetectorField())
        return false;
    }
  }
  if (j<fVecFieldManagers.size()) return false;

 return true;
}
*/

////////////////////////////////////////////////////////////////

void G4BlineTracer::ResetChordFinders()
{
  for (size_t i=0; i<fVecEquationOfMotion.size();i++)
  {
    delete fVecEquationOfMotion[i];
    delete fVecChordFinders[i];
  } 

  fVecChordFinders.clear();
  fVecFieldManagers.clear();
  fVecMagneticFields.clear();
  fVecEquationOfMotion.clear();

  // global field

  fVecChordFinders.push_back(0);
  fVecMagneticFields.push_back(0);
  fVecEquationOfMotion.push_back(0);
  fVecFieldManagers.push_back(G4TransportationManager::GetTransportationManager()
                             ->GetFieldManager());
  if (fVecFieldManagers[0])
  {
    fVecMagneticFields[0] =
      (G4MagneticField*) fVecFieldManagers[0]->GetDetectorField();
    if (fVecMagneticFields[0])
    {
      fVecEquationOfMotion[0] = new G4BlineEquation(fVecMagneticFields[0]);
      G4CashKarpRKF45* pStepper = new G4CashKarpRKF45(fVecEquationOfMotion[0]);
      G4MagInt_Driver* pIntgrDriver =
        new G4MagInt_Driver(0.01*mm,pStepper,pStepper->GetNumberOfVariables());
      fVecChordFinders[0] = new G4ChordFinder(pIntgrDriver);
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
      fVecFieldManagers.push_back(((*theVolumeStore)[i])->GetFieldManager());
      fVecMagneticFields.push_back((G4MagneticField*)
                                   fVecFieldManagers[j]->GetDetectorField());
      fVecEquationOfMotion.push_back(0);
      fVecChordFinders.push_back(0);
      if (fVecMagneticFields[j])
      {
        fVecEquationOfMotion[j]= new G4BlineEquation(fVecMagneticFields[j]);
        G4CashKarpRKF45* pStepper = new G4CashKarpRKF45(fVecEquationOfMotion[j]);
        G4MagInt_Driver* pIntgrDriver =
          new G4MagInt_Driver(.01*mm,pStepper,pStepper->GetNumberOfVariables());
        fVecChordFinders[j] = new G4ChordFinder(pIntgrDriver);
      } 
    }
  }             
}
