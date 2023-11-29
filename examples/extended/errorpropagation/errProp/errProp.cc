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
/// \file errProp/errProp.cc
/// \brief Main program of the errorpropagation example
//
// ------------------------------------------------------------
//      GEANT 4 example main
// ------------------------------------------------------------
//
// History:
// - Created:   P. Arce   May 2007
//

#include "ExErrorDetectorConstruction.hh"
#include "G4SteppingVerbose.hh"

#include "G4ErrorPropagator.hh"
#include "G4ErrorPropagatorData.hh"
#include "G4ErrorPropagatorManager.hh"
#include "G4ErrorPlaneSurfaceTarget.hh"
#include "G4ErrorCylSurfaceTarget.hh"
#include "G4ErrorGeomVolumeTarget.hh"
#include "G4ErrorTrackLengthTarget.hh"
#include "G4ErrorFreeTrajState.hh"

#include "G4UImanager.hh"
#include "G4SystemOfUnits.hh"

void Initialize();
G4ErrorTarget* BuildTarget( G4int iTarget );
void ProcessEvent( G4int iProp, size_t nEv );
void Finalize();

G4ErrorTarget* theTarget;
G4ErrorMode theG4ErrorMode;

//-------------------------------------------------------------
int main() 
{

  Initialize();

  //----- Choose propagation mode
  // 0: propagate until target, all steps in one go
  // 1: propagate until target, returning control to the user at each step
  G4int iProp = 0;
  char* prop = std::getenv("G4ERROR_PROP");
  if( prop ) {
    if( G4String(prop) == G4String("UNTIL_TARGET") ){
      iProp = 0;
    } else if ( G4String(prop) == G4String("STEP_BY_STEP") ) {
      iProp = 1;
    } else {
      G4Exception("exG4eReco","Fatal error in Argument",
        FatalErrorInArgument,
        G4String("Variable G4ERROR_PROP = " + G4String(prop) + 
                 "   It must be: UNTIL_TARGET or STEP_BY_STEP").c_str());
    }
  } else {
    G4Exception("exG4eReco","Fatal error in Argument",
      JustWarning,
      "Variable G4ERROR_PROP not defined, taking it = UNTIL_TARGET");
  } 

  size_t nEvents = 3;
  for( size_t ii = 0; ii < nEvents; ii++ ){
    ProcessEvent( iProp, ii );
  }

  Finalize();

}


//-------------------------------------------------------------
void Initialize() 
{
  G4VSteppingVerbose::SetInstance(new G4SteppingVerbose);

  // Initialize the GEANT4e manager 
  G4ErrorPropagatorManager* g4emgr 
    = G4ErrorPropagatorManager::GetErrorPropagatorManager();
  G4ErrorPropagatorData* g4edata 
    = G4ErrorPropagatorData::GetErrorPropagatorData();

  g4emgr->SetUserInitialization(new ExErrorDetectorConstruction); 

  G4UImanager::GetUIpointer()->ApplyCommand("/exerror/setField -10. kilogauss");

  g4emgr->InitGeant4e();

  G4UImanager::GetUIpointer()->ApplyCommand("/control/verbose 1");
  G4UImanager::GetUIpointer()->ApplyCommand("/tracking/verbose 1");
  G4UImanager::GetUIpointer()->ApplyCommand("/geant4e/limits/stepLength 100 mm");

  //----- Choose target type:
  // 1: PlaneSurfaceTarget
  // 2: CylSurfaceTarget
  // 3: GeomVolumeTarget
  // 4: TrackLengthTarget
  G4int iTarget = 1;
  char* target = std::getenv("G4ERROR_TARGET");
  if( target ) {
    if( G4String(target) == G4String("PLANE_SURFACE") ) {
      iTarget = 1;
    }else if( G4String(target) == G4String("CYL_SURFACE") ) {
      iTarget = 2;
    }else if( G4String(target) == G4String("VOLUME") ) {
      iTarget = 3;
    }else if( G4String(target) == G4String("TRKLEN") ) {
      iTarget = 4;
    }else {
      G4Exception("exG4eReco","Fatal error in Argument",
        FatalErrorInArgument,
        G4String("Variable G4ERROR_TARGET = " + G4String(target) + 
                 "   It must be:  PLANE_SURFACE, CYL_SURFACE, VOLUME, TRKLEN").c_str());
    }
  } else {
    G4Exception("exG4eReco","Fatal error in Argument",
      JustWarning,"Variable G4ERROR_TARGET not defined, taking it = PLANE_SURFACE");
  } 

  theTarget = BuildTarget( iTarget );
  g4edata->SetTarget( theTarget );

  theG4ErrorMode = G4ErrorMode_PropBackwards;
  char* mode = std::getenv("G4ERROR_MODE");
  if( mode ) {
    if( G4String(mode) == G4String("FORWARDS") ) {
      theG4ErrorMode = G4ErrorMode_PropForwards;
    } else if( G4String(mode) == G4String("BACKWARDS") ) {
      theG4ErrorMode = G4ErrorMode_PropBackwards;
    } else {
      G4Exception("exG4eReco","Fatal error in Argument",
        FatalErrorInArgument,
        G4String("Variable G4ERROR_MODE = " + G4String(mode) + 
                 "   It must be:  FORWARDS or BACKWARDS").c_str());
    }
  } else {
    G4Exception("exG4eReco","Fatal error in Argument",
      JustWarning,"Variable G4ERROR_MODE not defined, taking it = BACKWARDS");
  } 

}


void ProcessEvent( G4int iProp, size_t )
{
  
// Set the starting trajectory.
  G4ThreeVector xv3( 0, 0, 0 );
  G4ThreeVector pv3( 20.0*GeV, 0.0, 0.0 );
  G4ErrorTrajErr error( 5, 0 );
  G4ErrorFreeTrajState* g4ErrorTrajState 
    = new G4ErrorFreeTrajState( "mu-", xv3, pv3, error );

  G4ErrorPropagatorManager* g4emgr 
    = G4ErrorPropagatorManager::GetErrorPropagatorManager();

  //int ierr = 0;

  G4Point3D surfPos(224.*cm,0.,0.);
  G4Normal3D surfNorm(1.,0.,0.);
  //-  G4ErrorTarget* theG4ErrorTarget 
  //     = new G4ErrorPlaneSurfaceTarget(surfNorm, surfPos );

  if( iProp == 0){
    // Propagate until G4ErrorTarget is found all in one go
     //ierr = 
     g4emgr->Propagate( g4ErrorTrajState, theTarget, theG4ErrorMode );
  } else if( iProp == 1){

    // Propagate until G4ErrorTarget is reached step by step
  
    g4emgr->InitTrackPropagation();

    //    G4Track* aTrack 
    //      = G4EventManager::GetEventManager()->GetTrackingManager()->GetTrack();
    bool moreEvt = TRUE;
    while( moreEvt ){
      
      //ierr = 
      g4emgr->PropagateOneStep( g4ErrorTrajState, theG4ErrorMode );
      
      //---- Check if target is reached
      if( g4emgr->GetPropagator()->CheckIfLastStep( g4ErrorTrajState->GetG4Track() )) {
        g4emgr->GetPropagator()
          ->InvokePostUserTrackingAction( g4ErrorTrajState->GetG4Track() );  
        moreEvt = 0;
        G4cout << "STEP_BY_STEP propagation: Last Step " << G4endl;
      }
    }
  }

  G4cout << " $$$ PROPAGATION ENDED " << G4endl;
  // extract current state info
  G4Point3D posEnd = g4ErrorTrajState->GetPosition();
  G4Normal3D momEnd = g4ErrorTrajState->GetMomentum();
  G4ErrorTrajErr errorEnd = g4ErrorTrajState->GetError();

  G4cout << " Position: " << posEnd << G4endl
         << " Momentum: " << momEnd << G4endl
         << " Error: " << errorEnd << G4endl; 
}


//-------------------------------------------------------------
G4ErrorTarget* BuildTarget( G4int iTarget )
{

  G4ErrorTarget* target = 0;
  if( iTarget == 1 ) {
    G4Point3D surfPos(221.*cm,0.,0.);
    G4Normal3D surfNorm(1.,0.,0.);
    target = new G4ErrorPlaneSurfaceTarget(surfNorm, surfPos );
  }else if( iTarget == 2 ) {
    G4double radius = 222*cm;
    target = new G4ErrorCylSurfaceTarget(radius);
  }else if( iTarget == 3 ) {
    target = new G4ErrorGeomVolumeTarget("MUON");
  }else if( iTarget == 4 ) {
    target = new G4ErrorTrackLengthTarget(223.*cm);
  }else {
    G4Exception("exG4eReco","Fatal error in Argument",
                FatalErrorInArgument,"Target type has to be between 1 and 4");
  }
  return target;
}


//-------------------------------------------------------------
void Finalize()
{
  G4ErrorPropagatorManager::GetErrorPropagatorManager()->CloseGeometry();

}
