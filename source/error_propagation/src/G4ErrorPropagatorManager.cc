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
// $Id: G4ErrorPropagatorManager.cc 107054 2017-11-01 14:35:18Z gcosmo $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file 
// ------------------------------------------------------------
//

#include "G4ErrorPropagatorManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4Mag_EqRhs.hh"
#include "G4MagIntegratorDriver.hh"

#include "G4ClassicalRK4.hh"
#include "G4ExactHelixStepper.hh"
#include "G4HelixExplicitEuler.hh"

#include "G4EventManager.hh"
#include "G4ErrorRunManagerHelper.hh"
#include "G4ErrorPropagator.hh"
#include "G4ErrorMag_UsualEqRhs.hh"

#include "G4VParticleChange.hh"
#include "G4ParticleChangeForMSC.hh"
#include "G4ParticleChange.hh"
#include "G4Track.hh"
#include "G4TransportationManager.hh"
#include "G4ErrorPropagationNavigator.hh"
#include "G4GeometryManager.hh"
#include "G4StateManager.hh"
#include "G4ChordFinder.hh"
#include "G4EquationOfMotion.hh"
#include "G4FieldManager.hh"
#include "G4PropagatorInField.hh"
#include "G4RunManager.hh"
#include "G4VParticleChange.hh"

G4ThreadLocal G4ErrorPropagatorManager*
G4ErrorPropagatorManager::theG4ErrorPropagatorManager = 0;

//-----------------------------------------------------------------------
G4ErrorPropagatorManager* G4ErrorPropagatorManager::GetErrorPropagatorManager()
{
  if( theG4ErrorPropagatorManager == NULL ) {
    theG4ErrorPropagatorManager = new G4ErrorPropagatorManager;
  }

  return theG4ErrorPropagatorManager;
}


//-----------------------------------------------------------------------
G4ErrorPropagatorManager::G4ErrorPropagatorManager()
{
  //----- Initialize a few things
  //o  theG4ErrorPropagatorManager = this;

  char* g4emverb = getenv("G4EVERBOSE");
  if( !g4emverb ) {
    G4ErrorPropagatorData::GetErrorPropagatorData()->SetVerbose( 0 );
  } else {
    G4ErrorPropagatorData::GetErrorPropagatorData()->SetVerbose( atoi( g4emverb ) );
  }

  thePropagator = 0;

  theEquationOfMotion = 0;

  StartG4ErrorRunManagerHelper(); 
  
  G4ErrorPropagatorData::GetErrorPropagatorData()->SetState( G4ErrorState_PreInit );

  theG4ErrorPropagationNavigator = 0;

  StartNavigator(); //navigator has to be initialized at the beggining !?!?!


}


//-----------------------------------------------------------------------
G4ErrorPropagatorManager::~G4ErrorPropagatorManager()
{
  delete theEquationOfMotion;
  delete theG4ErrorPropagationNavigator;
  delete thePropagator;
  delete theG4ErrorRunManagerHelper;
  delete theG4ErrorPropagatorManager;
}


//-----------------------------------------------------------------------
void G4ErrorPropagatorManager::StartG4ErrorRunManagerHelper()
{
  //----- Initialize G4ErrorRunManagerHelper
  theG4ErrorRunManagerHelper = G4ErrorRunManagerHelper::GetRunManagerKernel();

  if( theG4ErrorRunManagerHelper == 0 ) {
    theG4ErrorRunManagerHelper = new G4ErrorRunManagerHelper();
  }

  //----- User Initialization classes 
  //--- GEANT4e PhysicsList
  if( G4ErrorPropagatorData::verbose() >= 4 ) G4cout << " G4ErrorPropagatorManager::StartG4eRunManager() done " << theG4ErrorRunManagerHelper << G4endl;
  //-  theG4eRunManager->SetUserInitialization(new G4ErrorPhysicsList);

}


//-----------------------------------------------------------------------
void G4ErrorPropagatorManager::StartNavigator()
{
  if( theG4ErrorPropagationNavigator == 0 ) {
    G4TransportationManager*transportationManager = G4TransportationManager::GetTransportationManager();
    
    G4Navigator* g4navi = transportationManager->GetNavigatorForTracking();
   
    G4VPhysicalVolume* world = g4navi->GetWorldVolume();
    G4int verb = g4navi->GetVerboseLevel();
    delete g4navi;

    theG4ErrorPropagationNavigator = new G4ErrorPropagationNavigator();

    if( world != 0 ) {
      theG4ErrorPropagationNavigator->SetWorldVolume( world );
    }
    theG4ErrorPropagationNavigator->SetVerboseLevel( verb );   
    
    transportationManager->SetNavigatorForTracking(theG4ErrorPropagationNavigator);
    transportationManager->GetPropagatorInField()->GetIntersectionLocator()
                         ->SetNavigatorFor(theG4ErrorPropagationNavigator);
    G4EventManager::GetEventManager()->GetTrackingManager()->GetSteppingManager()
                         ->SetNavigator(theG4ErrorPropagationNavigator);
    //  G4ThreeVector center(0,0,0);
    //  theG4ErrorPropagationNavigator->LocateGlobalPointAndSetup(center,0,false);
    
  }

  if( G4ErrorPropagatorData::verbose() >= 2 ) G4cout << " theState at StartNavigator " << PrintG4ErrorState() << G4endl;

}
  

//-----------------------------------------------------------------------
void G4ErrorPropagatorManager::InitGeant4e()
{
  if( G4ErrorPropagatorData::verbose() >= 1 ) G4cout << "InitGeant4e GEANT4e State= " << PrintG4ErrorState() << " GEANT4 State= " << PrintG4State() << G4endl;
  G4ApplicationState currentState = G4StateManager::GetStateManager()->GetCurrentState();
  //----- Initialize run
  //  if( G4StateManager::GetStateManager()->GetCurrentState() == G4State_PreInit) {
  
  if( G4ErrorPropagatorData::GetErrorPropagatorData()->GetState() == G4ErrorState_PreInit ) {
    if ( currentState == G4State_PreInit || currentState == G4State_Idle) {
      //    G4eRunManager::GetRunManager()->Initialize();
      theG4ErrorRunManagerHelper->InitializeGeometry();
      theG4ErrorRunManagerHelper->InitializePhysics();
    }
    
    InitFieldForBackwards();
    
    //-    G4StateManager::GetStateManager()->SetNewState(G4State_Idle);
    
    if( G4ErrorPropagatorData::verbose() >= 4 )   G4cout << " bef  theG4ErrorPropagatorManager->RunInitialization() " <<  G4StateManager::GetStateManager()->GetCurrentState() << G4endl;
    theG4ErrorRunManagerHelper->RunInitialization();
    if( G4ErrorPropagatorData::verbose() >= 4 ) G4cout << " aft  theG4ErrorPropagatorManager->RunInitialization() " <<  G4StateManager::GetStateManager()->GetCurrentState() << G4endl;
    
    if( !thePropagator ) thePropagator = new G4ErrorPropagator();  // currently the only propagator possible
    
    InitTrackPropagation();
  } else {
    std::ostringstream message;
    message << "Illegal GEANT4e State= " << PrintG4ErrorState();
    G4Exception("G4ErrorPropagatorManager::InitGeant4e()",
                "IllegalState", JustWarning, message);
  }
  
  //----- Set the tracking geometry for this propagation
  //t  SetTrackingGeometry();
  //----- Set the physics list for this propagation
  //t  SetPhysicsList();
  //----- Set the field propagation parameters for this propagation
  //t  SetFieldPropagationParameters();
  G4ErrorPropagatorData::GetErrorPropagatorData()->SetState( G4ErrorState_Init );

  if( G4ErrorPropagatorData::verbose() >= 2 ) G4cout << "End InitGeant4e GEANT4e State= " << PrintG4ErrorState() << " GEANT4 State= " << PrintG4State() << G4endl;


}


//-----------------------------------------------------------------------
void G4ErrorPropagatorManager::InitTrackPropagation()
{
  thePropagator->SetStepN( 0 );

  G4ErrorPropagatorData::GetErrorPropagatorData()->SetState( G4ErrorState_Propagating );

}


//-----------------------------------------------------------------------
G4bool G4ErrorPropagatorManager::InitFieldForBackwards()
{

  if( G4ErrorPropagatorData::verbose() >= 4 ) G4cout << " G4ErrorPropagatorManager::InitFieldForBackwards() " << G4endl;
  //----- Gets the current equation of motion
  G4FieldManager* fieldMgr= G4TransportationManager::GetTransportationManager()->GetFieldManager();
  //  G4cout << " fieldMgr " << fieldMgr << G4endl;
  if( !fieldMgr ) return 0;

  //  G4Field* myfield = fieldMgr->GetDetectorField();
  G4ChordFinder* cf = fieldMgr ->GetChordFinder();
  if( !cf ) return 0;
  auto driver = cf->GetIntegrationDriver();
  if( !driver ) return 0;
  auto equation = driver->GetEquationOfMotion();

  //----- Replaces the equation by a G4ErrorMag_UsualEqRhs to handle backwards tracking
  if ( !dynamic_cast<G4ErrorMag_UsualEqRhs*>(equation) ) {

    G4MagneticField* myfield = (G4MagneticField*)fieldMgr->GetDetectorField();
    
    //    G4Mag_UsualEqRhs* fEquation_usual = dynamic_cast<G4Mag_UsualEqRhs*>(equation);
    if( theEquationOfMotion == 0 ) theEquationOfMotion = new G4ErrorMag_UsualEqRhs(myfield);
 
    //---- Pass the equation of motion to the G4MagIntegratorStepper
    driver->SetEquationOfMotion( theEquationOfMotion );

    //--- change stepper for speed tests
   G4MagIntegratorStepper* g4eStepper = new G4ClassicalRK4(theEquationOfMotion);
   // G4MagIntegratorStepper* g4eStepper = new G4ExactHelixStepper(theEquationOfMotion);
    
    //---- 
    G4MagneticField* field = static_cast<G4MagneticField*>(const_cast<G4Field*>(fieldMgr->GetDetectorField()));
    G4ChordFinder* pChordFinder = new G4ChordFinder(field, 1.0e-2*mm, g4eStepper);

    fieldMgr->SetChordFinder(pChordFinder);

  }

  return 1;
}


//-----------------------------------------------------------------------
G4int G4ErrorPropagatorManager::Propagate( G4ErrorTrajState* currentTS, const G4ErrorTarget* target, G4ErrorMode mode )
{
  G4ErrorPropagatorData::GetErrorPropagatorData()->SetMode( mode );
  if( !thePropagator ) thePropagator = new G4ErrorPropagator();  // currently the only propagator possible

  SetSteppingManagerVerboseLevel();
  InitTrackPropagation();

  G4int ierr = thePropagator->Propagate( currentTS, target, mode );

  EventTermination();

  return ierr;
}


//-----------------------------------------------------------------------
G4int G4ErrorPropagatorManager::PropagateOneStep( G4ErrorTrajState* currentTS, G4ErrorMode mode )
{
  G4ErrorPropagatorData::GetErrorPropagatorData()->SetMode( mode );

  if( !thePropagator ) thePropagator = new G4ErrorPropagator();  // currently the only propagator possible

  SetSteppingManagerVerboseLevel();

  return thePropagator->PropagateOneStep( currentTS );
}


//-----------------------------------------------------------------------
G4bool G4ErrorPropagatorManager::CloseGeometry()
{
  G4GeometryManager* geomManager = G4GeometryManager::GetInstance();
  geomManager->OpenGeometry();
  if(  G4StateManager::GetStateManager()->GetCurrentState() != G4State_GeomClosed) {
    G4StateManager::GetStateManager()->SetNewState(G4State_Quit);
  }

  return TRUE;
}


//---------------------------------------------------------------------------
void G4ErrorPropagatorManager::SetUserInitialization(G4VUserDetectorConstruction* userInit)
{
  theG4ErrorRunManagerHelper->SetUserInitialization( userInit); 
}


//---------------------------------------------------------------------------
void G4ErrorPropagatorManager::SetUserInitialization(G4VPhysicalVolume* userInit)
{ 
  theG4ErrorRunManagerHelper->SetUserInitialization( userInit); 
}
 

//---------------------------------------------------------------------------
void G4ErrorPropagatorManager::SetUserInitialization(G4VUserPhysicsList* userInit)
{ 
  theG4ErrorRunManagerHelper->SetUserInitialization( userInit); 
}


//---------------------------------------------------------------------------
void G4ErrorPropagatorManager::SetUserAction(G4UserTrackingAction* userAction)
{
  G4EventManager::GetEventManager()->SetUserAction( userAction ); 
}


//---------------------------------------------------------------------------
void G4ErrorPropagatorManager::SetUserAction(G4UserSteppingAction* userAction)
{
  G4EventManager::GetEventManager()->SetUserAction( userAction ); 
}


//---------------------------------------------------------------------------
void G4ErrorPropagatorManager::SetSteppingManagerVerboseLevel()
{
  G4TrackingManager* trkmgr = G4EventManager::GetEventManager()->GetTrackingManager();
  trkmgr->GetSteppingManager()->SetVerboseLevel( trkmgr->GetVerboseLevel() );
}


//---------------------------------------------------------------------------
void G4ErrorPropagatorManager::EventTermination()
{
  G4ErrorPropagatorData::GetErrorPropagatorData()->SetState( G4ErrorState_Init );
}


//---------------------------------------------------------------------------
void G4ErrorPropagatorManager::RunTermination()
{
G4ErrorPropagatorData::GetErrorPropagatorData()->SetState( G4ErrorState_PreInit );
  theG4ErrorRunManagerHelper->RunTermination(); 
}


//---------------------------------------------------------------------------
G4String G4ErrorPropagatorManager::PrintG4ErrorState() 
{
  return PrintG4ErrorState( G4ErrorPropagatorData::GetErrorPropagatorData()->GetState() );
}


//---------------------------------------------------------------------------
G4String G4ErrorPropagatorManager::PrintG4ErrorState( G4ErrorState state ) 
{
  G4String nam = "";
  switch (state){
  case G4ErrorState_PreInit: 
    nam = "G4ErrorState_PreInit"; 
    break;
  case G4ErrorState_Init: 
    nam = "G4ErrorState_Init"; 
    break;
  case G4ErrorState_Propagating:
    nam = "G4ErrorState_Propagating";
    break;
  case G4ErrorState_TargetCloserThanBoundary:
    nam = "G4ErrorState_TargetCloserThanBoundary";
    break;
  case G4ErrorState_StoppedAtTarget:
    nam = "G4ErrorState_StoppedAtTarget";
    break;
  }

  return nam;
}


//---------------------------------------------------------------------------
G4String G4ErrorPropagatorManager::PrintG4State()
{
  return PrintG4State(G4StateManager::GetStateManager()->GetCurrentState());
}


//---------------------------------------------------------------------------
G4String G4ErrorPropagatorManager::PrintG4State( G4ApplicationState state )
{
  G4String nam = "";
  switch ( state ){
  case G4State_PreInit:
    nam = "G4State_PreInit";
    break;
  case G4State_Init:
    nam = "G4State_Init";
    break;
  case G4State_Idle:
    nam = "G4State_Idle";
    break;
  case G4State_GeomClosed:
    nam = "G4State_GeomClosed";
    break;
  case G4State_EventProc:
    nam = "G4State_EventProc"; 
    break;
  case G4State_Quit:
    nam = "G4State_Quit";
    break;
  case G4State_Abort:
    nam = "G4State_Abort";
    break;
  }
  
  return nam;

}
