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
// $Id: G4MaterialScanner.cc 66892 2013-01-17 10:57:59Z gunter $
//
//
//

#include "G4MaterialScanner.hh"

#include "G4SystemOfUnits.hh"
#include "G4EventManager.hh"
#include "G4MatScanMessenger.hh"
#include "G4RayShooter.hh"
#include "G4MSSteppingAction.hh"
#include "G4GeometryManager.hh"
#include "G4StateManager.hh"
#include "G4Event.hh"
#include "G4TransportationManager.hh"
#include "G4RunManagerKernel.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4SDManager.hh"


G4MaterialScanner::G4MaterialScanner()
{
  theRayShooter = new G4RayShooter();
  theMessenger = new G4MatScanMessenger(this);
  theEventManager = G4EventManager::GetEventManager();

  theUserEventAction = 0;
  theUserStackingAction = 0;
  theUserTrackingAction = 0;
  theUserSteppingAction = 0;

  theMatScannerEventAction = 0;
  theMatScannerStackingAction = 0;
  theMatScannerTrackingAction = 0;
  theMatScannerSteppingAction = 0;

  eyePosition = G4ThreeVector(0.,0.,0.);
  nTheta = 91;
  thetaMin = 0.*deg;
  thetaSpan = 90.*deg;
  nPhi = 37;
  phiMin = 0.*deg;
  phiSpan = 360.*deg;

  regionSensitive = false;
  regionName = "notDefined";
  theRegion = 0;
}

G4MaterialScanner::~G4MaterialScanner()
{
  delete theRayShooter;
  delete theMatScannerSteppingAction;
  delete theMessenger;
}

void G4MaterialScanner::Scan()
{
  G4StateManager* theStateMan = G4StateManager::GetStateManager();
  G4ApplicationState currentState = theStateMan->GetCurrentState();
  if(currentState!=G4State_Idle)
  {
    G4cerr << "Illegal application state - Scan() ignored." << G4endl;
    return;
  }

  if(!theMatScannerSteppingAction)
  { theMatScannerSteppingAction = new G4MSSteppingAction(); }
  StoreUserActions();
  DoScan();
  RestoreUserActions();
}

void G4MaterialScanner::StoreUserActions()
{ 
  theUserEventAction = theEventManager->GetUserEventAction();
  theUserStackingAction = theEventManager->GetUserStackingAction();
  theUserTrackingAction = theEventManager->GetUserTrackingAction();
  theUserSteppingAction = theEventManager->GetUserSteppingAction();

  theEventManager->SetUserAction(theMatScannerEventAction);
  theEventManager->SetUserAction(theMatScannerStackingAction);
  theEventManager->SetUserAction(theMatScannerTrackingAction);
  theEventManager->SetUserAction(theMatScannerSteppingAction);

  G4SDManager* theSDMan = G4SDManager::GetSDMpointerIfExist();
  if(theSDMan)
  { theSDMan->Activate("/",false); }

  G4GeometryManager* theGeomMan = G4GeometryManager::GetInstance();
  theGeomMan->OpenGeometry();
  theGeomMan->CloseGeometry(true);
}

void G4MaterialScanner::RestoreUserActions()
{
  theEventManager->SetUserAction(theUserEventAction);
  theEventManager->SetUserAction(theUserStackingAction);
  theEventManager->SetUserAction(theUserTrackingAction);
  theEventManager->SetUserAction(theUserSteppingAction);

  G4SDManager* theSDMan = G4SDManager::GetSDMpointerIfExist();
  if(theSDMan)
  { theSDMan->Activate("/",true); }
}

void G4MaterialScanner::DoScan()
{
// Confirm material table is updated
  G4RunManagerKernel::GetRunManagerKernel()->UpdateRegion();

/////  // Make sure Geantino has been initialized
/////  G4ProcessVector* pVector
/////    = G4Geantino::GeantinoDefinition()->GetProcessManager()->GetProcessList();
/////  for (G4int j=0; j < pVector->size(); ++j) {
/////      (*pVector)[j]->BuildPhysicsTable(*(G4Geantino::GeantinoDefinition()));
/////  }

// Close geometry and set the application state
  G4GeometryManager* geomManager = G4GeometryManager::GetInstance();
  geomManager->OpenGeometry();
  geomManager->CloseGeometry(1,0);
  
  G4ThreeVector center(0,0,0);
  G4Navigator* navigator =
      G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  navigator->LocateGlobalPointAndSetup(center,0,false);

  G4StateManager* theStateMan = G4StateManager::GetStateManager();
  theStateMan->SetNewState(G4State_GeomClosed); 

// Event loop
  G4int iEvent = 0;
  for(G4int iTheta=0;iTheta<nTheta;iTheta++)
  {
   G4double theta = thetaMin;
   if(iTheta>0) theta += G4double(iTheta)*thetaSpan/G4double(nTheta-1);
   G4double aveLength = 0.;
   G4double aveX0 = 0.;
   G4double aveLambda = 0.;
   G4cout << G4endl;
   G4cout << "         Theta(deg)    Phi(deg)  Length(mm)          x0     lambda0" << G4endl;
   G4cout << G4endl;
   for(G4int iPhi=0;iPhi<nPhi;iPhi++)
   {
    G4Event* anEvent = new G4Event(iEvent++);
    G4double phi = phiMin;
    if(iPhi>0) phi += G4double(iPhi)*phiSpan/G4double(nPhi-1);
    eyeDirection = G4ThreeVector(std::cos(theta)*std::cos(phi),
                                 std::cos(theta)*std::sin(phi),
                                 std::sin(theta));
    theRayShooter->Shoot(anEvent,eyePosition,eyeDirection);
    theMatScannerSteppingAction->Initialize(regionSensitive,theRegion);
    theEventManager->ProcessOneEvent(anEvent);
    G4double length = theMatScannerSteppingAction->GetTotalStepLength();
    G4double x0 = theMatScannerSteppingAction->GetX0();
    G4double lambda = theMatScannerSteppingAction->GetLambda0();

    G4cout << "        "
           << std::setw(11) << theta/deg << " "
           << std::setw(11) << phi/deg << " "
           << std::setw(11) << length/mm << " "
           << std::setw(11) << x0 << " "
           << std::setw(11) << lambda << G4endl;
    aveLength += length/mm;
    aveX0 += x0;
    aveLambda += lambda;
   }
   if(nPhi>1)
   {
    G4cout << G4endl;
    G4cout << " ave. for theta = " << std::setw(11) << theta/deg << " : "
           << std::setw(11) << aveLength/nPhi << " "
           << std::setw(11) << aveX0/nPhi << " "
           << std::setw(11) << aveLambda/nPhi << G4endl;
   }
  }

  theStateMan->SetNewState(G4State_Idle); 
  return;
}

G4bool G4MaterialScanner::SetRegionName(const G4String& val)
{
  G4Region* aRegion = G4RegionStore::GetInstance()->GetRegion(val);
  if(aRegion)
  {
    theRegion = aRegion;
    regionName = val;
    return true;
  }
  else
  {
    G4cerr << "Region <" << val << "> not found. Command ignored." << G4endl;
    G4cerr << "Defined regions are : " << G4endl;
    for(size_t i=0;i<G4RegionStore::GetInstance()->size();i++)
    { G4cerr << " " << (*(G4RegionStore::GetInstance()))[i]->GetName(); }
    G4cerr << G4endl;
    return false;
  }
}
