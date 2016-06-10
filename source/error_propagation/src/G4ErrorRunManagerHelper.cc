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
// $Id: G4ErrorRunManagerHelper.cc 78318 2013-12-11 15:02:40Z gcosmo $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file 
// ------------------------------------------------------------
//

#include "G4Navigator.hh"

#include "G4Timer.hh"

#include "G4ErrorRunManagerHelper.hh"

#include "G4RunManagerKernel.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4ErrorPhysicsList.hh"
#include "G4TransportationManager.hh"
#include "G4RunManager.hh"

//-----------------------------------------------------------------------

G4ThreadLocal G4ErrorRunManagerHelper*
G4ErrorRunManagerHelper::fRunManagerKernel = 0;

//-----------------------------------------------------------------------
G4ErrorRunManagerHelper* G4ErrorRunManagerHelper::GetRunManagerKernel()
{ return fRunManagerKernel; }


//-----------------------------------------------------------------------
G4ErrorRunManagerHelper::G4ErrorRunManagerHelper()
{
  if(fRunManagerKernel) {
    G4Exception("G4ErrorRunManagerHelper::G4ErrorRunManagerHelper()",
                "InvalidSetup", FatalException,
                "G4eRunManageKernel constructed twice.");
  }
  fRunManagerKernel = this;

  //----- Look if somebody has created a G4RunManagerKernel
  theG4RunManagerKernel = G4RunManagerKernel::GetRunManagerKernel();
  if( theG4RunManagerKernel == 0 ) {
    //--- if not create it
    theG4RunManagerKernel = new G4RunManagerKernel();
    G4cout << " creating G4RunManagerKernel " <<  theG4RunManagerKernel << G4endl;
  }
    
  theG4RunManagerKernel->SetVerboseLevel(2);
  theUserPhysicsList = 0;
  theUserWorld = 0;

}


//-----------------------------------------------------------------------
G4ErrorRunManagerHelper::~G4ErrorRunManagerHelper()
{
}


//-----------------------------------------------------------------------
void G4ErrorRunManagerHelper::SetUserInitialization(G4VUserDetectorConstruction* userInit)
{ 
  theUserWorld = userInit->Construct();
}


//-----------------------------------------------------------------------
void G4ErrorRunManagerHelper::SetUserInitialization(G4VPhysicalVolume* userInit)
{
  theUserWorld = userInit;
}


//-----------------------------------------------------------------------
void G4ErrorRunManagerHelper::SetUserInitialization(G4VUserPhysicsList* userInit)
{
  theUserPhysicsList = userInit;
}


//-----------------------------------------------------------------------
void G4ErrorRunManagerHelper::InitializeGeometry()
{
  //check if user world has been directly called or someone initialized the world volume already 
  //----- First option: geometry has been defined to GEANT4e 
  if( theUserWorld != 0 ) {
    theG4RunManagerKernel->DefineWorldVolume( theUserWorld );
    
    //----- Second option: geometry has been defined to GEANT4, do nothing GEANT4 should take care 
  } else {
    //--- Check that indeed geometry has been defined to GEANT4
    if ( G4TransportationManager::GetTransportationManager()
         ->GetNavigatorForTracking()->GetWorldVolume() == 0 ) {
      G4Exception("G4ErrorRunManagerHelper::InitializeGeometry()",
                  "InvalisSetup", FatalException,
                  "No world defined in your geometry!" );
    }
    
  }
}


//-----------------------------------------------------------------------
void G4ErrorRunManagerHelper::InitializePhysics()
{

  G4cout << "  G4ErrorRunManagerHelper::InitializePhysics "  << G4endl;

  //----- First option: physics list has been defined to GEANT4e 
  if( theUserPhysicsList != 0 ) {
    theG4RunManagerKernel->SetPhysics(theUserPhysicsList);
    theG4RunManagerKernel->InitializePhysics();
  }else {
  //----- Second option: physics list has been defined to GEANT4, do nothing GEANT4 should take care 
    if( G4RunManager::GetRunManager() != 0 && G4RunManager::GetRunManager()->GetUserPhysicsList() != 0 ){ 
      //--- Physics should be G4ErrorPhysicsList, else send a warning
      if( static_cast<const G4ErrorPhysicsList*>(G4RunManager::GetRunManager()->GetUserPhysicsList()) == 0 ) {
        std::ostringstream message;
        message << "Physics list is not G4ErrorPhysicsList. Are you sure?";
        G4Exception("G4ErrorRunManagerHelper::InitializePhysics()",
                    "GEANT4e-Notification", JustWarning, message);
      }
    } else {
      //----- Third option: no physics list has been defined, define a G4ErrorPhysicsList
      theG4RunManagerKernel->SetPhysics(new G4ErrorPhysicsList);
      //    theG4RunManagerKernel->SetPhysics(new ExN02PhysicsList);
      theG4RunManagerKernel->InitializePhysics();
    }
  }
 
}


//-----------------------------------------------------------------------
void G4ErrorRunManagerHelper::RunInitialization()
{
  theG4RunManagerKernel->RunInitialization();
}


//-----------------------------------------------------------------------
void G4ErrorRunManagerHelper::SetUserAction(G4UserTrackingAction* userAction)
{

  G4EventManager::GetEventManager()->SetUserAction( userAction );
}


//-----------------------------------------------------------------------
void G4ErrorRunManagerHelper::SetUserAction(G4UserSteppingAction* userAction)
{
  G4EventManager::GetEventManager()->SetUserAction( userAction );
}


//-----------------------------------------------------------------------
void G4ErrorRunManagerHelper::RunTermination()
{
  theG4RunManagerKernel->RunTermination();
}

