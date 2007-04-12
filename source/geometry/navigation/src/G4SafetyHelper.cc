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
// $Id: G4SafetyHelper.cc,v 1.8 2007-04-12 11:51:48 vnivanch Exp $
// GEANT4 tag $ Name:  $


#include "G4SafetyHelper.hh"
#include "G4PathFinder.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"

#include "globals.hh"

//G4bool G4SafetyHelper::fUseParallelGeometries= false;  
                                            // By default, one geometry only

G4SafetyHelper::G4SafetyHelper()
{
  first = true;
  factor = 0.2;   
  fUseParallelGeometries = false;
}

void G4SafetyHelper::InitialiseNavigator()
{
  fpPathFinder= G4PathFinder::GetInstance();
 
  G4TransportationManager* pTransportMgr= 
    G4TransportationManager::GetTransportationManager();

  fpMassNavigator = pTransportMgr->GetNavigatorForTracking(); 

  // Check
  G4VPhysicalVolume* worldPV= fpMassNavigator->GetWorldVolume(); 
  if( worldPV == 0 ) { 
    G4Exception("G4SafetyHelper::ComputeMassStep",
		"InvalidNavigatorWorld", 
		FatalException, 
		"Found that existing mass Navigator has null world"); 
  }

  fMassNavigatorId = pTransportMgr->ActivateNavigator( fpMassNavigator ); 
}

void G4SafetyHelper::InitialiseHelper()
{
  lastSafetyPosition = G4ThreeVector(0.0,0.0,0.0);
  lastSafety         = 0.0;
  if(first) InitialiseNavigator();    
  first = false;
}

G4SafetyHelper::~G4SafetyHelper()
{}

G4double   
G4SafetyHelper::CheckNextStep(const G4ThreeVector &position, 
			      const G4ThreeVector &direction,
			      const G4double &currentMaxStep,
			      G4double &newSafety )
{
  // G4cout << "pos= " << position << "  dir= " << direction 
  //       << "  safety= " << newSafety<< " minStep= " << currentMaxStep << G4endl;
  
  // Distance in the Mass geometry
  G4double linstep = fpMassNavigator->CheckNextStep( position,
						     direction,
						     currentMaxStep,
						     newSafety);
  lastSafetyPosition = position;
  lastSafety         = newSafety;
  //    G4cout << "linStep= " << linstep << "  safety= " 
  //   << newSafety << " minStep= " << currentMaxStep << G4endl;

  // TO-DO: Can replace this with a call to PathFinder 
  //        giving id of Mass Geometry --> this avoid doing the work twice
  return linstep;
}

G4double G4SafetyHelper::ComputeSafety( const G4ThreeVector& position ) 
{
  G4double newSafety;

  // return last value if position is not significantly changed 
  G4double moveLen=(position-lastSafetyPosition).mag();
  if(moveLen >= factor*lastSafety) {
    
    lastSafetyPosition = position;
 
    if( !fUseParallelGeometries) {
      // Safety for mass geometry
      lastSafety = fpMassNavigator->ComputeSafety(position); 
    } else {
      // Safety for all geometries
      lastSafety= fpPathFinder->ComputeSafety( position ); 
    } 
    newSafety= lastSafety;
  }else{
    newSafety= lastSafety-moveLen;
  } 
  // G4cout << "G4SafetyHelper::ComputeSafety: newSafety= " << newSafety << G4endl; 
  return newSafety;
}

void  G4SafetyHelper::ReLocateWithinVolume( const G4ThreeVector &newPosition )
{
  if( !fUseParallelGeometries) {
    fpMassNavigator->LocateGlobalPointWithinVolume( newPosition ); 
  }else{
    fpPathFinder->ReLocate( newPosition ); 
  }
}

void  G4SafetyHelper::Locate( const G4ThreeVector& newPosition, 
                              const G4ThreeVector& newDirection)
{
  if( !fUseParallelGeometries) {
    fpMassNavigator->LocateGlobalPointAndSetup( newPosition, &newDirection, true, false); 
  }else{
    fpPathFinder->Locate( newPosition, newDirection ); 
  }
}
