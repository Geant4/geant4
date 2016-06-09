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
// $Id: G4SafetyHelper.cc,v 1.7 2006/11/14 10:22:12 japost Exp $
// GEANT4 tag $ Name:  $


#include "G4SafetyHelper.hh"
#include "G4PathFinder.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"

// #include "G4Exception.hh"
#include "globals.hh"

G4bool G4SafetyHelper::fUseParallelGeometries= false;  
                                            // By default, one geometry only

G4SafetyHelper::G4SafetyHelper()
{
   fpPathFinder= G4PathFinder::GetInstance();
   
   InitialiseNavigator();    
}

void G4SafetyHelper::InitialiseNavigator()
{
   // G4Navigator* 
   static G4TransportationManager* pTransportMgr= 
     G4TransportationManager::GetTransportationManager();

   fpMassNavigator = pTransportMgr->GetNavigatorForTracking(); 

   fMassNavigatorId = pTransportMgr->ActivateNavigator( fpMassNavigator ); 
} 

G4SafetyHelper::~G4SafetyHelper()
{}

G4double   
G4SafetyHelper::ComputeMassStep( const G4ThreeVector &position, 
				 const G4ThreeVector &direction,
				 G4double  &newSafety )
{
   // Step for mass geometry
    G4double linearStep;
    const G4double proposedStep = DBL_MAX;

    // Check
    G4VPhysicalVolume* worldPV= fpMassNavigator->GetWorldVolume(); 
    if( worldPV == 0 ) { 
       G4Exception("G4SafetyHelper::ComputeMassStep",
		   "InvalidNavigatorWorld", 
		   FatalException, 
		   "Found that existing mass Navigator has null world"); 
    }

    fpMassNavigator->LocateGlobalPointWithinVolume(position);
        // Potentially dangerous to relocate the point.  
        //   Safe in PostStepDoIt -- and possibly in AlongStepGPIL

    G4cout << "G4SafetyHelper::ComputeMassStep " 
	   << " trial step size = " << proposedStep << " ." << G4endl; 

    // Distance in the Mass geometry
    linearStep = fpMassNavigator->ComputeStep( position,
					 direction,
					 proposedStep,
					 newSafety);

    fpMassNavigator->LocateGlobalPointWithinVolume(position);
    G4cout << "G4UrbanMscModel relocates mass Navigator back to " 
	   << position  << G4endl;

    // TO-DO: Can replace this with a call to PathFinder 
    //        giving id of Mass Geometry --> this avoid doing the work twice

    return linearStep; 
}

G4double   G4SafetyHelper::ComputeSafety( const G4ThreeVector& position ) 
{
   // Safety for all geometries
   G4double newsafety= 0.0; 

   if( !fUseParallelGeometries) {
      // Old code: safety for mass geometry
      fpMassNavigator->LocateGlobalPointWithinVolume(position);
      newsafety = fpMassNavigator->ComputeSafety(position); 
   }else{

      // fpPathFinder->ReLocate( position );   // Safe in PostStepDoIt only ??
      newsafety= fpPathFinder->ComputeSafety( position ); 

#ifdef CHECK_WITH_ONE_GEOM 
      // Check against mass safety
      fpMassNavigator->LocateGlobalPointWithinVolume(position);
      G4double mass_safety = fpMassNavigator->ComputeSafety(position); 

      // For initial tests check assume that mass is only geometry
      if( (mass_safety - newsafety) > 1e-4 * newsafety ){
         G4cerr << " ERROR in G4SafetyHelper " << G4endl
	     << "   Safety from PathFinder is " << newsafety << " "
	     << "    not equal to " <<  mass_safety << "  " << G4endl;
	 G4Exception("G4SafetyHelper::ComputeSafety", "SafetyError",
                   FatalException, 
		  "Incompatible safeties between navigator and pathfinder" );
	 exit(1); 
      }
#endif
   }
   return newsafety;
}


void  G4SafetyHelper::ReLocateWithinVolume( const G4ThreeVector &newPosition )
{

#ifdef G4VERBOSE_HELPER
   G4int oldPrec= G4cout.precision( 10 ); 
   G4cout << "  G4SafetyHelper::ReLocateWithinVolume " 
	  << " calling PathFinder->ReLocating at position " << newPosition << G4endl;
   G4cout.precision( oldPrec ); 
#endif

   if( !fUseParallelGeometries) {
      fpMassNavigator->LocateGlobalPointWithinVolume( newPosition ); 
   }else{
      fpPathFinder->ReLocate( newPosition ); 
   }
}
