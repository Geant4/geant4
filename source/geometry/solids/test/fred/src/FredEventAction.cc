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
// FredEventAction.cc
//
// Implementation of Fred's user action
//

#include "FredEventAction.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4SDManager.hh"

#include "FredHit.hh"
#include "FredMessenger.hh"

#include "FredTrackCheck.hh"

//
// Constructor
//
FredEventAction::FredEventAction( FredMessenger *ourMessenger )
{
	hitID = -1;
	hitIDMother = -1;
	messenger = ourMessenger;
}

//
// Destructor
//
FredEventAction::~FredEventAction() {;}

//
// BeginOfEventAction
//
// Do something interesting at the beginning of each event
//
void FredEventAction::BeginOfEventAction(const G4Event*)
{
	//
	// Get the hit collection id for our test hits
	//
	if (hitID < 0) {
		G4SDManager *sensitiveMan = G4SDManager::GetSDMpointer();
		hitID = sensitiveMan->GetCollectionID( "fredsStuff" );
		hitIDMother = sensitiveMan->GetCollectionID( "fredsMotherStuff" );
	}
}

//
// EndOfEventAction
//
// Do something interesting at the end of each event
//
void FredEventAction::EndOfEventAction(const G4Event*)
{

	G4UImanager *UI = G4UImanager::GetUIpointer();
	
	//
	// Draw the event automatically
	//
	G4VVisManager *visManager = G4VVisManager::GetConcreteInstance();
	
	if (visManager) {
	
		// Prepare
		UI->ApplyCommand( "/vis~/clear/view" );
		
		// Draw detector
		UI->ApplyCommand( "/vis~/draw/current" );
		
		// Draw stuff
		if (messenger->SelectedDrawing() == NORMAL) 
			DrawNormal();
		else
			DrawShadow();

		// Finish
		UI->ApplyCommand( "/vis~/show/view" );
	}
}


//
// DrawNormal
//
// Standard type GEANT4 display
//
void FredEventAction::DrawNormal()
{
	const G4Event	*evt = fpEventManager->GetConstCurrentEvent();

	// Draw trajectories, all in this case
	G4TrajectoryContainer *trajectoryContainer = evt->GetTrajectoryContainer();
	G4int numTraj = (trajectoryContainer) ? trajectoryContainer->entries() : 0;
	
	while(numTraj--) {
	 	G4VTrajectory *trajectory = (*trajectoryContainer)[numTraj];
	  	trajectory->DrawTrajectory();
	}
	
	// Draw hits
	G4HCofThisEvent *HCE = evt->GetHCofThisEvent();
	FredHitCollection *hits = (FredHitCollection *)HCE->GetHC( hitID );
	if (hits) hits->DrawAllHits();
}


//
// DrawShadow
//
// Draw tracks that miss or have errors
//
void FredEventAction::DrawShadow()
{
	const G4Event	*evt = fpEventManager->GetConstCurrentEvent();
	G4HCofThisEvent *HCE = evt->GetHCofThisEvent();
	FredHitCollection *hits = (FredHitCollection *)HCE->GetHC( hitID );

	//
	// Build a list of hit/track associations
	//
	FredTrackCheck	checkTracks;
	
	G4int numHit = (hits) ? hits->entries() : 0;
	while(numHit--) {
		FredHit *hit = (*hits)[numHit];
		if (hit->GetEnters()) 
			checkTracks.AddInHit( hit->GetTrack() );
		else
			checkTracks.AddOutHit( hit->GetTrack() );
	}
	
	//
	// Now, loop over mother hits
	//
	hits = (FredHitCollection *)HCE->GetHC( hitIDMother );

	numHit = (hits) ? hits->entries() : 0;
	while(numHit--) {
		G4int nIn, nOut;
		
		FredHit *hit = (*hits)[numHit];

	  	checkTracks.GetTrackStat( hit->GetTrack(), &nIn, &nOut );
		
		if ((nIn == 0) && (nOut == 0)) {
		
			// Pristeen track
	    		hit->DrawShadowHit();
		}
		else if (nIn != nOut) {

			// Oh oh! Track error
	  		hit->DrawErrorHit();
		}
		else {
		
			// Hit track
	  		hit->DrawMaskHit();
		}
	}
	
}
