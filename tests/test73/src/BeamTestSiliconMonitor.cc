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
#include "BeamTestSiliconMonitor.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"    
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"

BeamTestSiliconMonitor::BeamTestSiliconMonitor(const G4String& name)
   :G4VSensitiveDetector(name)
{
   collectionName.insert( "MonitorCollection" );
   fHitsCollectionID = -1;
}

BeamTestSiliconMonitor::~BeamTestSiliconMonitor() {}

void BeamTestSiliconMonitor::Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent)
{
   // HandsOn4: Creating hit collection
   // Create a new collection
   fHitsCollection =
      new BeamTestSiliconMonitorHitsCollection(SensitiveDetectorName, collectionName[0]);
 
  if ( fHitsCollectionID < 0 )
       fHitsCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
 
  // Add collection to the event
  hitsCollectionOfThisEvent->AddHitsCollection(fHitsCollectionID, fHitsCollection);
 
}

G4bool BeamTestSiliconMonitor::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
	G4ParticleDefinition* pd = aStep->GetTrack()->GetDefinition() ;
	
	//if ( postStepPoint->GetStepStatus() == fGeomBoundary /*&& (pd->GetParticleName() == "e-")*/)
	{
		// Create Hit 
		BeamTestSiliconMonitorHit* aHit = new BeamTestSiliconMonitorHit();
		
		// Get Transportaion Matrix
		G4TouchableHistory* theTouchable = (G4TouchableHistory*)(postStepPoint->GetTouchable());
		
		G4ThreeVector worldPosition = postStepPoint->GetPosition();
		G4ThreeVector localPosition = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPosition);
		
		G4ThreeVector momentumD = postStepPoint->GetMomentumDirection();
		G4ThreeVector momentum = aStep->GetTrack()->GetMomentum();
		
		//G4cout << "test that it does print out z values: " << worldPosition.z() << G4endl;

		//G4ParticleDefinition* pd = aStep->GetTrack()->GetDefinition();
		G4double ke = postStepPoint->GetKineticEnergy();

		aHit->SetExitDefinition(pd);
		aHit->SetExitKineticEnergy(ke);
		//aHit->SetExitPosition(localPosition);
		aHit->SetExitPosition(worldPosition);
		aHit->SetExitMomentumDirection(momentumD);
		aHit->SetExitMomentum(momentum);
        //This is the chamber number
        G4int replica = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber( 0 );
		aHit->SetChamberNumber(replica);
        //This is the particle id of this track (i.e. ==1 for beam)
        aHit->SetTrackId( aStep->GetTrack()->GetTrackID() );
        //Energy deposited
        aHit->SetEDep( aStep->GetTotalEnergyDeposit() );
        aHit->SetStepLength( aStep->GetStepLength() );
        aHit->SetStepStatus( postStepPoint->GetStepStatus() );
        
		fHitsCollection->insert( aHit );
	}

    return true;
}

void BeamTestSiliconMonitor::EndOfEvent(G4HCofThisEvent*) {}
