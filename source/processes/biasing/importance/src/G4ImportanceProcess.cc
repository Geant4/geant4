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
// $Id: G4ImportanceProcess.cc 88804 2015-03-10 17:12:21Z gcosmo $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ImportanceProcess.cc
//
// ----------------------------------------------------------------------

#include "G4ImportanceProcess.hh"
#include "G4VImportanceAlgorithm.hh"
#include "G4GeometryCell.hh"
#include "G4SamplingPostStepAction.hh"
#include "G4VTrackTerminator.hh"
#include "G4VIStore.hh"

#include "G4Step.hh"
#include "G4Navigator.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleChange.hh"
#include "G4PathFinder.hh"
#include "G4TransportationManager.hh"
#include "G4StepPoint.hh"
#include "G4FieldTrackUpdator.hh"


G4ImportanceProcess::
G4ImportanceProcess(const G4VImportanceAlgorithm &aImportanceAlgorithm,
                        const G4VIStore &aIstore,
                        const G4VTrackTerminator *TrackTerminator,
                        const G4String &aName, G4bool para)
 : G4VProcess(aName),
   fParticleChange(new G4ParticleChange),
   fImportanceAlgorithm(aImportanceAlgorithm),
   fIStore(aIstore),
   fPostStepAction(0),
   fGhostWorldName("NoParallelWorld"),fGhostWorld(0),
   fGhostNavigator(0), fNavigatorID(-1), fFieldTrack('0'),
   fParaflag(para), fEndTrack('0'), feLimited(kDoNot)  
{
  G4cout << G4endl << G4endl << G4endl;
  G4cout << "G4ImportanceProcess:: Creating " << G4endl;
  if (TrackTerminator)
  {
    fPostStepAction = new G4SamplingPostStepAction(*TrackTerminator);
  }
  else
  {
    fPostStepAction = new G4SamplingPostStepAction(*this);
  }
  if (!fParticleChange)
  {
    G4Exception("G4ImportanceProcess::G4ImportanceProcess()",
                "FatalError", FatalException,
                "Failed allocation of G4ParticleChange !");
  }
  G4VProcess::pParticleChange = fParticleChange;


  fGhostStep = new G4Step();
  fGhostPreStepPoint = fGhostStep->GetPreStepPoint();
  fGhostPostStepPoint = fGhostStep->GetPostStepPoint();

  fTransportationManager = G4TransportationManager::GetTransportationManager();
  fPathFinder = G4PathFinder::GetInstance();

  if (verboseLevel>0)
  {
    G4cout << GetProcessName() << " is created " << G4endl;
  }

  G4cout << "G4ImportanceProcess:: importance process paraflag is: " << fParaflag << G4endl;

}

G4ImportanceProcess::~G4ImportanceProcess()
{

  delete fPostStepAction;
  delete fParticleChange;
  //  delete fGhostStep;
  // delete fGhostWorld;
  // delete fGhostNavigator;

}



//------------------------------------------------------
//
// SetParallelWorld 
//
//------------------------------------------------------
void G4ImportanceProcess::
SetParallelWorld(const G4String &parallelWorldName)
{
  G4cout << G4endl << G4endl << G4endl;
  G4cout << "G4ImportanceProcess:: SetParallelWorld name = " << parallelWorldName << G4endl;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Get pointers of the parallel world and its navigator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  fGhostWorldName = parallelWorldName;
  fGhostWorld = fTransportationManager->GetParallelWorld(fGhostWorldName);
  fGhostNavigator = fTransportationManager->GetNavigator(fGhostWorld);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
}

// void G4ImportanceProcess::
// SetParallelWorld(const G4VPhysicalVolume* parallelWorld)
// {
// //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// // Get pointer of navigator
// //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//   // G4cout << " G4ImportanceProcess:: Got here 1 " << G4endl;
//   // fGhostWorldName = parallelWorld->GetName();
//   // G4cout << " G4ImportanceProcess:: Got here 2 ghostName:" << fGhostWorldName << G4endl;
//   fGhostWorld = parallelWorld;
//   G4cout << " G4ImportanceProcess:: Got here 3 " << G4endl;
//   fGhostNavigator = fTransportationManager->GetNavigator(parallelWorld);
//   G4cout << " G4ImportanceProcess:: Got here 4 " << G4endl;
// //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// }

//------------------------------------------------------
//
// StartTracking
//
//------------------------------------------------------
void G4ImportanceProcess::StartTracking(G4Track* trk)
{
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Activate navigator and get the navigator ID
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// G4cout << " G4ParallelWorldScoringProcess::StartTracking" << G4endl;

  if(fParaflag) {
    if(fGhostNavigator)
      { fNavigatorID = fTransportationManager->ActivateNavigator(fGhostNavigator); }
    else
      {
	G4Exception("G4ImportanceProcess::StartTracking",
		    "ProcParaWorld000",FatalException,
		    "G4ImportanceProcess is used for tracking without having a parallel world assigned");
      }
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// G4cout << "G4ParallelWorldScoringProcess::StartTracking <<<<<<<<<<<<<<<<<< " << G4endl;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Let PathFinder initialize
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fPathFinder->PrepareNewTrack(trk->GetPosition(),trk->GetMomentumDirection());
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Setup initial touchables for the first step
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fOldGhostTouchable = fPathFinder->CreateTouchableHandle(fNavigatorID);
    fGhostPreStepPoint->SetTouchableHandle(fOldGhostTouchable);
    fNewGhostTouchable = fOldGhostTouchable;
    fGhostPostStepPoint->SetTouchableHandle(fNewGhostTouchable);

  // Initialize 
    fGhostSafety = -1.;
    fOnBoundary = false;
  }

}


G4double G4ImportanceProcess::
PostStepGetPhysicalInteractionLength(const G4Track& ,
                                     G4double   ,
                                     G4ForceCondition* condition)
{
//   *condition = Forced;
//   return kInfinity;

//  *condition = StronglyForced;
  *condition = Forced;
  return DBL_MAX;
}
  
G4VParticleChange *
G4ImportanceProcess::PostStepDoIt(const G4Track &aTrack,
                                      const G4Step &aStep)
{
  fParticleChange->Initialize(aTrack);

  if(fParaflag) {
    fOldGhostTouchable = fGhostPostStepPoint->GetTouchableHandle();
    //xbug?    fOnBoundary = false;
    CopyStep(aStep);
    
    if(fOnBoundary)
      {
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Locate the point and get new touchable
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //??  fPathFinder->Locate(step.GetPostStepPoint()->GetPosition(),
  //??                      step.GetPostStepPoint()->GetMomentumDirection());
	fNewGhostTouchable = fPathFinder->CreateTouchableHandle(fNavigatorID);
	//AH	G4cout << " on boundary " << G4endl;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      }
    else
      {
	// Do I need this ??????????????????????????????????????????????????????????
	// fGhostNavigator->LocateGlobalPointWithinVolume(track.GetPosition());
	// ?????????????????????????????????????????????????????????????????????????
	
	// fPathFinder->ReLocate(track.GetPosition());
	
	// reuse the touchable
	fNewGhostTouchable = fOldGhostTouchable;
	//AH	G4cout << " NOT on boundary " << G4endl;
      }
    
    fGhostPreStepPoint->SetTouchableHandle(fOldGhostTouchable);
    fGhostPostStepPoint->SetTouchableHandle(fNewGhostTouchable);
    
//   if ( (aStep.GetPostStepPoint()->GetStepStatus() == fGeomBoundary)
//     && (aStep.GetStepLength() > kCarTolerance) )
//   {
//     if (aTrack.GetTrackStatus()==fStopAndKill)
//     {
//       G4cout << "WARNING - G4ImportanceProcess::PostStepDoIt()"
//              << "          StopAndKill track." << G4endl;
//     }

//     G4StepPoint *prepoint  = aStep.GetPreStepPoint();
//     G4StepPoint *postpoint = aStep.GetPostStepPoint();
//     G4GeometryCell prekey(*(prepoint->GetPhysicalVolume()), 
//                          prepoint->GetTouchable()->GetReplicaNumber());
//     G4GeometryCell postkey(*(postpoint->GetPhysicalVolume()), 
//                           postpoint->GetTouchable()->GetReplicaNumber());


    if ( (fGhostPostStepPoint->GetStepStatus() == fGeomBoundary)
	 && (aStep.GetStepLength() > kCarTolerance) )
      {
	if (aTrack.GetTrackStatus()==fStopAndKill)
	  {
	    G4cout << "WARNING - G4ImportanceProcess::PostStepDoIt()"
		   << "          StopAndKill track. on boundary" << G4endl;
	  }
	
	G4GeometryCell prekey(*(fGhostPreStepPoint->GetPhysicalVolume()), 
			      fGhostPreStepPoint->GetTouchable()->GetReplicaNumber());
	G4GeometryCell postkey(*(fGhostPostStepPoint->GetPhysicalVolume()), 
			       fGhostPostStepPoint->GetTouchable()->GetReplicaNumber());
	
	//AH
	/*
	G4cout << G4endl;
	G4cout << G4endl;
	G4cout << " inside parallel importance process " << aTrack.GetCurrentStepNumber() << G4endl;
	G4cout << G4endl;
	G4cout << G4endl;
	G4cout << " prekey: " << fGhostPreStepPoint->GetPhysicalVolume()->GetName() << " replica:" 
	       << fGhostPreStepPoint->GetTouchable()->GetReplicaNumber() << G4endl;
	G4cout << " prekey ISTORE: " << fIStore.GetImportance(prekey) << G4endl;
	G4cout << " postkey: " << G4endl;
	G4cout << " postkey ISTORE: " << fIStore.GetImportance(postkey) << G4endl;
	*/
	//AH
	G4Nsplit_Weight nw = fImportanceAlgorithm.
	  Calculate(fIStore.GetImportance(prekey),
		    fIStore.GetImportance(postkey), 
		    aTrack.GetWeight());
	//AH
	/*
	G4cout << " prekey weight: " << fIStore.GetImportance(prekey) 
	       << " postkey weight: " << fIStore.GetImportance(postkey) 
	       << " split weight: " << nw << G4endl;
	G4cout << " before poststepaction " << G4endl;
	*/
	//AH
	fPostStepAction->DoIt(aTrack, fParticleChange, nw);
	//AH	G4cout << " after post step do it " << G4endl;
      }
  } else {
    if ( (aStep.GetPostStepPoint()->GetStepStatus() == fGeomBoundary)
	 && (aStep.GetStepLength() > kCarTolerance) )
      {
	//AH	G4cout << " inside non-parallel importance process " << G4endl;
	if (aTrack.GetTrackStatus()==fStopAndKill)
	  {
	    G4cout << "WARNING - G4ImportanceProcess::PostStepDoIt()"
		   << "          StopAndKill track. on boundary non-parallel" << G4endl;
	  }
	
	G4StepPoint *prepoint  = aStep.GetPreStepPoint();
	G4StepPoint *postpoint = aStep.GetPostStepPoint();
	
	G4GeometryCell prekey(*(prepoint->GetPhysicalVolume()), 
			      prepoint->GetTouchable()->GetReplicaNumber());
	G4GeometryCell postkey(*(postpoint->GetPhysicalVolume()), 
			       postpoint->GetTouchable()->GetReplicaNumber());
	
	G4Nsplit_Weight nw = fImportanceAlgorithm.
	  Calculate(fIStore.GetImportance(prekey),
		    fIStore.GetImportance(postkey), 
		    aTrack.GetWeight());
	//AH
	/*
	G4cout << " prekey weight: " << fIStore.GetImportance(prekey) 
	       << " postkey weight: " << fIStore.GetImportance(postkey) 
	       << " split weight: " << nw << G4endl;
	G4cout << " before poststepaction 2 " << G4endl;
	*/
	//AH
	fPostStepAction->DoIt(aTrack, fParticleChange, nw);
	//AH	G4cout << " after poststepaction 2 " << G4endl;
      }
  }
  return fParticleChange;
}

void G4ImportanceProcess::KillTrack() const
{
  fParticleChange->ProposeTrackStatus(fStopAndKill);
}

const G4String &G4ImportanceProcess::GetName() const
{
  return theProcessName;
}

G4double G4ImportanceProcess::
AlongStepGetPhysicalInteractionLength(
            const G4Track& track, G4double  previousStepSize, G4double  currentMinimumStep,
            G4double& proposedSafety, G4GPILSelection* selection)
{
  if(fParaflag) {    
    *selection = NotCandidateForSelection;
    G4double returnedStep = DBL_MAX;
    
    if (previousStepSize > 0.)
      { fGhostSafety -= previousStepSize; }
    //  else
    //  { fGhostSafety = -1.; }
    if (fGhostSafety < 0.) fGhostSafety = 0.0;
    
    // ------------------------------------------
    // Determination of the proposed STEP LENGTH:
    // ------------------------------------------
    if (currentMinimumStep <= fGhostSafety && currentMinimumStep > 0.)
      {
	// I have no chance to limit
	returnedStep = currentMinimumStep;
	fOnBoundary = false;
	proposedSafety = fGhostSafety - currentMinimumStep;
	//AH	G4cout << " step not limited, why? " << G4endl;
      }
    else // (currentMinimumStep > fGhostSafety: I may limit the Step)
      {
	G4FieldTrackUpdator::Update(&fFieldTrack,&track);
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// ComputeStep
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	returnedStep
	  = fPathFinder->ComputeStep(fFieldTrack,currentMinimumStep,fNavigatorID,
				     track.GetCurrentStepNumber(),fGhostSafety,feLimited,
				     fEndTrack,track.GetVolume());
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if(feLimited == kDoNot)
	  {
	    //AH	    G4cout << " computestep came back with not-boundary " << G4endl;
	    // Track is not on the boundary
	    fOnBoundary = false;
	    fGhostSafety = fGhostNavigator->ComputeSafety(fEndTrack.GetPosition());
	  }
	else
	  {
	    // Track is on the boundary
	    //AH	    G4cout << " FOUND A BOUNDARY ! " << track.GetCurrentStepNumber() << G4endl;
	    fOnBoundary = true;
	    // proposedSafety = fGhostSafety;
	  }
	proposedSafety = fGhostSafety;
	if(feLimited == kUnique || feLimited == kSharedOther) {
	  *selection = CandidateForSelection;
	}else if (feLimited == kSharedTransport) { 
	  returnedStep *= (1.0 + 1.0e-9);  
	  // Expand to disable its selection in Step Manager comparison
	}
	
      }

  // ----------------------------------------------
  // Returns the fGhostSafety as the proposedSafety
  // The SteppingManager will take care of keeping
  // the smallest one.
  // ----------------------------------------------
    return returnedStep;

  } else {

    return DBL_MAX;

  }

}

G4double G4ImportanceProcess::
AtRestGetPhysicalInteractionLength(const G4Track& ,
                                   G4ForceCondition*) 
{
  return -1.0;
}
  
G4VParticleChange* G4ImportanceProcess::
AtRestDoIt(const G4Track&, const G4Step&) 
{
  return 0;
}

G4VParticleChange* G4ImportanceProcess::
AlongStepDoIt(const G4Track& aTrack, const G4Step& )
{
  // Dummy ParticleChange ie: does nothing
  // Expecting G4Transportation to move the track
  //AH  G4cout << " along step do it " << G4endl;
  pParticleChange->Initialize(aTrack);
  return pParticleChange; 
  //  return 0;
}

void G4ImportanceProcess::CopyStep(const G4Step & step)
{
  fGhostStep->SetTrack(step.GetTrack());
  fGhostStep->SetStepLength(step.GetStepLength());
  fGhostStep->SetTotalEnergyDeposit(step.GetTotalEnergyDeposit());
  fGhostStep->SetControlFlag(step.GetControlFlag());

  *fGhostPreStepPoint = *(step.GetPreStepPoint());
  *fGhostPostStepPoint = *(step.GetPostStepPoint());

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Set StepStatus for ghost world
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if(fOnBoundary)
  { fGhostPostStepPoint->SetStepStatus(fGeomBoundary); }
  else if(fGhostPostStepPoint->GetStepStatus()==fGeomBoundary)
  { fGhostPostStepPoint->SetStepStatus(fPostStepDoItProc); }
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
}
