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
// $Id: G4ScoreProcess.cc,v 1.1 2006/11/20 10:02:18 ahoward Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ScoreProcess.cc
//
// ----------------------------------------------------------------------

#include "G4ScoreProcess.hh"
#include "G4VScorer.hh"

#include "G4GeometryCellStep.hh"

#include "G4Step.hh"
#include "G4Navigator.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleChange.hh"
#include "G4PathFinder.hh"
#include "G4TransportationManager.hh"
#include "G4StepPoint.hh"
#include "G4FieldTrackUpdator.hh"


G4ScoreProcess::G4ScoreProcess(G4VScorer &aScorer,
                                 const G4String &aName, G4bool para)
 : G4VProcess(aName), 
   fParticleChange(new G4ParticleChange),
   fScorer(aScorer),
   fKillTrack(false),
   fGhostNavigator(0), fNavigatorID(-1), fFieldTrack('0')
{

  if (!fParticleChange)
  {
    G4Exception("G4ScoreProcess::G4ScoreProcess()", "Fatal Error",
                FatalException, "Failed to allocate G4ParticleChange !");
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

  paraflag = para;


}

G4ScoreProcess::~G4ScoreProcess()
{
  delete fParticleChange;
  delete fGhostStep;
}



//------------------------------------------------------
//
// SetParallelWorld 
//
//------------------------------------------------------
void G4ScoreProcess::
SetParallelWorld(G4String parallelWorldName)
{
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Get pointers of the parallel world and its navigator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  fGhostWorldName = parallelWorldName;
  fGhostWorld = fTransportationManager->GetParallelWorld(fGhostWorldName);
  fGhostNavigator = fTransportationManager->GetNavigator(fGhostWorld);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
}

void G4ScoreProcess::
SetParallelWorld(G4VPhysicalVolume* parallelWorld)
{
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Get pointer of navigator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  G4cout << " Getting Name: " << G4endl;
  fGhostWorldName = parallelWorld->GetName();
  G4cout << " Get Name :" << fGhostWorldName << G4endl;
  fGhostWorld = parallelWorld;
  fGhostNavigator = fTransportationManager->GetNavigator(fGhostWorld);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
}


//------------------------------------------------------
//
// StartTracking
//
//------------------------------------------------------
void G4ScoreProcess::StartTracking(G4Track* trk)
{
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Activate navigator and get the navigator ID
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// G4cout << " G4ParallelWorldScoringProcess::StartTracking" << G4endl;

  if(paraflag) {
    if(fGhostNavigator)
      { fNavigatorID = fTransportationManager->ActivateNavigator(fGhostNavigator); }
    else
      {
	G4Exception("G4ScoreProcess::StartTracking",
		    "ProcParaWorld000",FatalException,
		    "G4ScoreProcess is used for tracking without having a parallel world assigned");
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



G4double G4ScoreProcess::
PostStepGetPhysicalInteractionLength(const G4Track&,
                                     G4double  ,
                                     G4ForceCondition* condition)
{
//   *condition = Forced;
//   return kInfinity;

//  *condition = StronglyForced;
  *condition = Forced;
  return DBL_MAX;
}
  
G4VParticleChange * 
G4ScoreProcess::PostStepDoIt(const G4Track& aTrack, const G4Step &aStep)
{

  //AH  G4cout << " inside scoring post-step, track: " << aTrack.GetTrackID() << " step: " << aTrack.GetCurrentStepNumber() << G4endl;
  //AH  G4cout << " and fOnBoundary " << fOnBoundary << G4endl;
  fParticleChange->Initialize(aTrack);

  if(paraflag) {
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
	//AH	G4cout << " GOT HERE? inside post step " << G4endl;
	fNewGhostTouchable = fPathFinder->CreateTouchableHandle(fNavigatorID);
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
      }
        
    fGhostPreStepPoint->SetTouchableHandle(fOldGhostTouchable);
    fGhostPostStepPoint->SetTouchableHandle(fNewGhostTouchable);
  }
    

	//     G4StepPoint *prepoint = aStep.GetPreStepPoint();
	//     G4StepPoint *postpoint = aStep.GetPostStepPoint();
	
	//     G4GeometryCell prekey(*(prepoint->GetPhysicalVolume()), 
	//                            prepoint->GetTouchable()->GetReplicaNumber());
	//     G4GeometryCell postkey(*(postpoint->GetPhysicalVolume()), 
//                             postpoint->GetTouchable()->GetReplicaNumber());

  if(paraflag) {
    if (fGhostStep->GetStepLength() > kCarTolerance)
      {
	G4GeometryCell prekey(*(fGhostPreStepPoint->GetPhysicalVolume()), 
			      fGhostPreStepPoint->GetTouchable()->GetReplicaNumber());
	G4GeometryCell postkey(*(fGhostPostStepPoint->GetPhysicalVolume()), 
				 fGhostPostStepPoint->GetTouchable()->GetReplicaNumber());
	G4GeometryCellStep pstep(prekey, postkey);
	pstep.SetCrossBoundary(false);
	//xxtest	  pstep.SetCrossBoundary(true);
	
	if (prekey != postkey)
	  {
	    pstep.SetCrossBoundary(true);
	  } 
	fScorer.Score(aStep, pstep); 
      } else {	  
	if (aStep.GetStepLength() > kCarTolerance) 
	  {
	    G4StepPoint *prepoint = aStep.GetPreStepPoint();
	    G4StepPoint *postpoint = aStep.GetPostStepPoint();
	    
	    G4GeometryCell prekey(*(prepoint->GetPhysicalVolume()), 
				  prepoint->GetTouchable()->GetReplicaNumber());
	    G4GeometryCell postkey(*(postpoint->GetPhysicalVolume()), 
				   postpoint->GetTouchable()->GetReplicaNumber());
	    G4GeometryCellStep pstep(prekey, postkey);
	    pstep.SetCrossBoundary(false);
	    
	    if (prekey != postkey)
	      {
		pstep.SetCrossBoundary(true);
	      } 
	    fScorer.Score(aStep, pstep); 
	  }
      }
  }  

  if (fKillTrack)
    {
      //AH      G4cout << " score process is killing track - why? " << G4endl;
      fKillTrack = false;
      fParticleChange->ProposeTrackStatus(fStopAndKill);
    }
  
  return fParticleChange;
  //  return G4VProcess::pParticleChange;
}

const G4String &G4ScoreProcess::GetName() const
{
  return theProcessName;
}


G4double G4ScoreProcess::
AlongStepGetPhysicalInteractionLength(
            const G4Track& track, G4double  previousStepSize, G4double  currentMinimumStep,
            G4double& proposedSafety, G4GPILSelection* selection)
{
  //AH  G4cout << " inside scoring along-step GPIL, track : " << track.GetTrackID() << " step#: " << track.GetCurrentStepNumber() << G4endl;
  if(paraflag) {
    static G4FieldTrack endTrack('0');
    static ELimited eLimited;
    
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
      }
    else // (currentMinimumStep > fGhostSafety: I may limit the Step)
      {
	G4FieldTrackUpdator::Update(&fFieldTrack,&track);
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// ComputeStep
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	returnedStep
	  = fPathFinder->ComputeStep(fFieldTrack,currentMinimumStep,fNavigatorID,
				     track.GetCurrentStepNumber(),fGhostSafety,eLimited,
				     endTrack,track.GetVolume());
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	if(eLimited == kDoNot)
	  {
	    // Track is not on the boundary
	    fOnBoundary = false;
	    fGhostSafety = fGhostNavigator->ComputeSafety(endTrack.GetPosition());
	  }
	else
	  {
	    // Track is on the boundary
	    //AH	    G4cout << " GOT HERE " << G4endl;
	    fOnBoundary = true;
	    proposedSafety = fGhostSafety;
	  }
	//xbug?	proposedSafety = fGhostSafety;
	if(eLimited == kUnique || eLimited == kSharedOther) {
	  *selection = CandidateForSelection;
	}else if (eLimited == kSharedTransport) { 
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
    // not sensible!    return -1.0;
  }

}

G4double G4ScoreProcess::
AtRestGetPhysicalInteractionLength(const G4Track&,
                                   G4ForceCondition*)
{
  return -1.0;
}

G4VParticleChange* G4ScoreProcess::AtRestDoIt(const G4Track&,
                                               const G4Step&)
{
  return 0;
}

G4VParticleChange* G4ScoreProcess::AlongStepDoIt(const G4Track& track,
                                                  const G4Step&)
{
  // Dummy ParticleChange ie: does nothing
  // Expecting G4Transportation to move the track
  pParticleChange->Initialize(track);
  return pParticleChange; 

  //  return 0;
}

void G4ScoreProcess::KillTrack() const
{
  fKillTrack = true;
}

void G4ScoreProcess::CopyStep(const G4Step & step)
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
