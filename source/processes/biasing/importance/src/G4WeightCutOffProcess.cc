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
// $Id: G4WeightCutOffProcess.cc 88804 2015-03-10 17:12:21Z gcosmo $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4WeightCutOffProcess.cc
//
// ----------------------------------------------------------------------

#include "G4WeightCutOffProcess.hh"
//#include "G4VScorer.hh"
#include "G4GeometryCellStep.hh"
//#include "G4GCellFinder.hh"
#include "G4TouchableHandle.hh"
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


G4WeightCutOffProcess::
G4WeightCutOffProcess(G4double wsurvival,
                      G4double wlimit,
                      G4double isource,
                      G4VIStore *istore,
		      //                      const G4VGCellFinder &aGCellFinder,
                      const G4String &aName, G4bool para)
  : G4VProcess(aName), 
    fParticleChange(new G4ParticleChange),
    fWeightSurvival(wsurvival),
    fWeightLimit(wlimit),
    fSourceImportance(isource), 
    fIStore(istore),
    //    fGCellFinder(aGCellFinder),
    fGhostWorldName("NoParallelWorld"), fGhostWorld(0),
    fGhostNavigator(0), fNavigatorID(-1), fFieldTrack('0'),
    fParaflag(para), fEndTrack('0'), feLimited(kDoNot)
{
  if (!fParticleChange)
  {
    G4Exception("G4WeightCutOffProcess::G4WeightCutOffProcess()",
                "FatalError", FatalException,
                "Failed to allocate G4ParticleChange !");
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
}

G4WeightCutOffProcess::~G4WeightCutOffProcess()
{
  delete fParticleChange;
  //  delete fGhostStep;
}


//------------------------------------------------------
//
// SetParallelWorld 
//
//------------------------------------------------------
void G4WeightCutOffProcess::
SetParallelWorld(const G4String &parallelWorldName)
{
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Get pointers of the parallel world and its navigator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  fGhostWorldName = parallelWorldName;
  fGhostWorld = fTransportationManager->GetParallelWorld(fGhostWorldName);
  fGhostNavigator = fTransportationManager->GetNavigator(fGhostWorld);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
}

void G4WeightCutOffProcess::
SetParallelWorld(G4VPhysicalVolume* parallelWorld)
{
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Get pointer of navigator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  fGhostWorldName = parallelWorld->GetName();
  fGhostWorld = parallelWorld;
  fGhostNavigator = fTransportationManager->GetNavigator(fGhostWorld);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
}

//------------------------------------------------------
//
// StartTracking
//
//------------------------------------------------------
void G4WeightCutOffProcess::StartTracking(G4Track* trk)
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
	G4Exception("G4WeightCutOffProcess::StartTracking",
		    "ProcParaWorld000",FatalException,
		    "G4WeightCutOffProcess is used for tracking without having a parallel world assigned");
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


G4double G4WeightCutOffProcess::
PostStepGetPhysicalInteractionLength(const G4Track &,
                                     G4double, G4ForceCondition* condition)
{
//   *condition = Forced;
//   return kInfinity;

//  *condition = StronglyForced;
  *condition = Forced;
  return DBL_MAX;
}
  
G4VParticleChange * 
G4WeightCutOffProcess::PostStepDoIt(const G4Track& aTrack, 
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

  if(fParaflag) {
    G4GeometryCell postCell(*(fGhostPostStepPoint->GetPhysicalVolume()), 
			    fGhostPostStepPoint->GetTouchable()->GetReplicaNumber());


    //  G4GeometryCell postCell = fGCellFinder.GetPostGeometryCell(aStep);
    //  G4GeometryCell postCell = fGCellFinder.GetPostGeometryCell(fGhostStep);
    G4double R = fSourceImportance;
    if (fIStore)
      {
	G4double i = fIStore->GetImportance(postCell);
	if (i>0)
	  {
	    R/=i;
	  }
      }
    G4double w = aTrack.GetWeight();
    if (w<R*fWeightLimit)
      {
	G4double ws = fWeightSurvival*R;
	G4double p = w/(ws);
	if (G4UniformRand()<p)
	  {
	    fParticleChange->ProposeTrackStatus(fStopAndKill);
	  }
	else
	  {
	    fParticleChange->ProposeWeight(ws);
	  }                  
      }
  } else {

    G4GeometryCell postCell(*(aStep.GetPostStepPoint()->GetPhysicalVolume()), 
			    aStep.GetPostStepPoint()->GetTouchable()->GetReplicaNumber());

    //  G4GeometryCell postCell = fGCellFinder.GetPostGeometryCell(aStep);
    //  G4GeometryCell postCell = fGCellFinder.GetPostGeometryCell(fGhostStep);
    G4double R = fSourceImportance;
    if (fIStore)
      {
	G4double i = fIStore->GetImportance(postCell);
	if (i>0)
	  {
	    R/=i;
	  }
      }
    G4double w = aTrack.GetWeight();
    if (w<R*fWeightLimit)
      {
	G4double ws = fWeightSurvival*R;
	G4double p = w/(ws);
	if (G4UniformRand()<p)
	  {
	    fParticleChange->ProposeTrackStatus(fStopAndKill);
	  }
	else
	  {
	    fParticleChange->ProposeWeight(ws);
	  }                  
      }
  }

  return fParticleChange;

}

const G4String &G4WeightCutOffProcess::GetName() const
{
  return theProcessName;
}

G4double G4WeightCutOffProcess::
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
	    // Track is not on the boundary
	    fOnBoundary = false;
	    fGhostSafety = fGhostNavigator->ComputeSafety(fEndTrack.GetPosition());
	  }
	else
	  {
	    // Track is on the boundary
	    fOnBoundary = true;
	    proposedSafety = fGhostSafety;
	  }
	//xbug?	proposedSafety = fGhostSafety;
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
    //not sensible!    return -1.0;
  }

}

  
G4double G4WeightCutOffProcess::
AtRestGetPhysicalInteractionLength(const G4Track& ,
                                   G4ForceCondition*)
{
  return -1.0;
}
  
G4VParticleChange*
G4WeightCutOffProcess::AtRestDoIt(const G4Track&, const G4Step&)
{
  return 0;
}

G4VParticleChange*
G4WeightCutOffProcess::AlongStepDoIt(const G4Track& track, const G4Step&)
{
  // Dummy ParticleChange ie: does nothing
  // Expecting G4Transportation to move the track
  pParticleChange->Initialize(track);
  return pParticleChange; 

  //  return 0;
}

void G4WeightCutOffProcess::CopyStep(const G4Step & step)
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
