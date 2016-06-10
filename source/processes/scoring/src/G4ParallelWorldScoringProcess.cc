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
// $Id: G4ParallelWorldScoringProcess.cc 88349 2015-02-16 08:47:16Z gcosmo $
//
//

#include "G4ParallelWorldScoringProcess.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4Step.hh"
#include "G4Navigator.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleChange.hh"
#include "G4PathFinder.hh"
#include "G4TransportationManager.hh"
#include "G4ParticleChange.hh"
#include "G4StepPoint.hh"
#include "G4FieldTrackUpdator.hh"
#include "G4ParticleDefinition.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"

//--------------------------------
// Constructor with name and type:
//--------------------------------
G4ParallelWorldScoringProcess::
G4ParallelWorldScoringProcess(const G4String& processName,G4ProcessType theType)
:G4VProcess(processName,theType), fGhostNavigator(0), fNavigatorID(-1), fFieldTrack('0')
{
  pParticleChange = &aDummyParticleChange;

  fGhostStep = new G4Step();
  fGhostPreStepPoint = fGhostStep->GetPreStepPoint();
  fGhostPostStepPoint = fGhostStep->GetPostStepPoint();

  fTransportationManager = G4TransportationManager::GetTransportationManager();
  fPathFinder = G4PathFinder::GetInstance();

  fGhostWorld = 0;
  fGhostSafety = 0.;
  fOnBoundary = false;

  if (verboseLevel>0)
  {
    G4cout << GetProcessName() << " is created " << G4endl;
  }
}

// -----------
// Destructor:
// -----------
G4ParallelWorldScoringProcess::~G4ParallelWorldScoringProcess()
{
  delete fGhostStep;
}

//------------------------------------------------------
//
// SetParallelWorld 
//
//------------------------------------------------------
void G4ParallelWorldScoringProcess::
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

void G4ParallelWorldScoringProcess::
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

G4bool G4ParallelWorldScoringProcess::
IsAtRestRequired(G4ParticleDefinition* partDef)
{
  G4int pdgCode = partDef->GetPDGEncoding();
  if(pdgCode==0)
  {
    G4String partName = partDef->GetParticleName();
    if(partName=="opticalphoton") return false;
    if(partName=="geantino") return false;
    if(partName=="chargedgeantino") return false;
  }
  else
  {
    if(pdgCode==22) return false; // gamma
    if(pdgCode==11) return false; // electron
    if(pdgCode==2212) return false; // proton
    if(pdgCode==-12) return false; // anti_nu_e
    if(pdgCode==12) return false; // nu_e
    if(pdgCode==-14) return false; // anti_nu_mu
    if(pdgCode==14) return false; // nu_mu
    if(pdgCode==-16) return false; // anti_nu_tau
    if(pdgCode==16) return false; // nu_tau
  }
  return true;
}


//------------------------------------------------------
//
// StartTracking
//
//------------------------------------------------------
void G4ParallelWorldScoringProcess::StartTracking(G4Track* trk)
{
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Activate navigator and get the navigator ID
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// G4cout << " G4ParallelWorldScoringProcess::StartTracking" << G4endl;
  if(fGhostNavigator)
  { fNavigatorID = fTransportationManager->ActivateNavigator(fGhostNavigator); }
  else
  {
    G4Exception("G4ParallelWorldScoringProcess::StartTracking",
       "ProcParaWorld000",FatalException,
       "G4ParallelWorldScoringProcess is used for tracking without having a parallel world assigned");
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
  fGhostPreStepPoint->SetStepStatus(fUndefined);
  fGhostPostStepPoint->SetStepStatus(fUndefined);
}

//----------------------------------------------------------
//
//  AtRestGetPhysicalInteractionLength()
//
//----------------------------------------------------------
G4double 
G4ParallelWorldScoringProcess::AtRestGetPhysicalInteractionLength(
         const G4Track& /*track*/, 
         G4ForceCondition* condition)
{
  *condition = Forced;
  return DBL_MAX;
}

//------------------------------------
//
//             AtRestDoIt()
//
//------------------------------------
G4VParticleChange* G4ParallelWorldScoringProcess::AtRestDoIt(
     const G4Track& track,
     const G4Step& step)
{ 
  fOldGhostTouchable = fGhostPostStepPoint->GetTouchableHandle();
  G4VSensitiveDetector* aSD = 0;
  if(fOldGhostTouchable->GetVolume())
  { aSD = fOldGhostTouchable->GetVolume()->GetLogicalVolume()->GetSensitiveDetector(); }
  fOnBoundary = false;
  CopyStep(step);
  fGhostPreStepPoint->SetSensitiveDetector(aSD);

  fNewGhostTouchable = fOldGhostTouchable;
  
  fGhostPreStepPoint->SetTouchableHandle(fOldGhostTouchable);
  fGhostPostStepPoint->SetTouchableHandle(fNewGhostTouchable);
  if(fNewGhostTouchable->GetVolume())
  {
    fGhostPostStepPoint->SetSensitiveDetector(
      fNewGhostTouchable->GetVolume()->GetLogicalVolume()->GetSensitiveDetector());
  }
  else
  { fGhostPostStepPoint->SetSensitiveDetector(0); }

  if (verboseLevel>1) Verbose(step);

  G4VSensitiveDetector* sd = fGhostPreStepPoint->GetSensitiveDetector();
  if(sd)
  {
    sd->Hit(fGhostStep);
  }

  pParticleChange->Initialize(track);
  return pParticleChange;
}

//----------------------------------------------------------
//
//  PostStepGetPhysicalInteractionLength()
//
//----------------------------------------------------------
G4double 
G4ParallelWorldScoringProcess::PostStepGetPhysicalInteractionLength(
         const G4Track& /*track*/, 
         G4double   /*previousStepSize*/, 
         G4ForceCondition* condition)
{
  // I must be invoked anyway to score the hit.
  *condition = StronglyForced;
  return DBL_MAX;
}

//------------------------------------
//
//             PostStepDoIt()
//
//------------------------------------
G4VParticleChange* G4ParallelWorldScoringProcess::PostStepDoIt(
     const G4Track& track,
     const G4Step& step)
{ 
  fOldGhostTouchable = fGhostPostStepPoint->GetTouchableHandle();
  G4VSensitiveDetector* aSD = 0;
  if(fOldGhostTouchable->GetVolume())
  { aSD = fOldGhostTouchable->GetVolume()->GetLogicalVolume()->GetSensitiveDetector(); }
  CopyStep(step);
  fGhostPreStepPoint->SetSensitiveDetector(aSD);

  //  fPathFinder->Locate( track.GetPosition(),
  //                       track.GetMomentumDirection(),
  //                       true);

  //  fPathFinder->Locate(step.GetPostStepPoint()->GetPosition(),
  //                      step.GetPostStepPoint()->GetMomentumDirection());

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

  if(fNewGhostTouchable->GetVolume())
  {
    fGhostPostStepPoint->SetSensitiveDetector(
      fNewGhostTouchable->GetVolume()->GetLogicalVolume()->GetSensitiveDetector());
  }
  else
  { fGhostPostStepPoint->SetSensitiveDetector(0); }

  if (verboseLevel>1) Verbose(step);

  G4VSensitiveDetector* sd = fGhostPreStepPoint->GetSensitiveDetector();
  if(sd)
  {
    sd->Hit(fGhostStep);
  }

  pParticleChange->Initialize(track); // Does not change the track properties
  return pParticleChange;
}


//---------------------------------------
//
//  AlongStepGetPhysicalInteractionLength
//
//---------------------------------------
G4double G4ParallelWorldScoringProcess::AlongStepGetPhysicalInteractionLength(
            const G4Track& track, G4double  previousStepSize, G4double  currentMinimumStep,
            G4double& proposedSafety, G4GPILSelection* selection)
{
  static G4ThreadLocal G4FieldTrack *endTrack_G4MT_TLS_ = 0 ; if (!endTrack_G4MT_TLS_) endTrack_G4MT_TLS_ = new  G4FieldTrack ('0') ;  G4FieldTrack &endTrack = *endTrack_G4MT_TLS_;
  static G4ThreadLocal ELimited *eLimited_G4MT_TLS_ = 0 ; if (!eLimited_G4MT_TLS_) eLimited_G4MT_TLS_ = new  ELimited  ;  ELimited &eLimited = *eLimited_G4MT_TLS_;
  
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
      fOnBoundary = true;
      // fEndGhostSafety = 0.0;  // Will apply at the end of the step ...
    }
    proposedSafety = fGhostSafety;
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
}

G4VParticleChange* G4ParallelWorldScoringProcess::AlongStepDoIt(
    const G4Track& track, const G4Step& )
{
  // Dummy ParticleChange ie: does nothing
  // Expecting G4Transportation to move the track
  pParticleChange->Initialize(track);
  return pParticleChange; 
}


void G4ParallelWorldScoringProcess::CopyStep(const G4Step & step)
{
  G4StepStatus prevStat = fGhostPostStepPoint->GetStepStatus();

  fGhostStep->SetTrack(step.GetTrack());
  fGhostStep->SetStepLength(step.GetStepLength());
  fGhostStep->SetTotalEnergyDeposit(step.GetTotalEnergyDeposit());
  fGhostStep->SetNonIonizingEnergyDeposit(step.GetNonIonizingEnergyDeposit());
  fGhostStep->SetControlFlag(step.GetControlFlag());

  *fGhostPreStepPoint = *(step.GetPreStepPoint());
  *fGhostPostStepPoint = *(step.GetPostStepPoint());

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Set StepStatus for ghost world
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  fGhostPreStepPoint->SetStepStatus(prevStat);
  if(fOnBoundary)
  { fGhostPostStepPoint->SetStepStatus(fGeomBoundary); }
  else if(fGhostPostStepPoint->GetStepStatus()==fGeomBoundary)
  { fGhostPostStepPoint->SetStepStatus(fPostStepDoItProc); }
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
}

void G4ParallelWorldScoringProcess::Verbose(const G4Step& step) const
{
  G4cout << "In mass geometry ------------------------------------------------" << G4endl;
  G4cout << " StepLength : " << step.GetStepLength()/mm << "      TotalEnergyDeposit : "
         << step.GetTotalEnergyDeposit()/MeV << G4endl;
  G4cout << " PreStepPoint : "
         << step.GetPreStepPoint()->GetPhysicalVolume()->GetName() << " - ";
  if(step.GetPreStepPoint()->GetProcessDefinedStep())
  { G4cout << step.GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName(); }
  else
  { G4cout << "NoProcessAssigned"; }
  G4cout << G4endl;
  G4cout << "                " << step.GetPreStepPoint()->GetPosition() << G4endl;
  G4cout << " PostStepPoint : ";
  if(step.GetPostStepPoint()->GetPhysicalVolume()) 
  { G4cout << step.GetPostStepPoint()->GetPhysicalVolume()->GetName(); }
  else
  { G4cout << "OutOfWorld"; }
  G4cout << " - ";
  if(step.GetPostStepPoint()->GetProcessDefinedStep())
  { G4cout << step.GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName(); }
  else
  { G4cout << "NoProcessAssigned"; }
  G4cout << G4endl;
  G4cout << "                 " << step.GetPostStepPoint()->GetPosition() << G4endl;

  G4cout << "In ghost geometry ------------------------------------------------" << G4endl;
  G4cout << " StepLength : " << fGhostStep->GetStepLength()/mm
         << "      TotalEnergyDeposit : "
         << fGhostStep->GetTotalEnergyDeposit()/MeV << G4endl;
  G4cout << " PreStepPoint : "
         << fGhostStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() << " ["
           << fGhostStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber()
           << " ]" << " - ";
  if(fGhostStep->GetPreStepPoint()->GetProcessDefinedStep())
  { G4cout << fGhostStep->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName(); }
  else
  { G4cout << "NoProcessAssigned"; }
  G4cout << G4endl;
  G4cout << "                " << fGhostStep->GetPreStepPoint()->GetPosition() << G4endl;
  G4cout << " PostStepPoint : ";
  if(fGhostStep->GetPostStepPoint()->GetPhysicalVolume()) 
  {
    G4cout << fGhostStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() << " ["
           << fGhostStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber()
           << " ]";
  }
  else
  { G4cout << "OutOfWorld"; }
  G4cout << " - ";
  if(fGhostStep->GetPostStepPoint()->GetProcessDefinedStep())
  { G4cout << fGhostStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName(); }
  else
  { G4cout << "NoProcessAssigned"; }
  G4cout << G4endl;
  G4cout << "                 " << fGhostStep->GetPostStepPoint()->GetPosition() << " == "
         << fGhostStep->GetTrack()->GetMomentumDirection() 
         << G4endl;

}

