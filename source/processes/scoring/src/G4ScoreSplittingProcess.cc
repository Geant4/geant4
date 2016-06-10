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
// $Id: G4ScoreSplittingProcess.cc 68733 2013-04-05 09:45:28Z gcosmo $
//

#include "G4ScoreSplittingProcess.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleChange.hh"
#include "G4TransportationManager.hh"
#include "G4RegularNavigationHelper.hh"
#include "G4ParticleChange.hh"
#include "G4StepPoint.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"

#include "G4EnergySplitter.hh"
#include "G4TouchableHistory.hh"

//--------------------------------
// Constructor with name and type:
//--------------------------------
G4ScoreSplittingProcess::
G4ScoreSplittingProcess(const G4String& processName,G4ProcessType theType)
  :G4VProcess(processName,theType),
   fOldTouchableH(),  fNewTouchableH(), fInitialTouchableH(), fFinalTouchableH()
{
  pParticleChange = &xParticleChange;

  fSplitStep = new G4Step();
  fSplitPreStepPoint =  fSplitStep->GetPreStepPoint();
  fSplitPostStepPoint = fSplitStep->GetPostStepPoint();

  if (verboseLevel>0)
  {
    G4cout << GetProcessName() << " is created " << G4endl;
  }
  fpEnergySplitter = new G4EnergySplitter(); 
}

// -----------
// Destructor:
// -----------
G4ScoreSplittingProcess::~G4ScoreSplittingProcess()
{
  delete fSplitStep;
  delete fpEnergySplitter; 
}

//------------------------------------------------------
//
// StartTracking
//
//------------------------------------------------------
void G4ScoreSplittingProcess::StartTracking(G4Track* trk)
{
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Setup initial touchables for the first step
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  const G4Step* pStep= trk->GetStep();

  fOldTouchableH = trk->GetTouchableHandle();
  *fSplitPreStepPoint = *(pStep->GetPreStepPoint()); // Best to copy, so as to initialise
  fSplitPreStepPoint->SetTouchableHandle(fOldTouchableH);
  fNewTouchableH = fOldTouchableH;
  *fSplitPostStepPoint= *(pStep->GetPostStepPoint()); // Best to copy, so as to initialise
  fSplitPostStepPoint->SetTouchableHandle(fNewTouchableH);

  ///  Initialize 
  fSplitPreStepPoint ->SetStepStatus(fUndefined);
  fSplitPostStepPoint->SetStepStatus(fUndefined);
}


//----------------------------------------------------------
//
//  PostStepGetPhysicalInteractionLength()
//
//----------------------------------------------------------
G4double 
G4ScoreSplittingProcess::PostStepGetPhysicalInteractionLength(
         const G4Track& /*track*/, 
         G4double   /*previousStepSize*/, 
         G4ForceCondition* condition)
{
  // This process must be invoked anyway to score the hit
  //   - to do the scoring if the current volume is a regular structure, or
  //   - else to toggle the flag so that the SteppingManager does the scoring.
  *condition = StronglyForced;

  // Future optimisation: check whether in regular structure.  
  //  If it is in regular structure,  be StronglyForced
  //  If  not  in regular structure, 
  //         ask to be called only if SteppingControl is AvoidHitInvocation
  //         in order to reset it to NormalCondition

  return DBL_MAX;
}

//------------------------------------
//
//             PostStepDoIt()
//
//------------------------------------
G4VParticleChange* G4ScoreSplittingProcess::PostStepDoIt(
     const G4Track& track,
     const G4Step& step)
{ 
  G4VPhysicalVolume*    pCurrentVolume= track.GetVolume(); 
  G4LogicalVolume*      pLogicalVolume= pCurrentVolume->GetLogicalVolume(); 
  G4VSensitiveDetector* ptrSD = pLogicalVolume->GetSensitiveDetector();

  pParticleChange->Initialize(track); 
  if(  ( ! pCurrentVolume->IsRegularStructure() ) || ( !ptrSD ) 
    || G4RegularNavigationHelper::Instance()->GetStepLengths().size() <= 1) {
     // Set the flag to make sure that Stepping Manager does the scoring
     pParticleChange->ProposeSteppingControl( NormalCondition );     
  } else { 
     G4ThreeVector preStepPosition, postStepPosition, direction, finalPostStepPosition;
     pParticleChange->ProposeSteppingControl( AvoidHitInvocation );
 
     G4double totalEnergyDeposit=  step.GetTotalEnergyDeposit();
     G4StepStatus  fullStepStatus= step.GetPostStepPoint()->GetStepStatus(); 

     CopyStepStart(step);
     fSplitPreStepPoint->SetSensitiveDetector(ptrSD);
     fOldTouchableH = fInitialTouchableH;
     fNewTouchableH=  fOldTouchableH; 
     *fSplitPostStepPoint= *(step.GetPreStepPoint()); 
     
     // Split the energy
     // ----------------
     G4int numberVoxelsInStep= fpEnergySplitter->SplitEnergyInVolumes( &step );

     preStepPosition= step.GetPreStepPoint()->GetPosition();
     finalPostStepPosition= step.GetPostStepPoint()->GetPosition(); 
     direction= (finalPostStepPosition - preStepPosition).unit(); 

     fFinalTouchableH= track.GetNextTouchableHandle(); 

     postStepPosition= preStepPosition;
     // Loop over the sub-parts of this step
     G4int iStep;
     for ( iStep=0; iStep < numberVoxelsInStep; iStep++ ){ 
        G4int     idVoxel=  -1;  // Voxel ID
        G4double  stepLength=0.0, energyLoss= 0.0;

        *fSplitPreStepPoint  = *fSplitPostStepPoint;
        fOldTouchableH = fNewTouchableH; 

        preStepPosition= postStepPosition;
        fSplitPreStepPoint->SetPosition( preStepPosition );
        fSplitPreStepPoint->SetTouchableHandle(fOldTouchableH);
        
        fpEnergySplitter->GetLengthAndEnergyDeposited( iStep, idVoxel, stepLength, energyLoss);

        // Correct the material, so that the track->GetMaterial gives correct answer
        pLogicalVolume->SetMaterial( fpEnergySplitter->GetVoxelMaterial( iStep) );  // idVoxel) ); 

        postStepPosition= preStepPosition + stepLength * direction; 
        fSplitPostStepPoint->SetPosition(postStepPosition); 

        // Load the Step with the new values
        fSplitStep->SetStepLength(stepLength);
        fSplitStep->SetTotalEnergyDeposit(energyLoss);
        if( iStep < numberVoxelsInStep -1 ){ 
          fSplitStep->GetPostStepPoint()->SetStepStatus( fGeomBoundary );
          G4int  nextVoxelId= -1;
          fpEnergySplitter->GetVoxelID( iStep+1, nextVoxelId );

          // Create new "next" touchable for each section ??
          G4VTouchable* fNewTouchablePtr= 
              CreateTouchableForSubStep( nextVoxelId, postStepPosition );
          fNewTouchableH= G4TouchableHandle(fNewTouchablePtr); 
          fSplitPostStepPoint->SetTouchableHandle( fNewTouchableH );
        } else { 
          fSplitStep->GetPostStepPoint()->SetStepStatus( fullStepStatus );
          fSplitPostStepPoint->SetTouchableHandle( fFinalTouchableH );
        }


        // As first approximation, split the NIEL in the same fractions as the energy deposit
        G4double eLossFraction; 
        eLossFraction= (totalEnergyDeposit>0.0) ? energyLoss / totalEnergyDeposit : 1.0 ; 
        fSplitStep->SetNonIonizingEnergyDeposit(step.GetNonIonizingEnergyDeposit()*eLossFraction);
        
        fSplitPostStepPoint->SetSensitiveDetector( ptrSD ); 

        // Call the Sensitive Detector
        ptrSD->Hit(fSplitStep);

        if (verboseLevel>1) Verbose(step);
     }
  }

  // This must change the Stepping Control
  return pParticleChange;
}

G4TouchableHistory*
G4ScoreSplittingProcess::CreateTouchableForSubStep( G4int newVoxelNum, G4ThreeVector )
{
  // G4cout << " Creating touchable handle for voxel-no " << newVoxelNum << G4endl;

  G4TouchableHistory*  oldTouchableHistory= dynamic_cast<G4TouchableHistory*>(fOldTouchableH());
  G4TouchableHistory*  ptrTouchableHistory= G4TransportationManager::GetTransportationManager()->
                                            GetNavigatorForTracking()->CreateTouchableHistory(oldTouchableHistory->GetHistory());
  
  // Change the history
  G4NavigationHistory* ptrNavHistory= const_cast<G4NavigationHistory*>(ptrTouchableHistory->GetHistory());
  G4VPhysicalVolume*   curPhysicalVol= ptrNavHistory->GetTopVolume();

  EVolume curVolumeType=  ptrNavHistory->GetTopVolumeType();
  if( curVolumeType == kParameterised )
  { 
    ptrNavHistory->BackLevel(); 
    // G4VPVParameterised parameterisedPV= pNewMother
    G4VPVParameterisation* curParamstn=  curPhysicalVol->GetParameterisation();

    // From G4ParameterisedNavigation::IdentifyAndPlaceSolid() inline method
    G4VSolid* sampleSolid = curParamstn->ComputeSolid(newVoxelNum, curPhysicalVol);
    sampleSolid->ComputeDimensions(curParamstn, newVoxelNum, curPhysicalVol);
    curParamstn->ComputeTransformation(newVoxelNum, curPhysicalVol);

    ptrNavHistory->NewLevel( curPhysicalVol, kParameterised, newVoxelNum );
  }
  else
  {
    G4cout << " Current volume type is not Parameterised. " << G4endl; 
    G4Exception("G4ScoreSplittingProcess::CreateTouchableForSubStep",
          "ErrorRegularParamaterisation", JustWarning,
         "Score Splitting Process is used for Regular Structure - but did not find one here.");
  }
  return ptrTouchableHistory; 
}

void G4ScoreSplittingProcess::CopyStepStart(const G4Step & step)
{
  fSplitStep->SetTrack(step.GetTrack());
  fSplitStep->SetStepLength(step.GetStepLength());
  fSplitStep->SetTotalEnergyDeposit(step.GetTotalEnergyDeposit());
  fSplitStep->SetNonIonizingEnergyDeposit(step.GetNonIonizingEnergyDeposit());
  fSplitStep->SetControlFlag(step.GetControlFlag());

  *fSplitPreStepPoint  = *(step.GetPreStepPoint());

  fInitialTouchableH=  (step.GetPreStepPoint()) ->GetTouchableHandle();
  fFinalTouchableH =   (step.GetPostStepPoint())->GetTouchableHandle();
}

void G4ScoreSplittingProcess::Verbose(const G4Step& step) const
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
  G4cout << " StepLength : " << fSplitStep->GetStepLength()/mm
         << "      TotalEnergyDeposit : "
         << fSplitStep->GetTotalEnergyDeposit()/MeV << G4endl;
  G4cout << " PreStepPoint : "
         << fSplitStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() << " ["
           << fSplitStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber()
           << " ]" << " - ";
  if(fSplitStep->GetPreStepPoint()->GetProcessDefinedStep())
  { G4cout << fSplitStep->GetPreStepPoint()->GetProcessDefinedStep()->GetProcessName(); }
  else
  { G4cout << "NoProcessAssigned"; }
  G4cout << G4endl;
  G4cout << "                " << fSplitStep->GetPreStepPoint()->GetPosition() << G4endl;
  G4cout << " PostStepPoint : ";
  if(fSplitStep->GetPostStepPoint()->GetPhysicalVolume()) 
  {
    G4cout << fSplitStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() << " ["
           << fSplitStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber()
           << " ]";
  }
  else
  { G4cout << "OutOfWorld"; }
  G4cout << " - ";
  if(fSplitStep->GetPostStepPoint()->GetProcessDefinedStep())
  { G4cout << fSplitStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName(); }
  else
  { G4cout << "NoProcessAssigned"; }
  G4cout << G4endl;
  G4cout << "                 " << fSplitStep->GetPostStepPoint()->GetPosition() << " == "
         << fSplitStep->GetTrack()->GetMomentumDirection() 
         << G4endl;

}


//----------------------------------------------------------
//
//  AtRestGetPhysicalInteractionLength()
//
//----------------------------------------------------------
G4double 
G4ScoreSplittingProcess::AtRestGetPhysicalInteractionLength(
         const G4Track& /*track*/, 
         G4ForceCondition* condition)
{
  *condition = NotForced;  // Was Forced
  return DBL_MAX;
}


//---------------------------------------
//  AlongStepGetPhysicalInteractionLength
//---------------------------------------
G4double G4ScoreSplittingProcess::AlongStepGetPhysicalInteractionLength(
            const G4Track&   , // track, 
            G4double         , // previousStepSize, 
            G4double         , // currentMinimumStep,
            G4double&        , // proposedSafety, 
            G4GPILSelection* selection)
{
  *selection = NotCandidateForSelection;
  return DBL_MAX;
}

//------------------------------------
//             AlongStepDoIt()
//------------------------------------

G4VParticleChange* G4ScoreSplittingProcess::AlongStepDoIt(
    const G4Track& track, const G4Step& )
{
  // Dummy ParticleChange ie: does nothing
  // Expecting G4Transportation to move the track
  dummyParticleChange.Initialize(track);
  return &dummyParticleChange;
}

//------------------------------------
//             AtRestDoIt()
//------------------------------------
G4VParticleChange* G4ScoreSplittingProcess::AtRestDoIt(
     const G4Track& track,
     const G4Step&)
{ 
  pParticleChange->Initialize(track);
  return pParticleChange;
}
