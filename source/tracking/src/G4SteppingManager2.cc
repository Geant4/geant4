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
// $Id: G4SteppingManager2.cc 106991 2017-10-31 10:13:42Z gcosmo $
//
//---------------------------------------------------------------
//
// G4SteppingManager2.cc
//
// Description:
//   This class represents the manager who steers to move the give
//   particle from the TrackingManger by one Step.
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

//#define debug

#include "G4UImanager.hh"
#include "G4ForceCondition.hh"
#include "G4GPILSelection.hh"
#include "G4SteppingControl.hh"
#include "G4TransportationManager.hh"
#include "G4SteppingManager.hh"
#include "G4LossTableManager.hh"
#include "G4ParticleTable.hh"

/////////////////////////////////////////////////
void G4SteppingManager::GetProcessNumber()
/////////////////////////////////////////////////
{
#ifdef debug
  G4cout<<"G4SteppingManager::GetProcessNumber: is called track="
        <<fTrack<<G4endl;
#endif

  G4ProcessManager* pm= fTrack->GetDefinition()->GetProcessManager();
  if(!pm)
  {
    G4cerr << "ERROR - G4SteppingManager::GetProcessNumber()" << G4endl
           << "        ProcessManager is NULL for particle = "
           << fTrack->GetDefinition()->GetParticleName() << ", PDG_code = "
           << fTrack->GetDefinition()->GetPDGEncoding() << G4endl;
    G4Exception("G4SteppingManager::GetProcessNumber()", "Tracking0011",
                FatalException, "Process Manager is not found.");
    return;
  }

// AtRestDoits
   MAXofAtRestLoops =        pm->GetAtRestProcessVector()->entries();
   fAtRestDoItVector =       pm->GetAtRestProcessVector(typeDoIt);
   fAtRestGetPhysIntVector = pm->GetAtRestProcessVector(typeGPIL);
#ifdef debug
   G4cout << "G4SteppingManager::GetProcessNumber: #ofAtRest="
          << MAXofAtRestLoops << G4endl;
#endif

// AlongStepDoits
   MAXofAlongStepLoops = pm->GetAlongStepProcessVector()->entries();
   fAlongStepDoItVector = pm->GetAlongStepProcessVector(typeDoIt);
   fAlongStepGetPhysIntVector = pm->GetAlongStepProcessVector(typeGPIL);
#ifdef debug
   G4cout << "G4SteppingManager::GetProcessNumber:#ofAlongStp="
          << MAXofAlongStepLoops << G4endl;
#endif

// PostStepDoits
   MAXofPostStepLoops = pm->GetPostStepProcessVector()->entries();
   fPostStepDoItVector = pm->GetPostStepProcessVector(typeDoIt);
   fPostStepGetPhysIntVector = pm->GetPostStepProcessVector(typeGPIL);
#ifdef debug
   G4cout << "G4SteppingManager::GetProcessNumber: #ofPostStep="
          << MAXofPostStepLoops << G4endl;
#endif

   if (SizeOfSelectedDoItVector<MAXofAtRestLoops    ||
       SizeOfSelectedDoItVector<MAXofAlongStepLoops ||
       SizeOfSelectedDoItVector<MAXofPostStepLoops  )
   {
     G4cerr << "ERROR - G4SteppingManager::GetProcessNumber()" << G4endl
            << "        SizeOfSelectedDoItVector= " << SizeOfSelectedDoItVector
            << " ; is smaller then one of MAXofAtRestLoops= "
            << MAXofAtRestLoops << G4endl
            << "        or MAXofAlongStepLoops= " << MAXofAlongStepLoops
            << " or MAXofPostStepLoops= " << MAXofPostStepLoops << G4endl;
     G4Exception("G4SteppingManager::GetProcessNumber()",
                 "Tracking0012", FatalException,
                 "The array size is smaller than the actual No of processes.");
   }
}


// ************************************************************************
//
//  Private Member Functions
//
// ************************************************************************


/////////////////////////////////////////////////////////
 void G4SteppingManager::DefinePhysicalStepLength()
/////////////////////////////////////////////////////////
{

// ReSet the counter etc.
   PhysicalStep  = DBL_MAX;          // Initialize by a huge number    
   physIntLength = DBL_MAX;          // Initialize by a huge number    
#ifdef G4VERBOSE
                         // !!!!! Verbose
           if(verboseLevel>0) fVerbose->DPSLStarted();
#endif

// Obtain the user defined maximum allowed Step in the volume
//   1997.12.13 adds argument for  GetMaxAllowedStep by K.Kurashige
//   2004.01.20 This block will be removed by Geant4 7.0 
//   G4UserLimits* ul= fCurrentVolume->GetLogicalVolume()->GetUserLimits();
//   if (ul) {
//      physIntLength = ul->GetMaxAllowedStep(*fTrack);
//#ifdef G4VERBOSE
//                         // !!!!! Verbose
//           if(verboseLevel>0) fVerbose->DPSLUserLimit();
//#endif
//   }
//
//   if(physIntLength < PhysicalStep ){
//      PhysicalStep = physIntLength;
//      fStepStatus = fUserDefinedLimit;
//      fStep->GetPostStepPoint()
//           ->SetProcessDefinedStep(0);
//      // Take note that the process pointer is 'NULL' if the Step
//      // is defined by the user defined limit.
//   }
//   2004.01.20 This block will be removed by Geant4 7.0 

// GPIL for PostStep
   fPostStepDoItProcTriggered = MAXofPostStepLoops;

   for(size_t np=0; np < MAXofPostStepLoops; np++){
     fCurrentProcess = (*fPostStepGetPhysIntVector)(np);
     if (fCurrentProcess== 0) {
       (*fSelectedPostStepDoItVector)[np] = InActivated;
       continue;
     }   // NULL means the process is inactivated by a user on fly.

     physIntLength = fCurrentProcess->
                     PostStepGPIL( *fTrack,
                                                 fPreviousStepSize,
                                                      &fCondition );
#ifdef G4VERBOSE
                         // !!!!! Verbose
     if(verboseLevel>0) fVerbose->DPSLPostStep();
#endif

     switch (fCondition) {
	 case ExclusivelyForced:
	     (*fSelectedPostStepDoItVector)[np] = ExclusivelyForced;
	     fStepStatus = fExclusivelyForcedProc;
	     fStep->GetPostStepPoint()
		 ->SetProcessDefinedStep(fCurrentProcess);
	     break;
	 case Conditionally:
	   //	     (*fSelectedPostStepDoItVector)[np] = Conditionally;
	     G4Exception("G4SteppingManager::DefinePhysicalStepLength()", "Tracking1001", FatalException, "This feature no more supported");

	     break;
	 case Forced:
	     (*fSelectedPostStepDoItVector)[np] = Forced;
	     break;
	 case StronglyForced:
	     (*fSelectedPostStepDoItVector)[np] = StronglyForced;
	     break;
	 default:
	     (*fSelectedPostStepDoItVector)[np] = InActivated;
	     break;
     }
       


     if (fCondition==ExclusivelyForced) { 
 	 for(size_t nrest=np+1; nrest < MAXofPostStepLoops; nrest++){ 
 	     (*fSelectedPostStepDoItVector)[nrest] = InActivated; 
 	 } 
	 return;  // Take note the 'return' at here !!! 
     } 
     else{
	 if(physIntLength < PhysicalStep ){
	     PhysicalStep = physIntLength;
	     fStepStatus = fPostStepDoItProc;
	     fPostStepDoItProcTriggered = G4int(np);
	     fStep->GetPostStepPoint()
		 ->SetProcessDefinedStep(fCurrentProcess);
	 }
     }
     

   }

   if (fPostStepDoItProcTriggered<MAXofPostStepLoops) {
       if ((*fSelectedPostStepDoItVector)[fPostStepDoItProcTriggered] == 
	   InActivated) {
	   (*fSelectedPostStepDoItVector)[fPostStepDoItProcTriggered] = 
	       NotForced;
       }
   }

// GPIL for AlongStep
   proposedSafety = DBL_MAX;
   G4double safetyProposedToAndByProcess = proposedSafety;

   for(size_t kp=0; kp < MAXofAlongStepLoops; kp++){
     fCurrentProcess = (*fAlongStepGetPhysIntVector)[kp];
     if (fCurrentProcess== 0) continue;
         // NULL means the process is inactivated by a user on fly.

     physIntLength = fCurrentProcess->
                     AlongStepGPIL( *fTrack, fPreviousStepSize,
                                     PhysicalStep,
				     safetyProposedToAndByProcess,
                                    &fGPILSelection );
#ifdef G4VERBOSE
                         // !!!!! Verbose
     if(verboseLevel>0) fVerbose->DPSLAlongStep();
#endif
     if(physIntLength < PhysicalStep){
       PhysicalStep = physIntLength;

       // Check if the process wants to be the GPIL winner. For example,
       // multi-scattering proposes Step limit, but won't be the winner.
       if(fGPILSelection==CandidateForSelection){
          fStepStatus = fAlongStepDoItProc;
          fStep->GetPostStepPoint()
               ->SetProcessDefinedStep(fCurrentProcess);
       }

    	  // Transportation is assumed to be the last process in the vector
       if(kp == MAXofAlongStepLoops-1)
	  fStepStatus = fGeomBoundary;
     }

     // Make sure to check the safety, even if Step is not limited 
     //  by this process.                      J. Apostolakis, June 20, 1998
     // 
     if (safetyProposedToAndByProcess < proposedSafety)
        // proposedSafety keeps the smallest value:
        proposedSafety               = safetyProposedToAndByProcess;
     else
        // safetyProposedToAndByProcess always proposes a valid safety:
        safetyProposedToAndByProcess = proposedSafety;
      
   } 
} // void G4SteppingManager::DefinePhysicalStepLength() //


//////////////////////////////////////////////////////
void G4SteppingManager::InvokeAtRestDoItProcs()
//////////////////////////////////////////////////////
{
// Select the rest process which has the shortest time before
// it is invoked. In rest processes, GPIL()
// returns the time before a process occurs.
   G4double lifeTime, shortestLifeTime;

   fAtRestDoItProcTriggered = 0;
   shortestLifeTime = DBL_MAX;

   unsigned int NofInactiveProc=0;
   for( size_t ri=0 ; ri < MAXofAtRestLoops ; ri++ ){
     fCurrentProcess = (*fAtRestGetPhysIntVector)[ri];
     if (fCurrentProcess== 0) {
       (*fSelectedAtRestDoItVector)[ri] = InActivated;
       NofInactiveProc++;
       continue;
     }   // NULL means the process is inactivated by a user on fly.

     lifeTime =
       fCurrentProcess->AtRestGPIL( *fTrack, &fCondition );

     

     if(fCondition==Forced && fCurrentProcess){
       (*fSelectedAtRestDoItVector)[ri] = Forced;
     }
     else{
       (*fSelectedAtRestDoItVector)[ri] = InActivated;
       if(lifeTime < shortestLifeTime ){
          shortestLifeTime = lifeTime;
          fAtRestDoItProcTriggered =  G4int(ri);
       }
     }
   }

  (*fSelectedAtRestDoItVector)[fAtRestDoItProcTriggered] = NotForced;

// at least one process is necessary to destroy the particle  
// exit with warning 
   if(NofInactiveProc==MAXofAtRestLoops){ 
     G4Exception("G4SteppingManager::InvokeAtRestDoItProcs()", "Tracking0013",
                 JustWarning, "No AtRestDoIt process is active!" );
   }

   fStep->SetStepLength( 0. );  //the particle has stopped
   fTrack->SetStepLength( 0. );

// invoke selected process
   for(size_t np=0; np < MAXofAtRestLoops; np++){
   //
   // Note: DoItVector has inverse order against GetPhysIntVector
   //       and SelectedAtRestDoItVector.
   //
     if( (*fSelectedAtRestDoItVector)[MAXofAtRestLoops-np-1] != InActivated){

       fCurrentProcess = (*fAtRestDoItVector)[np];
       fParticleChange 
         = fCurrentProcess->AtRestDoIt( *fTrack, *fStep);
                               
       // Set the current process as a process which defined this Step length
       fStep->GetPostStepPoint()
              ->SetProcessDefinedStep(fCurrentProcess);

       // Update Step  
       fParticleChange->UpdateStepForAtRest(fStep);

       // Now Store the secondaries from ParticleChange to SecondaryList
       G4Track* tempSecondaryTrack;
       G4int    num2ndaries;

       num2ndaries = fParticleChange->GetNumberOfSecondaries();

       for(G4int DSecLoop=0 ; DSecLoop< num2ndaries; DSecLoop++){
         tempSecondaryTrack = fParticleChange->GetSecondary(DSecLoop);

         if(tempSecondaryTrack->GetDefinition()->GetApplyCutsFlag())
         { ApplyProductionCut(tempSecondaryTrack); }

         // Set parentID 
         tempSecondaryTrack->SetParentID( fTrack->GetTrackID() );

	 // Set the process pointer which created this track 
	 tempSecondaryTrack->SetCreatorProcess( fCurrentProcess );
	 
	 // If this 2ndry particle has 'zero' kinetic energy, make sure
	 // it invokes a rest process at the beginning of the tracking
	 if(tempSecondaryTrack->GetKineticEnergy() <= DBL_MIN){
	   G4ProcessManager* pm = tempSecondaryTrack->GetDefinition()->GetProcessManager();
           if(!pm && tempSecondaryTrack->GetDefinition()->IsGeneralIon())
             { pm = G4ParticleTable::GetParticleTable()->GetGenericIon()->GetProcessManager(); }
#ifdef G4MUATOMS_INUSE
           if(!pm && tempSecondaryTrack->GetDefinition()->IsMuonicAtom())
             { pm = G4ParticleTable::GetParticleTable()->GetGenericMuonicAtom()->GetProcessManager(); }
#endif
	   if (pm->GetAtRestProcessVector()->entries()>0){
	     tempSecondaryTrack->SetTrackStatus( fStopButAlive );
	     fSecondary->push_back( tempSecondaryTrack );
             fN2ndariesAtRestDoIt++;
	   } else {
	     delete tempSecondaryTrack;
	   }
	 } else {
	   fSecondary->push_back( tempSecondaryTrack );
           fN2ndariesAtRestDoIt++;
	 }	
       } //end of loop on secondary 


       // clear ParticleChange
       fParticleChange->Clear();

     } //if(fSelectedAtRestDoItVector[np] != InActivated){
   } //for(size_t np=0; np < MAXofAtRestLoops; np++){
   fStep->UpdateTrack();

   fTrack->SetTrackStatus( fStopAndKill );
}


/////////////////////////////////////////////////////////
void G4SteppingManager::InvokeAlongStepDoItProcs()
/////////////////////////////////////////////////////////
{

// If the current Step is defined by a 'ExclusivelyForced' 
// PostStepDoIt, then don't invoke any AlongStepDoIt
   if(fStepStatus == fExclusivelyForcedProc){
     return;               // Take note 'return' at here !!!
   }

// Invoke the all active continuous processes
   for( size_t ci=0 ; ci<MAXofAlongStepLoops ; ci++ ){
     fCurrentProcess = (*fAlongStepDoItVector)[ci];
     if (fCurrentProcess== 0) continue;
         // NULL means the process is inactivated by a user on fly.

     fParticleChange 
       = fCurrentProcess->AlongStepDoIt( *fTrack, *fStep );

     // Update the PostStepPoint of Step according to ParticleChange
     fParticleChange->UpdateStepForAlongStep(fStep);
#ifdef G4VERBOSE
                         // !!!!! Verbose
               if(verboseLevel>0) fVerbose->AlongStepDoItOneByOne();
#endif

     // Now Store the secondaries from ParticleChange to SecondaryList
     G4Track* tempSecondaryTrack;
     G4int    num2ndaries;

     num2ndaries = fParticleChange->GetNumberOfSecondaries();

     for(G4int DSecLoop=0 ; DSecLoop< num2ndaries; DSecLoop++){
         tempSecondaryTrack = fParticleChange->GetSecondary(DSecLoop);

         if(tempSecondaryTrack->GetDefinition()->GetApplyCutsFlag())
         { ApplyProductionCut(tempSecondaryTrack); }

         // Set parentID
         tempSecondaryTrack->SetParentID( fTrack->GetTrackID() );

	 // Set the process pointer which created this track 
	 tempSecondaryTrack->SetCreatorProcess( fCurrentProcess );

	 // If this 2ndry particle has 'zero' kinetic energy, make sure
	 // it invokes a rest process at the beginning of the tracking
	 if(tempSecondaryTrack->GetKineticEnergy() <= DBL_MIN){
	   G4ProcessManager* pm = tempSecondaryTrack->GetDefinition()->GetProcessManager();
           if(!pm && tempSecondaryTrack->GetDefinition()->IsGeneralIon())
             { pm = G4ParticleTable::GetParticleTable()->GetGenericIon()->GetProcessManager(); }
#ifdef G4MUATOMS_INUSE
           if(!pm && tempSecondaryTrack->GetDefinition()->IsMuonicAtom())
             { pm = G4ParticleTable::GetParticleTable()->GetGenericMuonicAtom()->GetProcessManager(); }
#endif
	   if (pm->GetAtRestProcessVector()->entries()>0){
	     tempSecondaryTrack->SetTrackStatus( fStopButAlive );
	     fSecondary->push_back( tempSecondaryTrack );
             fN2ndariesAlongStepDoIt++;
	   } else {
	     delete tempSecondaryTrack;
	   }
	 } else {
	   fSecondary->push_back( tempSecondaryTrack );
           fN2ndariesAlongStepDoIt++;
	 }
     } //end of loop on secondary 
     
     // Set the track status according to what the process defined
     // if kinetic energy >0, otherwise set  fStopButAlive
     fTrack->SetTrackStatus( fParticleChange->GetTrackStatus() );
     
     // clear ParticleChange
     fParticleChange->Clear();
   }

   fStep->UpdateTrack();
   G4TrackStatus fNewStatus = fTrack->GetTrackStatus();

   if ( fNewStatus == fAlive && fTrack->GetKineticEnergy() <= DBL_MIN ) {
     if(MAXofAtRestLoops>0) fNewStatus = fStopButAlive;
     else                   fNewStatus = fStopAndKill;
     fTrack->SetTrackStatus( fNewStatus );
   }

}

////////////////////////////////////////////////////////
void G4SteppingManager::InvokePostStepDoItProcs()
////////////////////////////////////////////////////////
{

// Invoke the specified discrete processes
   for(size_t np=0; np < MAXofPostStepLoops; np++){
   //
   // Note: DoItVector has inverse order against GetPhysIntVector
   //       and SelectedPostStepDoItVector.
   //
     G4int Cond = (*fSelectedPostStepDoItVector)[MAXofPostStepLoops-np-1];
     if(Cond != InActivated){
       if( ((Cond == NotForced) && (fStepStatus == fPostStepDoItProc)) ||
	   ((Cond == Forced) && (fStepStatus != fExclusivelyForcedProc)) ||
	   //	   ((Cond == Conditionally) && (fStepStatus == fAlongStepDoItProc)) ||
	   ((Cond == ExclusivelyForced) && (fStepStatus == fExclusivelyForcedProc)) || 
	   ((Cond == StronglyForced) ) 
	  ) {

	 InvokePSDIP(np);
         if ((np==0) && (fTrack->GetNextVolume() == 0)){
           fStepStatus = fWorldBoundary;
           fStep->GetPostStepPoint()->SetStepStatus( fStepStatus );
         }
       }
     } //if(*fSelectedPostStepDoItVector(np)........

     // Exit from PostStepLoop if the track has been killed,
     // but extra treatment for processes with Strongly Forced flag
     if(fTrack->GetTrackStatus() == fStopAndKill) {
       for(size_t np1=np+1; np1 < MAXofPostStepLoops; np1++){ 
	 G4int Cond2 = (*fSelectedPostStepDoItVector)[MAXofPostStepLoops-np1-1];
	 if (Cond2 == StronglyForced) {
	   InvokePSDIP(np1);
         }
       }
       break;
     }
   } //for(size_t np=0; np < MAXofPostStepLoops; np++){
}



void G4SteppingManager::InvokePSDIP(size_t np)
{
         fCurrentProcess = (*fPostStepDoItVector)[np];
         fParticleChange 
            = fCurrentProcess->PostStepDoIt( *fTrack, *fStep);

         // Update PostStepPoint of Step according to ParticleChange
	 fParticleChange->UpdateStepForPostStep(fStep);
#ifdef G4VERBOSE
                 // !!!!! Verbose
           if(verboseLevel>0) fVerbose->PostStepDoItOneByOne();
#endif
         // Update G4Track according to ParticleChange after each PostStepDoIt
         fStep->UpdateTrack();

         // Update safety after each invocation of PostStepDoIts
         fStep->GetPostStepPoint()->SetSafety( CalculateSafety() );

         // Now Store the secondaries from ParticleChange to SecondaryList
         G4Track* tempSecondaryTrack;
         G4int    num2ndaries;

         num2ndaries = fParticleChange->GetNumberOfSecondaries();

         for(G4int DSecLoop=0 ; DSecLoop< num2ndaries; DSecLoop++){
            tempSecondaryTrack = fParticleChange->GetSecondary(DSecLoop);
   
            if(tempSecondaryTrack->GetDefinition()->GetApplyCutsFlag())
            { ApplyProductionCut(tempSecondaryTrack); }

            // Set parentID 
            tempSecondaryTrack->SetParentID( fTrack->GetTrackID() );
	    
	    // Set the process pointer which created this track 
	    tempSecondaryTrack->SetCreatorProcess( fCurrentProcess );

            // If this 2ndry particle has 'zero' kinetic energy, make sure
            // it invokes a rest process at the beginning of the tracking
	    if(tempSecondaryTrack->GetKineticEnergy() <= DBL_MIN){
	      G4ProcessManager* pm = tempSecondaryTrack->GetDefinition()->GetProcessManager();
              if(!pm && tempSecondaryTrack->GetDefinition()->IsGeneralIon())
                { pm = G4ParticleTable::GetParticleTable()->GetGenericIon()->GetProcessManager(); }
#ifdef G4MUATOMS_INUSE
              if(!pm && tempSecondaryTrack->GetDefinition()->IsMuonicAtom())
                { pm = G4ParticleTable::GetParticleTable()->GetGenericMuonicAtom()->GetProcessManager(); }
#endif
	      if (pm->GetAtRestProcessVector()->entries()>0){
		tempSecondaryTrack->SetTrackStatus( fStopButAlive );
		fSecondary->push_back( tempSecondaryTrack );
                fN2ndariesPostStepDoIt++;
 	      } else {
		delete tempSecondaryTrack;
	      }
	    } else {
	      fSecondary->push_back( tempSecondaryTrack );
              fN2ndariesPostStepDoIt++;
	    }
         } //end of loop on secondary 

         // Set the track status according to what the process defined
         fTrack->SetTrackStatus( fParticleChange->GetTrackStatus() );

         // clear ParticleChange
         fParticleChange->Clear();
}

#include "G4EnergyLossTables.hh"
#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"

void G4SteppingManager::ApplyProductionCut(G4Track* aSecondary)
{
  G4bool tBelowCutEnergyAndSafety = false;
  G4int tPtclIdx
    = G4ProductionCuts::GetIndex(aSecondary->GetDefinition());
  if (tPtclIdx<0)  { return; }
  G4ProductionCutsTable* tCutsTbl
    = G4ProductionCutsTable::GetProductionCutsTable();
  G4int tCoupleIdx
    = tCutsTbl->GetCoupleIndex(fPreStepPoint->GetMaterialCutsCouple());
  G4double tProdThreshold
    = (*(tCutsTbl->GetEnergyCutsVector(tPtclIdx)))[tCoupleIdx];
  if( aSecondary->GetKineticEnergy()<tProdThreshold )
  {
    tBelowCutEnergyAndSafety = true;
    if(std::abs(aSecondary->GetDynamicParticle()->GetCharge()) > DBL_MIN)
    {
      G4double currentRange
        = G4LossTableManager::Instance()->GetRange(aSecondary->GetDefinition(),
                 aSecondary->GetKineticEnergy(),
                 fPreStepPoint->GetMaterialCutsCouple());
      tBelowCutEnergyAndSafety = (currentRange < CalculateSafety() );
    }
  }

  if( tBelowCutEnergyAndSafety )
  {
    if( !(aSecondary->IsGoodForTracking()) )
    {
//G4cout << "!! Warning - G4SteppingManager:" << G4endl
//<< " This physics process generated a secondary"
//<< " of which energy is below cut but"
//<< " GoodForTracking is off !!!!!" << G4endl;
    // Add kinetic energy to the total energy deposit
      fStep->AddTotalEnergyDeposit(
      aSecondary->GetKineticEnergy() );
      aSecondary->SetKineticEnergy(0.0);
    }
  }
}

