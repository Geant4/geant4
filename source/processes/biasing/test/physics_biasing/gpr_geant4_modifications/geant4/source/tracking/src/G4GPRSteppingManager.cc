#include "G4GPRSteppingManager.hh"
#include "G4GPRManager.hh"
#include "G4GPRTriggering.hh"
#include "G4VSensitiveDetector.hh"

// try to just switch in gpr stuff & not change processing logic too much
G4StepStatus G4GPRSteppingManager::Stepping()
{
// Store last PostStepPoint to PreStepPoint, and swap current and nex
// volume information of G4Track. Reset total energy deposit in one Step. 
   fStep->CopyPostToPreStepPoint();
   fStep->ResetTotalEnergyDeposit();

// Switch next touchable in track to current one
   fTrack->SetTouchableHandle(fTrack->GetNextTouchableHandle());

// Reset the secondary particles
   fN2ndariesAtRestDoIt = 0;
   fN2ndariesAlongStepDoIt = 0;
   fN2ndariesPostStepDoIt = 0;

//JA Set the volume before it is used (in DefineStepLength() for User Limit) 
   fCurrentVolume = fStep->GetPreStepPoint()->GetPhysicalVolume();

// Reset the step's auxiliary points vector pointer
   fStep->SetPointerToVectorOfAuxiliaryPoints(0);

//-----------------
// AtRest Processes
//-----------------
   G4GPRManager* gprManager = fTrack->GetDefinition()->GetGPRManager();
   gprManager->Fire<G4GPRTriggering::Stepping::StartStep>(*fTrack, *fStep);

   gprManager->GetList<G4GPRProcessLists::AtRestGPIL>(pAtRestGPIL);
   MAXofAtRestLoops = pAtRestGPIL->size();

   if( fTrack->GetTrackStatus() == fStopButAlive ){
     if( MAXofAtRestLoops>0 ){
        InvokeAtRestDoItProcs();
        fStepStatus = fAtRestDoItProc;
        fStep->GetPostStepPoint()->SetStepStatus( fStepStatus );
       
#ifdef G4VERBOSE
            // !!!!! Verbose
             if(verboseLevel>0) fVerbose->AtRestDoItInvoked();
#endif 

     }
     // Make sure the track is killed
     fTrack->SetTrackStatus( fStopAndKill );
   }

//---------------------------------
// AlongStep and PostStep Processes
//---------------------------------


   else{
     // Find minimum Step length demanded by active disc./cont. processes
     DefinePhysicalStepLength();
     G4cout<<"jane gpil "<<fTrack->GetDefinition()->GetParticleName()<<" "<<fTrack->GetTrackID()<<" "<<physIntLength<<" "<<PhysicalStep<<G4endl;
     // Store the Step length (geometrical length) to G4Step and G4Track
     fStep->SetStepLength( PhysicalStep );
     fTrack->SetStepLength( PhysicalStep );
     G4double GeomStepLength = PhysicalStep;

     // Store StepStatus to PostStepPoint
     fStep->GetPostStepPoint()->SetStepStatus( fStepStatus );

     // Invoke AlongStepDoIt 
     InvokeAlongStepDoItProcs();

     // Update track by taking into account all changes by AlongStepDoIt
     fStep->UpdateTrack();

     // Update safety after invocation of all AlongStepDoIts
     endpointSafOrigin= fPostStepPoint->GetPosition();
     endpointSafety=  std::max( proposedSafety - GeomStepLength, 0.);

     fStep->GetPostStepPoint()->SetSafety( endpointSafety );

#ifdef G4VERBOSE
                         // !!!!! Verbose
           if(verboseLevel>0) fVerbose->AlongStepDoItAllDone();
#endif

     // Invoke PostStepDoIt
     InvokePostStepDoItProcs();

#ifdef G4VERBOSE
                 // !!!!! Verbose
     if(verboseLevel>0) fVerbose->PostStepDoItAllDone();
#endif
   }

//-------
// Finale
//-------

// Update 'TrackLength' and remeber the Step length of the current Step
   fTrack->AddTrackLength(fStep->GetStepLength());
   fPreviousStepSize = fStep->GetStepLength();
#ifdef G4VERBOSE
                         // !!!!! Verbose

           if(verboseLevel>0) fVerbose->StepInfo();
#endif
// Send G4Step information to Hit/Dig if the volume is sensitive
   fCurrentVolume = fStep->GetPreStepPoint()->GetPhysicalVolume();
   StepControlFlag =  fStep->GetControlFlag();
   if( fCurrentVolume != 0 && StepControlFlag != AvoidHitInvocation) {
      fSensitive = fStep->GetPreStepPoint()->
                                   GetSensitiveDetector();
      if( fSensitive != 0 ) {
        fSensitive->Hit(fStep);
      }
   }

// User intervention process.
   fStep->SetTrack(fTrack);
   if( fUserSteppingAction != NULL ) {
      fUserSteppingAction->UserSteppingAction(fStep);
   }

// Stepping process finish. Return the value of the StepStatus.
   return fStepStatus;

}

void G4GPRSteppingManager::DefinePhysicalStepLength()
{
  G4cout<<"jane G4GPRSteppingManager::DefinePhysicalStepLength()"<<G4endl;

// ReSet the counter etc.
   PhysicalStep  = DBL_MAX;          // Initialize by a huge number    
   physIntLength = DBL_MAX;          // Initialize by a huge number    

   G4GPRManager* gprManager = fTrack->GetDefinition()->GetGPRManager();
   gprManager->GetList<G4GPRProcessLists::DiscreteGPIL>(pDiscreteGPIL);

// GPIL for PostStep
   MAXofPostStepLoops = pDiscreteGPIL->size();
   fPostStepDoItProcTriggered = MAXofPostStepLoops;

   for(size_t np=0; np < MAXofPostStepLoops; np++){

     //     fCurrentProcess = (*fPostStepGetPhysIntVector)(np);
     //  pCurrentProcess = (*pDiscreteGPIL)(np);
     
     //  jane fixme somehow or figure out logic for loading inactive 
//       processes into gpr
  //   if (CurrentProcess== NULL) {
    //   (*fSelectedPostStepDoItVector)[np] = InActivated;
      // continue;
     //}   // NULL means the process is inactivated by a user on fly.

     physIntLength = (*pDiscreteGPIL)[np].operator()( *fTrack,
						      fPreviousStepSize,
						      &fCondition );

     switch (fCondition) {
     case ExclusivelyForced:
       (*fSelectedPostStepDoItVector)[np] = ExclusivelyForced;
       
       fStepStatus = fExclusivelyForcedProc;

       //	 jane fixme somehow 
       //	   fStep->GetPostStepPoint()
       //		 ->SetProcessDefinedStep(fCurrentProcess);
       
       break;
     case Conditionally:
       (*fSelectedPostStepDoItVector)[np] = Conditionally;
       break;
     case Forced:
       (*fSelectedPostStepDoItVector)[np] = Forced;
       break;
     case StronglyForced:
       (*fSelectedPostStepDoItVector)[np] = StronglyForced;
       break;
     default:
       (*fSelectedPostStepDoItVector)[np] = InActivated;
       if(physIntLength < PhysicalStep ){
         PhysicalStep = physIntLength;
	 fStepStatus = fPostStepDoItProc;
	 fPostStepDoItProcTriggered = G4int(np);
	 
	 //  jane fixme somehow
         //fStep->GetPostStepPoint()
           //   ->SetProcessDefinedStep(fCurrentProcess);
	 
       }
     }
     if (fCondition==ExclusivelyForced) { 
	 for(size_t nrest=np+1; nrest < MAXofPostStepLoops; nrest++){ 
	     (*fSelectedPostStepDoItVector)[nrest] = InActivated; 
	 } 
	 return;  // Take note the 'return' at here !!! 
     } 
   }

   if(fPostStepDoItProcTriggered<MAXofPostStepLoops)
     (*fSelectedPostStepDoItVector)[fPostStepDoItProcTriggered] = NotForced;


// GPIL for AlongStep
   proposedSafety = DBL_MAX;
   G4double safetyProposedToAndByProcess = proposedSafety;

   gprManager->GetList<G4GPRProcessLists::ContinuousGPIL>(pContinuousGPIL);

   MAXofAlongStepLoops = pContinuousGPIL->size();

   for(size_t kp=0; kp < MAXofAlongStepLoops; kp++){
     //     fCurrentProcess = (*fAlongStepGetPhysIntVector)[kp];
     // pCurrentProcess = (*fAlongStepGetPhysIntVector)[kp];
     //     if (fCurrentProcess== NULL) continue; jane fixme somehow
         // NULL means the process is inactivated by a user on fly.

     physIntLength = (*pContinuousGPIL)[kp].operator()(*fTrack,
						       fPreviousStepSize,
						       PhysicalStep,
						       safetyProposedToAndByProcess,
						       &fGPILSelection );

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
       if(kp == MAXofAlongStepLoops-1) {
	   if (fTrack->GetNextVolume() != 0)
	       fStepStatus = fGeomBoundary;
	   else
	       fStepStatus = fWorldBoundary;	
       }
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
}

//////////////////////////////////////////////////////
void G4GPRSteppingManager::InvokeAtRestDoItProcs()
//////////////////////////////////////////////////////
{

// Select the rest process which has the shortest time before
// it is invoked. In rest processes, GPIL()
// returns the time before a process occurs.
   G4double lifeTime, shortestLifeTime;

   fAtRestDoItProcTriggered = 0;
   shortestLifeTime = DBL_MAX;

   //   G4GPRManager* gprManager = fTrack->GetDefinition()->GetGPRManager();
   //   gprManager->GetList<G4GPRProcessLists::AtRestGPIL>(pAtRestGPIL);
   //   MAXofAtRestLoops = pAtRestGPIL->size();

   unsigned int NofInactiveProc=0;
   for( size_t ri=0 ; ri < MAXofAtRestLoops ; ri++ ){
     
     //     fCurrentProcess = (*fAtRestGetPhysIntVector)[ri];
     //jane fixme 
     //if (fCurrentProcess== NULL) {
     //  (*fSelectedAtRestDoItVector)[ri] = InActivated;
     // NofInactiveProc++;
     // continue;
     //}   // NULL means the process is inactivated by a user on fly.
     G4cout<<"jane rest gpil"<<G4endl;

     lifeTime = (*pAtRestGPIL)[ri].operator()(*fTrack,
					      &fCondition );
     G4cout<<"jane done rest gpil"<<G4endl;
     if(fCondition==Forced){
       (*fSelectedAtRestDoItVector)[ri] = Forced;
     }
     else{
       (*fSelectedAtRestDoItVector)[ri] = InActivated;
      if(lifeTime < shortestLifeTime ){
         shortestLifeTime = lifeTime;
	 fAtRestDoItProcTriggered =  G4int(int(ri)); 
       }
     }
   }

// at least one process is necessary to destory the particle  
// exit with warning 
   if(NofInactiveProc==MAXofAtRestLoops){ 
     //     G4Exception("G4SteppingManager::InvokeAtRestDoItProcs: No AtRestDoIt process is active. " );
     G4cerr << "G4SteppingManager::InvokeAtRestDoItProcs: No AtRestDoIt process is active. " << G4endl;
   } else {
        (*fSelectedAtRestDoItVector)[fAtRestDoItProcTriggered] = NotForced;
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
       G4cout<<"jane rest doit"<<G4endl;
       G4GPRManager* gprManager = fTrack->GetDefinition()->GetGPRManager();
       gprManager->GetList<G4GPRProcessLists::AtRestDoIt>(pAtRestDoIt);
       fParticleChange = (*pAtRestDoIt)[np](*fTrack, *fStep);
       G4cout<<"jane done rest doit"<<G4endl;
	 //fCurrentProcess = (*fAtRestDoItVector)[np];
	 //       fParticleChange 
	 //         = fCurrentProcess->AtRestDoIt( *fTrack, *fStep);
                               
       // Set the current process as a process which defined this Step length
       //jane fixme
       //fStep->GetPostStepPoint()
       //     ->SetProcessDefinedStep(fCurrentProcess);

       // Update Step  
       fParticleChange->UpdateStepForAtRest(fStep);

       // Now Store the secondaries from ParticleChange to SecondaryList
       G4Track* tempSecondaryTrack;
       G4int    num2ndaries;

       num2ndaries = fParticleChange->GetNumberOfSecondaries();
       
       fN2ndariesAtRestDoIt += num2ndaries;

       for(G4int DSecLoop=0 ; DSecLoop< num2ndaries; DSecLoop++){
         tempSecondaryTrack = fParticleChange->GetSecondary(DSecLoop);

         if(tempSecondaryTrack->GetDefinition()->GetApplyCutsFlag())
         { ApplyProductionCut(tempSecondaryTrack); }

         // Set parentID 
         tempSecondaryTrack->SetParentID( fTrack->GetTrackID() );

	 // Set the process pointer which created this track 
	 //jane fixme
	 //tempSecondaryTrack->SetCreatorProcess( fCurrentProcess );
	 
	 // If this 2ndry particle has 'zero' kinetic energy, make sure
	 // it invokes a rest process at the beginning of the tracking
	 // jane fixme, huh ?
	 if(tempSecondaryTrack->GetKineticEnergy() <= DBL_MIN){
	   G4ProcessManager* pm = tempSecondaryTrack->GetDefinition()->GetProcessManager();
	   if (pm->GetAtRestProcessVector()->entries()>0){
	     tempSecondaryTrack->SetTrackStatus( fStopButAlive );
	     fSecondary->push_back( tempSecondaryTrack );
	   } else {
	     delete tempSecondaryTrack;
	   }
	 } else {
	   fSecondary->push_back( tempSecondaryTrack );
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
void G4GPRSteppingManager::InvokeAlongStepDoItProcs()
/////////////////////////////////////////////////////////
{

// If the current Step is defined by a 'ExclusivelyForced' 
// PostStepDoIt, then don't invoke any AlongStepDoIt
   if(fStepStatus == fExclusivelyForcedProc){
     return;               // Take note 'return' at here !!!
   }

   G4GPRManager* gprManager = fTrack->GetDefinition()->GetGPRManager();
   gprManager->GetList<G4GPRProcessLists::ContinuousDoIt>(pContinuousDoIt);

// Invoke the all active continuous processes
   for( size_t ci=0 ; ci<MAXofAlongStepLoops ; ci++ ){
     //     fCurrentProcess = (*fAlongStepDoItVector)[ci];
     //jane fixmeif (fCurrentProcess== NULL) continue;
         // NULL means the process is inactivated by a user on fly.
     
     //     fParticleChange 
     //= fCurrentProcess->AlongStepDoIt( *fTrack, *fStep );
     G4cout<<"jane gpr along step "<<(*pContinuousDoIt)[ci].GetIdentifier()<<G4endl;
     fParticleChange = (*pContinuousDoIt)[ci].operator()( *fTrack, *fStep );
     //       = fCurrentProcess->AlongStepDoIt( *fTrack, *fStep );

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
     fN2ndariesAlongStepDoIt += num2ndaries;

     for(G4int DSecLoop=0 ; DSecLoop< num2ndaries; DSecLoop++){
         tempSecondaryTrack = fParticleChange->GetSecondary(DSecLoop);

         if(tempSecondaryTrack->GetDefinition()->GetApplyCutsFlag())
         { ApplyProductionCut(tempSecondaryTrack); }

         // Set parentID
         tempSecondaryTrack->SetParentID( fTrack->GetTrackID() );

	 //jane fixme
	 // Set the process pointer which created this track 
	 //	 tempSecondaryTrack->SetCreatorProcess( fCurrentProcess );

	 // If this 2ndry particle has 'zero' kinetic energy, make sure
	 // it invokes a rest process at the beginning of the tracking
	 if(tempSecondaryTrack->GetKineticEnergy() <= DBL_MIN){
	   //jane fixme
	   G4ProcessManager* pm = tempSecondaryTrack->GetDefinition()->GetProcessManager();
	   if (pm->GetAtRestProcessVector()->entries()>0){
	     tempSecondaryTrack->SetTrackStatus( fStopButAlive );
	     fSecondary->push_back( tempSecondaryTrack );
	   } else {
	     delete tempSecondaryTrack;
	   }
	 } else {
	   fSecondary->push_back( tempSecondaryTrack );
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
void G4GPRSteppingManager::InvokePostStepDoItProcs()
////////////////////////////////////////////////////////
{
  G4GPRManager* gprManager = fTrack->GetDefinition()->GetGPRManager();
  gprManager->GetList<G4GPRProcessLists::DiscreteDoIt>(pDiscreteDoIt);

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
	   ((Cond == Conditionally) && (fStepStatus == fAlongStepDoItProc)) ||
	   ((Cond == ExclusivelyForced) && (fStepStatus == fExclusivelyForcedProc)) || 
	   ((Cond == StronglyForced) ) 
	  ) {

	 InvokePSDIP(np);
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


void G4GPRSteppingManager::InvokePSDIP(size_t np)
{
  // fCurrentProcess = (*fPostStepDoItVector)[np];
  //     fParticleChange 
  //        = fCurrentProcess->PostStepDoIt( *fTrack, *fStep);

  //         fCurrentProcess = (*fPostStepDoItVector)[np];
  G4cout<<"jane gpr discrete "<<(*pDiscreteDoIt)[np].GetIdentifier()<<G4endl;
         fParticleChange = (*pDiscreteDoIt)[np].operator()( *fTrack, *fStep);
	 //            = fCurrentProcess->PostStepDoIt( *fTrack, *fStep);

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
	 G4cout<<"jane gpr npoststep secondaries "<<fTrack->GetDefinition()->GetParticleName()<<" "<<fTrack->GetTrackID()<<" "<<num2ndaries<<G4endl;

         fN2ndariesPostStepDoIt += num2ndaries;

         for(G4int DSecLoop=0 ; DSecLoop< num2ndaries; DSecLoop++){
            tempSecondaryTrack = fParticleChange->GetSecondary(DSecLoop);
   
            if(tempSecondaryTrack->GetDefinition()->GetApplyCutsFlag())
            { ApplyProductionCut(tempSecondaryTrack); }

            // Set parentID 
            tempSecondaryTrack->SetParentID( fTrack->GetTrackID() );
	    
	    // Set the process pointer which created this track 
	    //jane fixme	    tempSecondaryTrack->SetCreatorProcess( fCurrentProcess );

            // If this 2ndry particle has 'zero' kinetic energy, make sure
            // it invokes a rest process at the beginning of the tracking
	    if(tempSecondaryTrack->GetKineticEnergy() <= DBL_MIN){
	      G4ProcessManager* pm = tempSecondaryTrack->GetDefinition()->GetProcessManager();
	      //jane fixme
	      if (pm->GetAtRestProcessVector()->entries()>0){
		tempSecondaryTrack->SetTrackStatus( fStopButAlive );
		fSecondary->push_back( tempSecondaryTrack );
 	      } else {
		delete tempSecondaryTrack;
	      }
	    } else {
	      fSecondary->push_back( tempSecondaryTrack );
	    }
         } //end of loop on secondary 

         // Set the track status according to what the process defined
         fTrack->SetTrackStatus( fParticleChange->GetTrackStatus() );

         // clear ParticleChange
         fParticleChange->Clear();
}
////////////////////////////////////////////////////////


