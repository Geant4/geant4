// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SteppingVerbose.cc,v 1.1 1999-01-07 16:14:31 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
// G4SteppingVerbose.cc
//
// Description:
//    Implementation of  the G4SteppingVerbose class
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

#include "G4SteppingVerbose.hh"
#include "G4SteppingManager.hh"

///#define G4_USE_G4BESTUNIT_FOR_VERBOSE 1

#ifdef G4_USE_G4BESTUNIT_FOR_VERBOSE
#include "G4UnitsTable.hh"
#else
#define G4BestUnit(a,b) a
#endif

//////////////////////////////////////////////////
G4SteppingVerbose::G4SteppingVerbose()
//////////////////////////////////////////////////
{
}
//////////////////////////////////////////////////
G4SteppingVerbose::G4SteppingVerbose(G4SteppingManager* fMan)
:fManager(fMan)
//////////////////////////////////////////////////
{
//get the pointer of SteppingManager Object


}
//////////////////////////////////////////////////
G4SteppingVerbose::~G4SteppingVerbose()
//////////////////////////////////////////////////
{
}

//////////////////////////////////////////////////
void G4SteppingVerbose::NewStep()
//////////////////////////////////////////////////
{
}
//////////////////////////////////////////////////
void G4SteppingVerbose::CopyState()
//////////////////////////////////////////////////
{

   fUserSteppingAction = fManager->GetUserAction();
   fVerbose = this;

   PhysicalStep = fManager->GetPhysicalStep();
   GeometricalStep = fManager->GetGeometricalStep();
   CorrectedStep = fManager->GetCorrectedStep();
   PreStepPointIsGeom = fManager->GetPreStepPointIsGeom();
   FirstStep = fManager->GetFirstStep();
   fStepStatus = fManager->GetfStepStatus();

   TempInitVelocity = fManager->GetTempInitVelocity();
   TempVelocity = fManager->GetTempVelocity();
   Mass = fManager->GetMass();

   sumEnergyChange = fManager->GetsumEnergyChange();

   fParticleChange = fManager->GetfParticleChange();
   fTrack = fManager->GetfTrack(); 
   fSecondary = fManager->GetfSecondary();
   fStep = fManager->GetfStep();
   fPreStepPoint = fManager->GetfPreStepPoint();
   fPostStepPoint = fManager->GetfPostStepPoint();

   fCurrentVolume = fManager->GetfCurrentVolume();
   fSensitive = fManager->GetfSensitive();
   fCurrentProcess = fManager->GetfCurrentProcess();

   fAtRestDoItVector = fManager->GetfAtRestDoItVector(); 
   fAlongStepDoItVector = fManager->GetfAlongStepDoItVector();
   fPostStepDoItVector = fManager->GetfPostStepDoItVector();

   fAtRestGetPhysIntVector = fManager->GetfAtRestGetPhysIntVector();
   fAlongStepGetPhysIntVector = fManager->GetfAlongStepGetPhysIntVector();
   fPostStepGetPhysIntVector = fManager->GetfPostStepGetPhysIntVector();

   MAXofAtRestLoops = fManager->GetMAXofAtRestLoops();
   MAXofAlongStepLoops = fManager->GetMAXofAlongStepLoops();
   MAXofPostStepLoops = fManager->GetMAXofPostStepLoops();

   currentMinimumStep = fManager->GetcurrentMinimumStep();
   numberOfInteractionLengthLeft = fManager->GetnumberOfInteractionLengthLeft();

   fAtRestDoItProcTriggered = fManager->GetfAtRestDoItProcTriggered();
   fAlongStepDoItProcTriggered = fManager->GetfAlongStepDoItProcTriggered();
   fPostStepDoItProcTriggered = fManager->GetfPostStepDoItProcTriggered();

   fN2ndariesAtRestDoIt = fManager->GetfN2ndariesAtRestDoIt();
   fN2ndariesAlongStepDoIt = fManager->GetfN2ndariesAlongStepDoIt();
   fN2ndariesPostStepDoIt = fManager->GetfN2ndariesPostStepDoIt();

   fNavigator = fManager->GetfNavigator();

   verboseLevel = fManager->GetverboseLevel();

   fSelectedAtRestDoItVector = fManager->GetfSelectedAtRestDoItVector();
   fSelectedAlongStepDoItVector = fManager->GetfSelectedAlongStepDoItVector();
   fSelectedPostStepDoItVector = fManager->GetfSelectedPostStepDoItVector();

   fPreviousStepSize = fManager->GetfPreviousStepSize();

   fTouchable1 = fManager->GetfTouchable1();
   fTouchable2 = fManager->GetfTouchable2();
   fIsTouchable1Free = fManager->GetfIsTouchable1Free();
   fIsTouchable2Free = fManager->GetfIsTouchable2Free();

   StepControlFlag = fManager->GetStepControlFlag();

   physIntLength = fManager->GetphysIntLength();
   fCondition = fManager->GetfCondition();
   fGPILSelection = fManager->GetfGPILSelection();
}


//////////////////////////////////////////////////
void G4SteppingVerbose::AtRestDoItInvoked()
//////////////////////////////////////////////////
 {
   G4VProcess* ptProcManager;

   CopyState();

   if(verboseLevel >= 3 ){     
     G4int npt=0;
     G4cout << " **List of AtRestDoIt invoked:" << endl;
     for(size_t np=0; np < MAXofAtRestLoops; np++){
       size_t npGPIL = MAXofAtRestLoops-np-1;
       if( (*fSelectedAtRestDoItVector)(npGPIL) == 2 ){
	 npt++;                
	 ptProcManager = (*fAtRestDoItVector)(np);
	 G4cout << "   # " << npt << " : " 
	   << ptProcManager->GetProcessName() 
	     << " (Forced)" << endl;
       } else if ( (*fSelectedAtRestDoItVector)(npGPIL) == 1 ){
	 npt++;                
	 ptProcManager = (*fAtRestDoItVector)(np);
	 G4cout << "   # " << npt << " : " 
	   << ptProcManager->GetProcessName() << endl;
       }
     }
     
     G4cout << "   Generated secondries # : " << fN2ndariesAtRestDoIt << endl;
     
     if( fN2ndariesAtRestDoIt > 0 ){
       G4cout << "   -- List of secondaries generated : " 
	 << "(x,y,z,kE,t,PID) --" << endl; 
       for( G4int lp1=(*fSecondary).entries()-fN2ndariesAtRestDoIt; 
	   lp1<(*fSecondary).entries(); lp1++){
	 G4cout << "      "
		<< setw( 9)
		<< G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length") << " "
		<< setw( 9)
		<< G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length") << " "
		<< setw( 9)
		<< G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length") << " "
		<< setw( 9)
		<< G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy") << " "
		<< setw( 9)
		<< G4BestUnit((*fSecondary)[lp1]->GetGlobalTime(),"Time") << " "
		<< setw(18)
		<< (*fSecondary)[lp1]->GetDefinition()
	                             ->GetParticleName();
	 G4cout << endl;
       }
     }
   }
   
   if( verboseLevel >= 4 ){ 
     fStep->ShowStep();
     G4cout << endl;
   }
}
/////////////////////////////////////////////////////
void G4SteppingVerbose::AlongStepDoItAllDone()
/////////////////////////////////////////////////////
{
   G4VProcess* ptProcManager;

   CopyState();

     if(verboseLevel >= 3){ 
        G4cout << endl;
        G4cout << " >>AlongStepDoIt (after all invocations):" << endl;
        G4cout << "    ++List of invoked processes " << endl;
 
        for(size_t ci=0; ci<MAXofAlongStepLoops; ci++){
            ptProcManager = (*fAlongStepDoItVector)(ci);
            G4cout << "      " << ci+1 << ") ";
            if(ptProcManager != NULL){
               G4cout << ptProcManager->GetProcessName() << endl;
            }
        }         

        fStep->ShowStep();
        G4cout << endl;
        G4cout << "    ++List of secondaries generated " 
             << "(x,y,z,kE,t,PID):"
             << "  No. of secodaries = " 
             << (*fSecondary).entries() << endl;

        if((*fSecondary).entries()>0){
           for(G4int lp1=0; lp1<(*fSecondary).entries(); lp1++){
               G4cout << "      "
                    << setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length") << " "
                    << setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length") << " "
                    << setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length") << " "
                    << setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy") << " "
                    << setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetGlobalTime(),"Time")  << " "
                    << setw(18)
                    << (*fSecondary)[lp1]->GetDefinition()
                                         ->GetParticleName();
               G4cout << endl;
	   }
	}
     }
}
////////////////////////////////////////////////////
void G4SteppingVerbose::PostStepDoItAllDone()
////////////////////////////////////////////////////
{
   G4VProcess* ptProcManager;

   CopyState();

   if(fStepStatus != fPostStepDoItProc) return;

   if(verboseLevel >= 3){ 
        G4int npt=0;
        G4cout << endl;
        G4cout << " **PostStepDoIt (after all invocations):" << endl;
        G4cout << "    ++List of invoked processes " << endl;

        for(size_t np=0; np < MAXofPostStepLoops; np++){
	    size_t npGPIL = MAXofPostStepLoops-np-1;
            if((*fSelectedPostStepDoItVector)(npGPIL) == 2){
               npt++;                
               ptProcManager = (*fPostStepDoItVector)(np);
               G4cout << "      " << npt << ") " 
                    << ptProcManager->GetProcessName()  
                    << " (Forced)" << endl;
	     } else if ((*fSelectedPostStepDoItVector)(npGPIL) == 1){
               npt++;                
               ptProcManager = (*fPostStepDoItVector)(np);
               G4cout << "      " << npt << ") " 
                    << ptProcManager->GetProcessName() << endl;
	     }
	  }

        fStep->ShowStep();
        G4cout << endl;
        G4cout << "    ++List of secondaries generated " 
             << "(x,y,z,kE,t,PID):"
             << "  No. of secodaries = " 
             << (*fSecondary).entries() << endl;
        G4cout << "      [Note]Secondaries from AlongStepDoIt included."
             << endl; 

        if((*fSecondary).entries()>0){
	  for(G4int lp1=0; lp1<(*fSecondary).entries(); lp1++){
               G4cout << "      "
                    << setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetPosition().x() , "Length") << " "
                    << setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetPosition().y() , "Length") << " "
                    << setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetPosition().z() , "Length") << " "
                    << setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy() , "Energy") << " "
                    << setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetGlobalTime() , "Time") << " "
                    << setw(18)
                    << (*fSecondary)[lp1]->GetDefinition()
                                         ->GetParticleName();
               G4cout << endl;
	     }
	}
      }

 }
/////////////////////////////////////////
void G4SteppingVerbose::StepInfo()
/////////////////////////////////////////
{
  CopyState();
  
  G4int prec = G4cout.precision(3);

  if( verboseLevel >= 1 ){
    if( verboseLevel >= 4 ) VerboseTrack();
    if( verboseLevel >= 3 ){
      G4cout << endl;
#ifdef G4_USE_G4BESTUNIT_FOR_VERBOSE      
      G4cout << setw( 5) << "#Step#" << " "
	     << setw( 8) << "X"      << "     "
	     << setw( 8) << "Y"      << "     "  
	     << setw( 8) << "Z"      << "     "
	     << setw( 9) << "KineE"  << "     "
	     << setw( 8) << "dE"     << "     "  
	     << setw(12) << "StepLeng"   << " "  
	     << setw(12) << "TrackLeng"  << " "
	     << setw(12) << "NextVolume" << " "
	     << setw( 8) << "ProcName"   << endl;	       
#else
      G4cout << setw( 5) << "#Step#"     << " "
	     << setw( 8) << "X(mm)"      << " "
	     << setw( 8) << "Y(mm)"      << " "  
	     << setw( 8) << "Z(mm)"      << " "
	     << setw( 9) << "KinE(MeV)"  << " "
	     << setw( 8) << "dE(MeV)"    << " "  
	     << setw( 8) << "StepLeng"   << " "  
	     << setw( 9) << "TrackLeng"  << " "  
	     << setw(11) << "NextVolume" << " "
	     << setw( 8) << "ProcName"   << endl;
#endif	     
    }

    G4cout << setw( 5) << fTrack->GetCurrentStepNumber() << " "
	   << setw( 8) << G4BestUnit(fTrack->GetPosition().x() , "Length") << " "
	   << setw( 8) << G4BestUnit(fTrack->GetPosition().y() , "Length") << " "
	   << setw( 8) << G4BestUnit(fTrack->GetPosition().z() , "Length") << " "
	   << setw( 9) << G4BestUnit(fTrack->GetKineticEnergy() , "Energy") << " "
	   << setw( 8) << G4BestUnit(fStep->GetTotalEnergyDeposit(), "Energy") << " "
	   << setw( 8) << G4BestUnit(fStep->GetStepLength() , "Length") << " "
	   << setw( 9) << G4BestUnit(fTrack->GetTrackLength() , "Length") << " ";

    // if( fStepStatus != fWorldBoundary){ 
    if( fTrack->GetNextVolume() != 0 ) { 
      G4cout << setw(11) << fTrack->GetNextVolume()->GetName() << " ";
    } else {
      G4cout << setw(11) << "OutOfWorld" << " ";
    }

    if(fStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL){
      G4cout << fStep->GetPostStepPoint()->GetProcessDefinedStep()
	->GetProcessName();
    } else {
      G4cout << "User Limit";
    }

    G4cout << endl;

    if( verboseLevel == 2 ){
      G4int tN2ndariesTot = fN2ndariesAtRestDoIt +
	                    fN2ndariesAlongStepDoIt +
	                    fN2ndariesPostStepDoIt;
      if(tN2ndariesTot>0){
	G4cout << "    :----- List of 2ndaries - "
	       << "#SpawnInStep=" << setw(3) << tN2ndariesTot 
	       << "(Rest=" << setw(2) << fN2ndariesAtRestDoIt
	       << ",Along=" << setw(2) << fN2ndariesAlongStepDoIt
	       << ",Post="  << setw(2) << fN2ndariesPostStepDoIt
	       << "), "
	       << "#SpawnTotal=" << setw(3) << (*fSecondary).entries()
	       << " ---------------"
	       << endl;

	for(G4int lp1=(*fSecondary).entries()-tN2ndariesTot; 
                        lp1<(*fSecondary).entries(); lp1++){
	  G4cout << "    : "
		 << setw( 9)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().x() , "Length")<< " "
		 << setw( 9)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().y() , "Length")<< " "
		 << setw( 9)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().z() , "Length") << " "
		 << setw( 9)
		 << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy() , "Energy")<< " "
		 << setw(18)
		 << (*fSecondary)[lp1]->GetDefinition()
 	                              ->GetParticleName();
	  G4cout << endl;
	}
              
	G4cout << "    :-----------------------------"
	       << "----------------------------------"
	       << "-- EndOf2ndaries Info ---------------"
	       << endl;
      }
    }
    
  }
  G4cout.precision(prec);
}

////////////////////////////////////////////////
void G4SteppingVerbose::TrackingStarted()
////////////////////////////////////////////////
{

  CopyState();
G4int prec = G4cout.precision(3);
  if( verboseLevel > 0 ){

#ifdef G4_USE_G4BESTUNIT_FOR_VERBOSE
    G4cout << setw( 5) << "Step#"  << " "
           << setw( 8) << "X"      << "     "
	   << setw( 8) << "Y"      << "     "  
	   << setw( 8) << "Z"      << "     "
	   << setw( 9) << "KineE"  << "     "
	   << setw( 8) << "dE"     << "     "  
	   << setw(12) << "StepLeng"   << " "  
	   << setw(12) << "TrackLeng"  << " "
	   << setw(12) << "NextVolume" << " "
	   << setw( 8) << "ProcName"   << endl;	     
#else
    G4cout << setw( 5) << "Step#"      << " "
	   << setw( 8) << "X(mm)"      << " "
	   << setw( 8) << "Y(mm)"      << " "  
	   << setw( 8) << "Z(mm)"      << " "
	   << setw( 9) << "KinE(MeV)"  << " "
	   << setw( 8) << "dE(MeV)"    << " "  
	   << setw( 8) << "StepLeng"   << " "  
	   << setw( 9) << "TrackLeng"  << " "
	   << setw(11) << "NextVolume" << " "
	   << setw( 8) << "ProcName"   << endl;	     
#endif

    G4cout << setw( 5) << fTrack->GetCurrentStepNumber() << " "
	   << setw( 8) << G4BestUnit(fTrack->GetPosition().x(),"Length")<< " "
	   << setw( 8) << G4BestUnit(fTrack->GetPosition().y(),"Length") << " "
	   << setw( 8) << G4BestUnit(fTrack->GetPosition().z(),"Length")<< " "
	   << setw( 9) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")<< " "
	   << setw( 8) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy") << " "
	   << setw( 8) << G4BestUnit(fStep->GetStepLength(),"Length")<< " "
	   << setw( 9) << G4BestUnit(fTrack->GetTrackLength(),"Length") << " ";

    if(fTrack->GetNextVolume()){
      G4cout << setw(11) << fTrack->GetNextVolume()->GetName() << " ";
    } else {
      G4cout << setw(11) << "OutOfWorld" << " ";
    }
    G4cout << "initStep" << endl;
  }
  G4cout.precision(prec);
}
////////////////////////////////////////////
void G4SteppingVerbose::DPSLStarted()
////////////////////////////////////////////
{
  CopyState();

  if( verboseLevel > 5 ){
    G4cout << endl;
    G4cout << " >>DefinePhysicalStepLength (List of proposed StepLengths): "
      << endl;
  }
}
//////////////////////////////////////////////
void G4SteppingVerbose::DPSLUserLimit()
//////////////////////////////////////////////
{
  CopyState();

  if( verboseLevel > 5 ){
    G4cout << endl << endl;
    G4cout << "=== Defined Physical Step Length (DPSL)" << endl;
    G4cout << "    ++ProposedStep(UserLimit) = " 
      << setw( 9) << physIntLength
	<< " : ProcName = User defined maximum allowed Step"
	  << endl;
  }
}
/////////////////////////////////////////////
void G4SteppingVerbose::DPSLPostStep()
/////////////////////////////////////////////
{
  CopyState();

  if( verboseLevel > 5 ){
    G4cout << "    ++ProposedStep(PostStep ) = " 
      << setw( 9) << physIntLength
	<< " : ProcName = "
	  << fCurrentProcess->GetProcessName() 
            << " (";
    if(fCondition==ExclusivelyForced){
      G4cout << "ExclusivelyForced)" << endl;
    }
    else if(fCondition==Conditionally){
      G4cout << "Conditionally)" << endl;
    }
    else if(fCondition==Forced){
      G4cout << "Forced)" << endl;
    }
    else{
      G4cout << "No ForceCondition)" << endl;
    }
  }
}
/////////////////////////////////////////////
void G4SteppingVerbose::DPSLAlongStep()
/////////////////////////////////////////////
{
  CopyState();
  if( verboseLevel > 5 ){
    G4cout << "    ++ProposedStep(AlongStep) = " 
	   << setw( 9) << G4BestUnit(physIntLength , "Length")
	   << " : ProcName = "
	   << fCurrentProcess->GetProcessName() 
	   << " (";
    if(fGPILSelection==CandidateForSelection){
      G4cout << "CandidateForSelection)" << endl;
    }
    else if(fGPILSelection==NotCandidateForSelection){
      G4cout << "NotCandidateForSelection)" << endl;
    }
    else{
      G4cout << "???)" << endl;
    }
  }
}
//////////////////////////////////////////////////////
void G4SteppingVerbose::AlongStepDoItOneByOne()
//////////////////////////////////////////////////////
{ 
  CopyState();
    if(verboseLevel >= 4){ 
        G4cout << endl;
        G4cout << " >>AlongStepDoIt (process by process): "
             << "   Process Name = " 
             << fCurrentProcess->GetProcessName() << endl;

        fStep->ShowStep();
        G4cout << "          "
	       << "!Note! Safety of PostStep is only valid "
	       << "after all DoIt invocations."
	       << endl; 

        VerboseParticleChange();    
        G4cout << endl;

        G4cout << "    ++List of secondaries generated " 
	       << "(x,y,z,kE,t,PID):"
	       << "  No. of secodaries = " 
	       << fN2ndariesAlongStepDoIt << endl;

        if(fN2ndariesAlongStepDoIt>0){
           for(G4int lp1=(*fSecondary).entries()-fN2ndariesAlongStepDoIt; 
                     lp1<(*fSecondary).entries(); lp1++){
               G4cout << "      "
                    << setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetPosition().x() , "Length")<< " "
                    << setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetPosition().y() , "Length")<< " "
                    << setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetPosition().z() , "Length")<< " "
                    << setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy() , "Energy")<< " "
                    << setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetGlobalTime() , "Time")<< " "
                    << setw(18)
                    << (*fSecondary)[lp1]->GetDefinition()
                                         ->GetParticleName();
               G4cout << endl;
	   }
	}
     }

}
//////////////////////////////////////////////////////
void G4SteppingVerbose::PostStepDoItOneByOne()
//////////////////////////////////////////////////////
{
  CopyState();
     if(fStepStatus != fPostStepDoItProc) return;

     if(verboseLevel >= 4){ 
        G4cout << endl;
        G4cout << " >>PostStepDoIt (process by process): "
             << "   Process Name = " 
             << fCurrentProcess->GetProcessName() << endl;

        fStep->ShowStep();
        G4cout << endl;
        VerboseParticleChange();    
        G4cout << endl;
         
        G4cout << "    ++List of secondaries generated " 
             << "(x,y,z,kE,t,PID):"
             << "  No. of secodaries = " 
             << fN2ndariesPostStepDoIt << endl;

        if(fN2ndariesPostStepDoIt>0){
           for(G4int lp1=(*fSecondary).entries()-fN2ndariesPostStepDoIt; 
                     lp1<(*fSecondary).entries(); lp1++){
               G4cout << "      "
                    << setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetPosition().x() , "Length")<< " "
                    << setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(), "Length") << " "
                    << setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(), "Length") << " "
                    << setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(), "Energy") << " "
                    << setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetGlobalTime(), "Time") << " "
                    << setw(18)
                    << (*fSecondary)[lp1]->GetDefinition()
                                         ->GetParticleName();
               G4cout << endl;
	   }
	}
     }

}


//////////////////////////////////////
void G4SteppingVerbose::VerboseTrack()
//////////////////////////////////////
{
  CopyState();
// Show header
  G4cout << endl;
  G4cout << "    ++G4Track Information " << endl;
  G4int prec = G4cout.precision(3);


  G4cout << "      -----------------------------------------------" 
       << endl;
  G4cout << "        G4Track Information  " << setw(20) << endl;
  G4cout << "      -----------------------------------------------" 
       << endl;

  G4cout << "        Step number         : " 
       << setw(20) << fTrack->GetCurrentStepNumber()
       << endl; 
#ifdef G4_USE_G4BESTUNIT_FOR_VERBOSE
  G4cout << "        Position - x        : " 
       << setw(20) << G4BestUnit(fTrack->GetPosition().x(), "Length")
       << endl; 
  G4cout << "        Position - y        : " 
       << setw(20) << G4BestUnit(fTrack->GetPosition().y(), "Length")
       << endl; 
  G4cout << "        Position - z        : " 
       << setw(20) << G4BestUnit(fTrack->GetPosition().z(), "Length")
       << endl;
  G4cout << "        Global Time         : " 
       << setw(20) << G4BestUnit(fTrack->GetGlobalTime(), "Time")
       << endl;
  G4cout << "        Local Time          : " 
       << setw(20) << G4BestUnit(fTrack->GetLocalTime(), "Time")
       << endl;
#else
  G4cout << "        Position - x (mm)   : " 
       << setw(20) << fTrack->GetPosition().x() /mm
       << endl; 
  G4cout << "        Position - y (mm)   : " 
       << setw(20) << fTrack->GetPosition().y() /mm
       << endl; 
  G4cout << "        Position - z (mm)   : " 
       << setw(20) << fTrack->GetPosition().z() /mm
       << endl;
  G4cout << "        Global Time (ns)    : " 
       << setw(20) << fTrack->GetGlobalTime() /ns
       << endl;
  G4cout << "        Local Time (ns)     : " 
       << setw(20) << fTrack->GetLocalTime() /ns
       << endl;
#endif
  G4cout << "        Momentum Direct - x : " 
       << setw(20) << fTrack->GetMomentumDirection().x()
       << endl;
  G4cout << "        Momentum Direct - y : " 
       << setw(20) << fTrack->GetMomentumDirection().y()
       << endl;
  G4cout << "        Momentum Direct - z : " 
       << setw(20) << fTrack->GetMomentumDirection().z()
       << endl;
#ifdef G4_USE_G4BESTUNIT_FOR_VERBOSE
  G4cout << "        Kinetic Energy      : " 
#else
  G4cout << "        Kinetic Energy (MeV): " 
#endif
       << setw(20) << G4BestUnit(fTrack->GetKineticEnergy(), "Energy")
       << endl;
  G4cout << "        Polarization - x    : " 
       << setw(20) << fTrack->GetPolarization().x()
       << endl;
  G4cout << "        Polarization - y    : " 
       << setw(20) << fTrack->GetPolarization().y()
       << endl;
  G4cout << "        Polarization - z    : " 
       << setw(20) << fTrack->GetPolarization().z()
       << endl;
  G4cout << "        Track Length        : " 
       << setw(20) << G4BestUnit(fTrack->GetTrackLength(), "Length")
       << endl;
  G4cout << "        Track ID #          : " 
       << setw(20) << fTrack->GetTrackID()
       << endl;
  G4cout << "        Parent Track ID #   : " 
       << setw(20) << fTrack->GetParentID()
       << endl;
  G4cout << "        Next Volume         : " 
       << setw(20);
       if( fTrack->GetNextVolume() != 0 ) { 
         G4cout << fTrack->GetNextVolume()->GetName() << " ";
       } else {
         G4cout << "OutOfWorld" << " ";
       }
       G4cout << endl;
  G4cout << "        Track Status        : " 
       << setw(20);
       if( fTrack->GetTrackStatus() == fAlive ){
         G4cout << " Alive";
       } else if( fTrack->GetTrackStatus() == fStopButAlive ){
           G4cout << " StopButAlive";
       } else if( fTrack->GetTrackStatus() == fStopAndKill ){
           G4cout << " StopAndKill";
       } else if( fTrack->GetTrackStatus() == fKillTrackAndSecondaries ){
           G4cout << " KillTrackAndSecondaries";
       } else if( fTrack->GetTrackStatus() == fSuspend ){
           G4cout << " Suspend";
       } else if( fTrack->GetTrackStatus() == fPostponeToNextEvent ){
           G4cout << " PostponeToNextEvent";
       }
       G4cout << endl;
#ifdef G4_USE_G4BESTUNIT_FOR_VERBOSE
  G4cout << "        Vertex - x          : " 
       << setw(20) << G4BestUnit(fTrack->GetVertexPosition().x(),"Length")
       << endl; 
  G4cout << "        Vertex - y          : " 
       << setw(20) << G4BestUnit(fTrack->GetVertexPosition().y(),"Length")
       << endl; 
  G4cout << "        Vertex - z          : " 
       << setw(20) << G4BestUnit(fTrack->GetVertexPosition().z(),"Length")
       << endl;
#else
  G4cout << "        Vertex - x (mm)     : " 
       << setw(20) << fTrack->GetVertexPosition().x()/mm
       << endl; 
  G4cout << "        Vertex - y (mm)     : " 
       << setw(20) << fTrack->GetVertexPosition().y()/mm
       << endl; 
  G4cout << "        Vertex - z (mm)     : " 
       << setw(20) << fTrack->GetVertexPosition().z()/mm
       << endl;
#endif
  G4cout << "        Vertex - Px (MomDir): " 
       << setw(20) << fTrack->GetVertexMomentumDirection().x()
       << endl;
  G4cout << "        Vertex - Py (MomDir): " 
       << setw(20) << fTrack->GetVertexMomentumDirection().y()
       << endl;
  G4cout << "        Vertex - Pz (MomDir): " 
       << setw(20) << fTrack->GetVertexMomentumDirection().z()
       << endl;
#ifdef G4_USE_G4BESTUNIT_FOR_VERBOSE
  G4cout << "        Vertex - KineE      : " 
#else
  G4cout << "        Vertex - KineE (MeV): " 
#endif
       << setw(20) << G4BestUnit(fTrack->GetVertexKineticEnergy(),"Energy")
       << endl;
  
  G4cout << "        Creator Process     : " 
       << setw(20);
  if( fTrack->GetCreatorProcess() == NULL){
    G4cout << " Event Generator" << endl;
  } else {
    G4cout << fTrack->GetCreatorProcess()->GetProcessName() << endl;
  }

  G4cout << "      -----------------------------------------------" 
       << endl;
       
 G4cout.precision(prec);      
}


///////////////////////////////////////////////
void G4SteppingVerbose::VerboseParticleChange()
///////////////////////////////////////////////
{

// Show header
  G4cout << endl;
  G4cout << "    ++G4ParticleChange Information " << endl;
  fParticleChange->DumpInfo();
}
