// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SteppingVerboseWithUnits.cc,v 1.1 1999-07-27 09:21:51 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
// G4SteppingVerboseWithUnits.cc
//
// Description:
//    Implementation of  the G4SteppingVerboseWithUnits class
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

#include "G4SteppingVerboseWithUnits.hh"
#include "G4SteppingManager.hh"

#include "G4UnitsTable.hh"

//////////////////////////////////////////////////
//G4SteppingVerboseWithUnits::G4SteppingVerboseWithUnits()
//////////////////////////////////////////////////
//{
//}
//////////////////////////////////////////////////
G4SteppingVerboseWithUnits::G4SteppingVerboseWithUnits(G4SteppingManager* const fMan)
:G4SteppingVerbose(fMan)
//////////////////////////////////////////////////
{
}
//////////////////////////////////////////////////
G4SteppingVerboseWithUnits::~G4SteppingVerboseWithUnits()
//////////////////////////////////////////////////
{
}

//////////////////////////////////////////////////
void G4SteppingVerboseWithUnits::AtRestDoItInvoked()
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
		<< setw(8)
		<< G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length")
		<< setw(8)
		<< G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length")
		<< setw(8)
		<< G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length")
		<< setw(8)
		<< G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy")
		<< setw(8)
		<< G4BestUnit((*fSecondary)[lp1]->GetGlobalTime(),"Time")
		<< setw(18)
		<< (*fSecondary)[lp1]->GetDefinition()->GetParticleName();
	 G4cout << endl;
       }
     }
   }
   
   if( verboseLevel >= 4 ){ 
     ShowStep();
     G4cout << endl;
   }
}
/////////////////////////////////////////////////////
void G4SteppingVerboseWithUnits::AlongStepDoItAllDone()
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

        ShowStep();
        G4cout << endl;
        G4cout << "    ++List of secondaries generated " 
             << "(x,y,z,kE,t,PID):"
             << "  No. of secodaries = " 
             << (*fSecondary).entries() << endl;

        if((*fSecondary).entries()>0){
           for(G4int lp1=0; lp1<(*fSecondary).entries(); lp1++){
               G4cout << "      "
                 << setw(8)
                 << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length")
                 << setw(8)
                 << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length")
                 << setw(8)
                 << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length")
                 << setw(8)
                 << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy")
                 << setw(8)
                 << G4BestUnit((*fSecondary)[lp1]->GetGlobalTime(),"Time")
                 << setw(18)
                 << (*fSecondary)[lp1]->GetDefinition()->GetParticleName();
               G4cout << endl;
	   }
	}
     }
}
////////////////////////////////////////////////////
void G4SteppingVerboseWithUnits::PostStepDoItAllDone()
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

        ShowStep();
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
                 << setw(8)
                 << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length")
                 << setw(8)
                 << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length")
                 << setw(8)
                 << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length")
                 << setw(8)
                 << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy")
                 << setw(8)
                 << G4BestUnit((*fSecondary)[lp1]->GetGlobalTime(),"Time")
                 << setw(18)
                 << (*fSecondary)[lp1]->GetDefinition()->GetParticleName();
               G4cout << endl;
	     }
	}
      }

 }
/////////////////////////////////////////
void G4SteppingVerboseWithUnits::StepInfo()
/////////////////////////////////////////
{
  CopyState();
  
  G4int prec = G4cout.precision(3);

  if( verboseLevel >= 1 ){
    if( verboseLevel >= 4 ) VerboseTrack();
    if( verboseLevel >= 3 ){
      G4cout << endl;    
      G4cout << setw( 5) << "#Step#"     << " "
	     << setw( 6) << "X"          << "    "
	     << setw( 6) << "Y"          << "    "  
	     << setw( 6) << "Z"          << "    "
	     << setw( 9) << "KineE"      << " "
	     << setw( 9) << "dEStep"     << " "  
	     << setw(10) << "StepLeng"     
	     << setw(10) << "TrakLeng" 
	     << setw(10) << "NextVolu" 
	     << setw(10) << "Process"   << endl;	          
    }

    G4cout << setw( 5) << fTrack->GetCurrentStepNumber() << " "
	   << setw( 6) << G4BestUnit(fTrack->GetPosition().x(),"Length")
	   << setw( 6) << G4BestUnit(fTrack->GetPosition().y(),"Length")
	   << setw( 6) << G4BestUnit(fTrack->GetPosition().z(),"Length")
	   << setw( 6) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
	   << setw( 6) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
	   << setw( 6) << G4BestUnit(fStep->GetStepLength(),"Length")
	   << setw( 6) << G4BestUnit(fTrack->GetTrackLength(),"Length");

    // if( fStepStatus != fWorldBoundary){ 
    if( fTrack->GetNextVolume() != 0 ) { 
      G4cout << setw(10) << fTrack->GetNextVolume()->GetName();
    } else {
      G4cout << setw(10) << "OutOfWorld";
    }

    if(fStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL){
      G4cout << setw(10) << fStep->GetPostStepPoint()->GetProcessDefinedStep()
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
	       << "(Rest="  << setw(2) << fN2ndariesAtRestDoIt
	       << ",Along=" << setw(2) << fN2ndariesAlongStepDoIt
	       << ",Post="  << setw(2) << fN2ndariesPostStepDoIt
	       << "), "
	       << "#SpawnTotal=" << setw(3) << (*fSecondary).entries()
	       << " ---------------"
	       << endl;

	for(G4int lp1=(*fSecondary).entries()-tN2ndariesTot; 
                        lp1<(*fSecondary).entries(); lp1++){
	  G4cout << "    : "
		 << setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length")
		 << setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length")
		 << setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length")
		 << setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy")
		 << setw(10)
		 << (*fSecondary)[lp1]->GetDefinition()->GetParticleName();
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

////////////////////////////////////////////
void G4SteppingVerboseWithUnits::DPSLStarted()
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
void G4SteppingVerboseWithUnits::DPSLUserLimit()
//////////////////////////////////////////////
{
  CopyState();

  if( verboseLevel > 5 ){
    G4cout << endl << endl;
    G4cout << "=== Defined Physical Step Length (DPSL)" << endl;
    G4cout << "    ++ProposedStep(UserLimit) = " 
      << setw(8) << G4BestUnit(physIntLength,"Length")
	<< " : ProcName = User defined maximum allowed Step"
	  << endl;
  }
}
/////////////////////////////////////////////
void G4SteppingVerboseWithUnits::DPSLPostStep()
/////////////////////////////////////////////
{
  CopyState();

  if( verboseLevel > 5 ){
    G4cout << "    ++ProposedStep(PostStep ) = " 
      << setw(8) << G4BestUnit(physIntLength,"Length")
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
void G4SteppingVerboseWithUnits::DPSLAlongStep()
/////////////////////////////////////////////
{
  CopyState();
  if( verboseLevel > 5 ){
    G4cout << "    ++ProposedStep(AlongStep) = " 
	   << setw(8) << G4BestUnit(physIntLength,"Length")
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


////////////////////////////////////////////////
void G4SteppingVerboseWithUnits::TrackingStarted()
////////////////////////////////////////////////
{

  CopyState();
G4int prec = G4cout.precision(3);
  if( verboseLevel > 0 ){

    G4cout << setw( 5) << "Step#"      << " "
           << setw( 6) << "X"          << "    "
	   << setw( 6) << "Y"          << "    "  
	   << setw( 6) << "Z"          << "    "
	   << setw( 9) << "KineE"      << " "
	   << setw( 9) << "dEStep"     << " "  
	   << setw(10) << "StepLeng"  
	   << setw(10) << "TrakLeng"
	   << setw(10) << "NextVolu"
	   << setw(10) << "Process"    << endl;	     

    G4cout << setw( 5) << fTrack->GetCurrentStepNumber() << " "
	   << setw( 6) << G4BestUnit(fTrack->GetPosition().x(),"Length")
	   << setw( 6) << G4BestUnit(fTrack->GetPosition().y(),"Length")
	   << setw( 6) << G4BestUnit(fTrack->GetPosition().z(),"Length")
	   << setw( 6) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
	   << setw( 6) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
	   << setw( 6) << G4BestUnit(fStep->GetStepLength(),"Length")
	   << setw( 6) << G4BestUnit(fTrack->GetTrackLength(),"Length");

    if(fTrack->GetNextVolume()){
      G4cout << setw(10) << fTrack->GetNextVolume()->GetName() << " ";
    } else {
      G4cout << setw(10) << "OutOfWorld" << " ";
    }
    G4cout << setw(10) << "initStep" << endl;
  }
  G4cout.precision(prec);
}
//////////////////////////////////////////////////////
void G4SteppingVerboseWithUnits::AlongStepDoItOneByOne()
//////////////////////////////////////////////////////
{ 
  CopyState();
    if(verboseLevel >= 4){ 
        G4cout << endl;
        G4cout << " >>AlongStepDoIt (process by process): "
             << "   Process Name = " 
             << fCurrentProcess->GetProcessName() << endl;

        ShowStep();
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
                 << setw(8)
                 << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length")
                 << setw(8)
                 << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length")
                 << setw(8)
                 << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length")
                 << setw(8)
                 << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy")
                 << setw(8)
                 << G4BestUnit((*fSecondary)[lp1]->GetGlobalTime(),"Time")
                 << setw(18)
                 << (*fSecondary)[lp1]->GetDefinition()->GetParticleName();
               G4cout << endl;
	   }
	}
     }

}
//////////////////////////////////////////////////////
void G4SteppingVerboseWithUnits::PostStepDoItOneByOne()
//////////////////////////////////////////////////////
{
  CopyState();
     if(fStepStatus != fPostStepDoItProc) return;

     if(verboseLevel >= 4){ 
        G4cout << endl;
        G4cout << " >>PostStepDoIt (process by process): "
             << "   Process Name = " 
             << fCurrentProcess->GetProcessName() << endl;

        ShowStep();
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
                 << setw(8)
                 << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length")
                 << setw(8)
                 << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length")
                 << setw(8)
                 << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length")
                 << setw(8)
                 << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy")
                 << setw(8)
                 << G4BestUnit((*fSecondary)[lp1]->GetGlobalTime(), "Time")
                 << setw(18)
                 << (*fSecondary)[lp1]->GetDefinition()->GetParticleName();
               G4cout << endl;
	   }
	}
     }

}


//////////////////////////////////////
void G4SteppingVerboseWithUnits::VerboseTrack()
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
  G4cout << "        Position - x        : " 
       << setw(20) << G4BestUnit(fTrack->GetPosition().x(),"Length")
       << endl; 
  G4cout << "        Position - y        : " 
       << setw(20) << G4BestUnit(fTrack->GetPosition().y(),"Length")
       << endl; 
  G4cout << "        Position - z        : " 
       << setw(20) << G4BestUnit(fTrack->GetPosition().z(),"Length")
       << endl;
  G4cout << "        Global Time         : " 
       << setw(20) << G4BestUnit(fTrack->GetGlobalTime(),"Time")
       << endl;
  G4cout << "        Local Time          : " 
       << setw(20) << G4BestUnit(fTrack->GetLocalTime(),"Time")
       << endl;
  G4cout << "        Momentum Direct - x : " 
       << setw(20) << fTrack->GetMomentumDirection().x()
       << endl;
  G4cout << "        Momentum Direct - y : " 
       << setw(20) << fTrack->GetMomentumDirection().y()
       << endl;
  G4cout << "        Momentum Direct - z : " 
       << setw(20) << fTrack->GetMomentumDirection().z()
       << endl;
  G4cout << "        Kinetic Energy      : " 
       << setw(20) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
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
       << setw(20) << G4BestUnit(fTrack->GetTrackLength(),"Length")
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
  G4cout << "        Vertex - x          : " 
       << setw(20) << G4BestUnit(fTrack->GetVertexPosition().x(),"Length")
       << endl; 
  G4cout << "        Vertex - y          : " 
       << setw(20) << G4BestUnit(fTrack->GetVertexPosition().y(),"Length")
       << endl; 
  G4cout << "        Vertex - z          : " 
       << setw(20) << G4BestUnit(fTrack->GetVertexPosition().z(),"Length")
       << endl;
  G4cout << "        Vertex - Px (MomDir): " 
       << setw(20) << fTrack->GetVertexMomentumDirection().x()
       << endl;
  G4cout << "        Vertex - Py (MomDir): " 
       << setw(20) << fTrack->GetVertexMomentumDirection().y()
       << endl;
  G4cout << "        Vertex - Pz (MomDir): " 
       << setw(20) << fTrack->GetVertexMomentumDirection().z()
       << endl;
  G4cout << "        Vertex - KineE      : " 
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

