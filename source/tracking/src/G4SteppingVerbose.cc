// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SteppingVerbose.cc,v 1.7 2001-02-08 07:39:52 tsasaki Exp $
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

#include "G4VSensitiveDetector.hh"    // Include from 'hits/digi'
#include "G4StepStatus.hh"    // Include from 'tracking'

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
#ifdef G4_TRACKING_DEBUG
   G4cout << "G4SteppingVerbose has instantiated" << G4endl;
#endif
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
void G4SteppingVerbose::AtRestDoItInvoked()
//////////////////////////////////////////////////
 {
   G4VProcess* ptProcManager;

   CopyState();

   if(verboseLevel >= 3 ){     
     G4int npt=0;
     G4cout << " **List of AtRestDoIt invoked:" << G4endl;
     for(size_t np=0; np < MAXofAtRestLoops; np++){
       size_t npGPIL = MAXofAtRestLoops-np-1;
       if( (*fSelectedAtRestDoItVector)[npGPIL] == 2 ){
	 npt++;                
	 ptProcManager = (*fAtRestDoItVector)[np];
	 G4cout << "   # " << npt << " : " 
	   << ptProcManager->GetProcessName() 
	     << " (Forced)" << G4endl;
       } else if ( (*fSelectedAtRestDoItVector)[npGPIL] == 1 ){
	 npt++;                
	 ptProcManager = (*fAtRestDoItVector)[np];
	 G4cout << "   # " << npt << " : " 
	   << ptProcManager->GetProcessName() << G4endl;
       }
     }
     
     G4cout << "   Generated secondries # : " << fN2ndariesAtRestDoIt << G4endl;
     
     if( fN2ndariesAtRestDoIt > 0 ){
       G4cout << "   -- List of secondaries generated : " 
	 << "(x,y,z,kE,t,PID) --" << G4endl; 
       for( G4int lp1=(*fSecondary).size()-fN2ndariesAtRestDoIt; 
	   lp1<(*fSecondary).size(); lp1++){
	 G4cout << "      "
		<< G4std::setw( 9)
		<< G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length") << " "
		<< G4std::setw( 9)
		<< G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length") << " "
		<< G4std::setw( 9)
		<< G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length") << " "
		<< G4std::setw( 9)
		<< G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy") << " "
		<< G4std::setw( 9)
		<< G4BestUnit((*fSecondary)[lp1]->GetGlobalTime(),"Time") << " "
		<< G4std::setw(18)
		<< (*fSecondary)[lp1]->GetDefinition()
	                             ->GetParticleName();
	 G4cout << G4endl;
       }
     }
   }
   
   if( verboseLevel >= 4 ){ 
     ShowStep();
     G4cout << G4endl;
   }
}
/////////////////////////////////////////////////////
void G4SteppingVerbose::AlongStepDoItAllDone()
/////////////////////////////////////////////////////
{
   G4VProcess* ptProcManager;

   CopyState();

     if(verboseLevel >= 3){ 
        G4cout << G4endl;
        G4cout << " >>AlongStepDoIt (after all invocations):" << G4endl;
        G4cout << "    ++List of invoked processes " << G4endl;
 
        for(size_t ci=0; ci<MAXofAlongStepLoops; ci++){
            ptProcManager = (*fAlongStepDoItVector)(ci);
            G4cout << "      " << ci+1 << ") ";
            if(ptProcManager != NULL){
               G4cout << ptProcManager->GetProcessName() << G4endl;
            }
        }         

        ShowStep();
        G4cout << G4endl;
        G4cout << "    ++List of secondaries generated " 
             << "(x,y,z,kE,t,PID):"
             << "  No. of secodaries = " 
             << (*fSecondary).size() << G4endl;

        if((*fSecondary).size()>0){
           for(G4int lp1=0; lp1<(*fSecondary).size(); lp1++){
               G4cout << "      "
                    << G4std::setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length") << " "
                    << G4std::setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length") << " "
                    << G4std::setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length") << " "
                    << G4std::setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy") << " "
                    << G4std::setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetGlobalTime(),"Time")  << " "
                    << G4std::setw(18)
                    << (*fSecondary)[lp1]->GetDefinition()
                                         ->GetParticleName();
               G4cout << G4endl;
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
        G4cout << G4endl;
        G4cout << " **PostStepDoIt (after all invocations):" << G4endl;
        G4cout << "    ++List of invoked processes " << G4endl;

        for(size_t np=0; np < MAXofPostStepLoops; np++){
	    size_t npGPIL = MAXofPostStepLoops-np-1;
            if( (*fSelectedPostStepDoItVector)[npGPIL] == 2){
               npt++;                
               ptProcManager = (*fPostStepDoItVector)[np];
               G4cout << "      " << npt << ") " 
                    << ptProcManager->GetProcessName()  
                    << " (Forced)" << G4endl;
	     } else if ( (*fSelectedPostStepDoItVector)[npGPIL] == 1){
               npt++;                
               ptProcManager = (*fPostStepDoItVector)[np];
               G4cout << "      " << npt << ") " 
                    << ptProcManager->GetProcessName() << G4endl;
	     }
	  }

        ShowStep();
        G4cout << G4endl;
        G4cout << "    ++List of secondaries generated " 
             << "(x,y,z,kE,t,PID):"
             << "  No. of secodaries = " 
             << (*fSecondary).size() << G4endl;
        G4cout << "      [Note]Secondaries from AlongStepDoIt included."
             << G4endl; 

        if((*fSecondary).size()>0){
	  for(G4int lp1=0; lp1<(*fSecondary).size(); lp1++){
               G4cout << "      "
                    << G4std::setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetPosition().x() , "Length") << " "
                    << G4std::setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetPosition().y() , "Length") << " "
                    << G4std::setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetPosition().z() , "Length") << " "
                    << G4std::setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy() , "Energy") << " "
                    << G4std::setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetGlobalTime() , "Time") << " "
                    << G4std::setw(18)
                    << (*fSecondary)[lp1]->GetDefinition()
                                         ->GetParticleName();
               G4cout << G4endl;
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
      G4cout << G4endl;
#ifdef G4_USE_G4BESTUNIT_FOR_VERBOSE      
      G4cout << G4std::setw( 5) << "#Step#" << " "
	     << G4std::setw( 8) << "X"      << "     "
	     << G4std::setw( 8) << "Y"      << "     "  
	     << G4std::setw( 8) << "Z"      << "     "
	     << G4std::setw( 9) << "KineE"  << "     "
	     << G4std::setw( 8) << "dE"     << "     "  
	     << G4std::setw(12) << "StepLeng"   << " "  
	     << G4std::setw(12) << "TrackLeng"  << " "
	     << G4std::setw(12) << "NextVolume" << " "
	     << G4std::setw( 8) << "ProcName"   << G4endl;	       
#else
      G4cout << G4std::setw( 5) << "#Step#"     << " "
	     << G4std::setw( 8) << "X(mm)"      << " "
	     << G4std::setw( 8) << "Y(mm)"      << " "  
	     << G4std::setw( 8) << "Z(mm)"      << " "
	     << G4std::setw( 9) << "KinE(MeV)"  << " "
	     << G4std::setw( 8) << "dE(MeV)"    << " "  
	     << G4std::setw( 8) << "StepLeng"   << " "  
	     << G4std::setw( 9) << "TrackLeng"  << " "  
	     << G4std::setw(11) << "NextVolume" << " "
	     << G4std::setw( 8) << "ProcName"   << G4endl;
#endif	     
    }

    G4cout << G4std::setw( 5) << fTrack->GetCurrentStepNumber() << " "
	   << G4std::setw( 8) << G4BestUnit(fTrack->GetPosition().x() , "Length") << " "
	   << G4std::setw( 8) << G4BestUnit(fTrack->GetPosition().y() , "Length") << " "
	   << G4std::setw( 8) << G4BestUnit(fTrack->GetPosition().z() , "Length") << " "
	   << G4std::setw( 9) << G4BestUnit(fTrack->GetKineticEnergy() , "Energy") << " "
	   << G4std::setw( 8) << G4BestUnit(fStep->GetTotalEnergyDeposit(), "Energy") << " "
	   << G4std::setw( 8) << G4BestUnit(fStep->GetStepLength() , "Length") << " "
	   << G4std::setw( 9) << G4BestUnit(fTrack->GetTrackLength() , "Length") << " ";

    // if( fStepStatus != fWorldBoundary){ 
    if( fTrack->GetNextVolume() != 0 ) { 
      G4cout << G4std::setw(11) << fTrack->GetNextVolume()->GetName() << " ";
    } else {
      G4cout << G4std::setw(11) << "OutOfWorld" << " ";
    }

    if(fStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL){
      G4cout << fStep->GetPostStepPoint()->GetProcessDefinedStep()
	->GetProcessName();
    } else {
      G4cout << "User Limit";
    }

    G4cout << G4endl;

    if( verboseLevel == 2 ){
      G4int tN2ndariesTot = fN2ndariesAtRestDoIt +
	                    fN2ndariesAlongStepDoIt +
	                    fN2ndariesPostStepDoIt;
      if(tN2ndariesTot>0){
	G4cout << "    :----- List of 2ndaries - "
	       << "#SpawnInStep=" << G4std::setw(3) << tN2ndariesTot 
	       << "(Rest=" << G4std::setw(2) << fN2ndariesAtRestDoIt
	       << ",Along=" << G4std::setw(2) << fN2ndariesAlongStepDoIt
	       << ",Post="  << G4std::setw(2) << fN2ndariesPostStepDoIt
	       << "), "
	       << "#SpawnTotal=" << G4std::setw(3) << (*fSecondary).size()
	       << " ---------------"
	       << G4endl;

	for(G4int lp1=(*fSecondary).size()-tN2ndariesTot; 
                        lp1<(*fSecondary).size(); lp1++){
	  G4cout << "    : "
		 << G4std::setw( 9)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().x() , "Length")<< " "
		 << G4std::setw( 9)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().y() , "Length")<< " "
		 << G4std::setw( 9)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().z() , "Length") << " "
		 << G4std::setw( 9)
		 << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy() , "Energy")<< " "
		 << G4std::setw(18)
		 << (*fSecondary)[lp1]->GetDefinition()
 	                              ->GetParticleName();
	  G4cout << G4endl;
	}
              
	G4cout << "    :-----------------------------"
	       << "----------------------------------"
	       << "-- EndOf2ndaries Info ---------------"
	       << G4endl;
      }
    }
    
  }
  G4cout.precision(prec);
}

////////////////////////////////////////////
void G4SteppingVerbose::DPSLStarted()
////////////////////////////////////////////
{
  CopyState();

  if( verboseLevel > 5 ){
    G4cout << G4endl;
    G4cout << " >>DefinePhysicalStepLength (List of proposed StepLengths): "
      << G4endl;
  }
}
//////////////////////////////////////////////
void G4SteppingVerbose::DPSLUserLimit()
//////////////////////////////////////////////
{
  CopyState();

  if( verboseLevel > 5 ){
    G4cout << G4endl << G4endl;
    G4cout << "=== Defined Physical Step Length (DPSL)" << G4endl;
    G4cout << "    ++ProposedStep(UserLimit) = " 
      << G4std::setw( 9) << physIntLength
	<< " : ProcName = User defined maximum allowed Step"
	  << G4endl;
  }
}
/////////////////////////////////////////////
void G4SteppingVerbose::DPSLPostStep()
/////////////////////////////////////////////
{
  CopyState();

  if( verboseLevel > 5 ){
    G4cout << "    ++ProposedStep(PostStep ) = " 
      << G4std::setw( 9) << physIntLength
	<< " : ProcName = "
	  << fCurrentProcess->GetProcessName() 
            << " (";
    if(fCondition==ExclusivelyForced){
      G4cout << "ExclusivelyForced)" << G4endl;
    }
    else if(fCondition==Conditionally){
      G4cout << "Conditionally)" << G4endl;
    }
    else if(fCondition==Forced){
      G4cout << "Forced)" << G4endl;
    }
    else{
      G4cout << "No ForceCondition)" << G4endl;
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
	   << G4std::setw( 9) << G4BestUnit(physIntLength , "Length")
	   << " : ProcName = "
	   << fCurrentProcess->GetProcessName() 
	   << " (";
    if(fGPILSelection==CandidateForSelection){
      G4cout << "CandidateForSelection)" << G4endl;
    }
    else if(fGPILSelection==NotCandidateForSelection){
      G4cout << "NotCandidateForSelection)" << G4endl;
    }
    else{
      G4cout << "???)" << G4endl;
    }
  }
}


////////////////////////////////////////////////
void G4SteppingVerbose::TrackingStarted()
////////////////////////////////////////////////
{

  CopyState();
G4int prec = G4cout.precision(3);
  if( verboseLevel > 0 ){

#ifdef G4_USE_G4BESTUNIT_FOR_VERBOSE
    G4cout << G4std::setw( 5) << "Step#"  << " "
           << G4std::setw( 8) << "X"      << "     "
	   << G4std::setw( 8) << "Y"      << "     "  
	   << G4std::setw( 8) << "Z"      << "     "
	   << G4std::setw( 9) << "KineE"  << "     "
	   << G4std::setw( 8) << "dE"     << "     "  
	   << G4std::setw(12) << "StepLeng"   << " "  
	   << G4std::setw(12) << "TrackLeng"  << " "
	   << G4std::setw(12) << "NextVolume" << " "
	   << G4std::setw( 8) << "ProcName"   << G4endl;	     
#else
    G4cout << G4std::setw( 5) << "Step#"      << " "
	   << G4std::setw( 8) << "X(mm)"      << " "
	   << G4std::setw( 8) << "Y(mm)"      << " "  
	   << G4std::setw( 8) << "Z(mm)"      << " "
	   << G4std::setw( 9) << "KinE(MeV)"  << " "
	   << G4std::setw( 8) << "dE(MeV)"    << " "  
	   << G4std::setw( 8) << "StepLeng"   << " "  
	   << G4std::setw( 9) << "TrackLeng"  << " "
	   << G4std::setw(11) << "NextVolume" << " "
	   << G4std::setw( 8) << "ProcName"   << G4endl;	     
#endif

    G4cout << G4std::setw( 5) << fTrack->GetCurrentStepNumber() << " "
	   << G4std::setw( 8) << G4BestUnit(fTrack->GetPosition().x(),"Length")<< " "
	   << G4std::setw( 8) << G4BestUnit(fTrack->GetPosition().y(),"Length") << " "
	   << G4std::setw( 8) << G4BestUnit(fTrack->GetPosition().z(),"Length")<< " "
	   << G4std::setw( 9) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")<< " "
	   << G4std::setw( 8) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy") << " "
	   << G4std::setw( 8) << G4BestUnit(fStep->GetStepLength(),"Length")<< " "
	   << G4std::setw( 9) << G4BestUnit(fTrack->GetTrackLength(),"Length") << " ";

    if(fTrack->GetNextVolume()){
      G4cout << G4std::setw(11) << fTrack->GetNextVolume()->GetName() << " ";
    } else {
      G4cout << G4std::setw(11) << "OutOfWorld" << " ";
    }
    G4cout << "initStep" << G4endl;
  }
  G4cout.precision(prec);
}
//////////////////////////////////////////////////////
void G4SteppingVerbose::AlongStepDoItOneByOne()
//////////////////////////////////////////////////////
{ 
  CopyState();
    if(verboseLevel >= 4){ 
        G4cout << G4endl;
        G4cout << " >>AlongStepDoIt (process by process): "
             << "   Process Name = " 
             << fCurrentProcess->GetProcessName() << G4endl;

        ShowStep();
        G4cout << "          "
	       << "!Note! Safety of PostStep is only valid "
	       << "after all DoIt invocations."
	       << G4endl; 

        VerboseParticleChange();    
        G4cout << G4endl;

        G4cout << "    ++List of secondaries generated " 
	       << "(x,y,z,kE,t,PID):"
	       << "  No. of secodaries = " 
	       << fN2ndariesAlongStepDoIt << G4endl;

        if(fN2ndariesAlongStepDoIt>0){
           for(G4int lp1=(*fSecondary).size()-fN2ndariesAlongStepDoIt; 
                     lp1<(*fSecondary).size(); lp1++){
               G4cout << "      "
                    << G4std::setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetPosition().x() , "Length")<< " "
                    << G4std::setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetPosition().y() , "Length")<< " "
                    << G4std::setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetPosition().z() , "Length")<< " "
                    << G4std::setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy() , "Energy")<< " "
                    << G4std::setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetGlobalTime() , "Time")<< " "
                    << G4std::setw(18)
                    << (*fSecondary)[lp1]->GetDefinition()
                                         ->GetParticleName();
               G4cout << G4endl;
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
        G4cout << G4endl;
        G4cout << " >>PostStepDoIt (process by process): "
             << "   Process Name = " 
             << fCurrentProcess->GetProcessName() << G4endl;

        ShowStep();
        G4cout << G4endl;
        VerboseParticleChange();    
        G4cout << G4endl;
         
        G4cout << "    ++List of secondaries generated " 
             << "(x,y,z,kE,t,PID):"
             << "  No. of secodaries = " 
             << fN2ndariesPostStepDoIt << G4endl;

        if(fN2ndariesPostStepDoIt>0){
           for(G4int lp1=(*fSecondary).size()-fN2ndariesPostStepDoIt; 
                     lp1<(*fSecondary).size(); lp1++){
               G4cout << "      "
                    << G4std::setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetPosition().x() , "Length")<< " "
                    << G4std::setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(), "Length") << " "
                    << G4std::setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(), "Length") << " "
                    << G4std::setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(), "Energy") << " "
                    << G4std::setw( 9)
                    << G4BestUnit((*fSecondary)[lp1]->GetGlobalTime(), "Time") << " "
                    << G4std::setw(18)
                    << (*fSecondary)[lp1]->GetDefinition()
                                         ->GetParticleName();
               G4cout << G4endl;
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
  G4cout << G4endl;
  G4cout << "    ++G4Track Information " << G4endl;
  G4int prec = G4cout.precision(3);


  G4cout << "      -----------------------------------------------" 
       << G4endl;
  G4cout << "        G4Track Information  " << G4std::setw(20) << G4endl;
  G4cout << "      -----------------------------------------------" 
       << G4endl;

  G4cout << "        Step number         : " 
       << G4std::setw(20) << fTrack->GetCurrentStepNumber()
       << G4endl; 
#ifdef G4_USE_G4BESTUNIT_FOR_VERBOSE
  G4cout << "        Position - x        : " 
       << G4std::setw(20) << G4BestUnit(fTrack->GetPosition().x(), "Length")
       << G4endl; 
  G4cout << "        Position - y        : " 
       << G4std::setw(20) << G4BestUnit(fTrack->GetPosition().y(), "Length")
       << G4endl; 
  G4cout << "        Position - z        : " 
       << G4std::setw(20) << G4BestUnit(fTrack->GetPosition().z(), "Length")
       << G4endl;
  G4cout << "        Global Time         : " 
       << G4std::setw(20) << G4BestUnit(fTrack->GetGlobalTime(), "Time")
       << G4endl;
  G4cout << "        Local Time          : " 
       << G4std::setw(20) << G4BestUnit(fTrack->GetLocalTime(), "Time")
       << G4endl;
#else
  G4cout << "        Position - x (mm)   : " 
       << G4std::setw(20) << fTrack->GetPosition().x() /mm
       << G4endl; 
  G4cout << "        Position - y (mm)   : " 
       << G4std::setw(20) << fTrack->GetPosition().y() /mm
       << G4endl; 
  G4cout << "        Position - z (mm)   : " 
       << G4std::setw(20) << fTrack->GetPosition().z() /mm
       << G4endl;
  G4cout << "        Global Time (ns)    : " 
       << G4std::setw(20) << fTrack->GetGlobalTime() /ns
       << G4endl;
  G4cout << "        Local Time (ns)     : " 
       << G4std::setw(20) << fTrack->GetLocalTime() /ns
       << G4endl;
#endif
  G4cout << "        Momentum Direct - x : " 
       << G4std::setw(20) << fTrack->GetMomentumDirection().x()
       << G4endl;
  G4cout << "        Momentum Direct - y : " 
       << G4std::setw(20) << fTrack->GetMomentumDirection().y()
       << G4endl;
  G4cout << "        Momentum Direct - z : " 
       << G4std::setw(20) << fTrack->GetMomentumDirection().z()
       << G4endl;
#ifdef G4_USE_G4BESTUNIT_FOR_VERBOSE
  G4cout << "        Kinetic Energy      : " 
#else
  G4cout << "        Kinetic Energy (MeV): " 
#endif
       << G4std::setw(20) << G4BestUnit(fTrack->GetKineticEnergy(), "Energy")
       << G4endl;
  G4cout << "        Polarization - x    : " 
       << G4std::setw(20) << fTrack->GetPolarization().x()
       << G4endl;
  G4cout << "        Polarization - y    : " 
       << G4std::setw(20) << fTrack->GetPolarization().y()
       << G4endl;
  G4cout << "        Polarization - z    : " 
       << G4std::setw(20) << fTrack->GetPolarization().z()
       << G4endl;
  G4cout << "        Track Length        : " 
       << G4std::setw(20) << G4BestUnit(fTrack->GetTrackLength(), "Length")
       << G4endl;
  G4cout << "        Track ID #          : " 
       << G4std::setw(20) << fTrack->GetTrackID()
       << G4endl;
  G4cout << "        Parent Track ID #   : " 
       << G4std::setw(20) << fTrack->GetParentID()
       << G4endl;
  G4cout << "        Next Volume         : " 
       << G4std::setw(20);
       if( fTrack->GetNextVolume() != 0 ) { 
         G4cout << fTrack->GetNextVolume()->GetName() << " ";
       } else {
         G4cout << "OutOfWorld" << " ";
       }
       G4cout << G4endl;
  G4cout << "        Track Status        : " 
       << G4std::setw(20);
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
       G4cout << G4endl;
#ifdef G4_USE_G4BESTUNIT_FOR_VERBOSE
  G4cout << "        Vertex - x          : " 
       << G4std::setw(20) << G4BestUnit(fTrack->GetVertexPosition().x(),"Length")
       << G4endl; 
  G4cout << "        Vertex - y          : " 
       << G4std::setw(20) << G4BestUnit(fTrack->GetVertexPosition().y(),"Length")
       << G4endl; 
  G4cout << "        Vertex - z          : " 
       << G4std::setw(20) << G4BestUnit(fTrack->GetVertexPosition().z(),"Length")
       << G4endl;
#else
  G4cout << "        Vertex - x (mm)     : " 
       << G4std::setw(20) << fTrack->GetVertexPosition().x()/mm
       << G4endl; 
  G4cout << "        Vertex - y (mm)     : " 
       << G4std::setw(20) << fTrack->GetVertexPosition().y()/mm
       << G4endl; 
  G4cout << "        Vertex - z (mm)     : " 
       << G4std::setw(20) << fTrack->GetVertexPosition().z()/mm
       << G4endl;
#endif
  G4cout << "        Vertex - Px (MomDir): " 
       << G4std::setw(20) << fTrack->GetVertexMomentumDirection().x()
       << G4endl;
  G4cout << "        Vertex - Py (MomDir): " 
       << G4std::setw(20) << fTrack->GetVertexMomentumDirection().y()
       << G4endl;
  G4cout << "        Vertex - Pz (MomDir): " 
       << G4std::setw(20) << fTrack->GetVertexMomentumDirection().z()
       << G4endl;
#ifdef G4_USE_G4BESTUNIT_FOR_VERBOSE
  G4cout << "        Vertex - KineE      : " 
#else
  G4cout << "        Vertex - KineE (MeV): " 
#endif
       << G4std::setw(20) << G4BestUnit(fTrack->GetVertexKineticEnergy(),"Energy")
       << G4endl;
  
  G4cout << "        Creator Process     : " 
       << G4std::setw(20);
  if( fTrack->GetCreatorProcess() == NULL){
    G4cout << " Event Generator" << G4endl;
  } else {
    G4cout << fTrack->GetCreatorProcess()->GetProcessName() << G4endl;
  }

  G4cout << "      -----------------------------------------------" 
       << G4endl;
       
 G4cout.precision(prec);      
}


///////////////////////////////////////////////
void G4SteppingVerbose::VerboseParticleChange()
///////////////////////////////////////////////
{

// Show header
  G4cout << G4endl;
  G4cout << "    ++G4ParticleChange Information " << G4endl;
  fParticleChange->DumpInfo();
}
/////////////////////////////////////////
void G4SteppingVerbose::ShowStep() const
////////////////////////////////////////
{

// Show header
   G4cout << G4endl;
   G4cout << "    ++G4Step Information " << G4endl;
   G4cout.precision(3);

// Show G4Step specific information
   G4cout << "      Address of G4Track    : " << fStep->GetTrack() << G4endl;
   G4cout << "      Step Length (mm)      : " << fStep->GetTrack()->GetStepLength() << G4endl;
   G4cout << "      Energy Deposit (MeV)  : " << fStep->GetTotalEnergyDeposit() << G4endl;

// Show G4StepPoint specific information
   G4cout << "      -------------------------------------------------------" 
        << "----------------" <<  G4endl;
   G4cout << "        StepPoint Information  " << G4std::setw(20) << "PreStep" 
                                             << G4std::setw(20) << "PostStep" << G4endl;
   G4cout << "      -------------------------------------------------------" 
        << "----------------" <<  G4endl;
   G4cout << "         Position - x (mm)   : " 
        << G4std::setw(20) << fStep->GetPreStepPoint()->GetPosition().x() 
        << G4std::setw(20) << fStep->GetPostStepPoint()->GetPosition().x() << G4endl;
   G4cout << "         Position - y (mm)   : " 
        << G4std::setw(20) << fStep->GetPreStepPoint()->GetPosition().y() 
        << G4std::setw(20) << fStep->GetPostStepPoint()->GetPosition().y() << G4endl;
   G4cout << "         Position - z (mm)   : " 
        << G4std::setw(20) << fStep->GetPreStepPoint()->GetPosition().z() 
        << G4std::setw(20) << fStep->GetPostStepPoint()->GetPosition().z() << G4endl;
   G4cout << "         Global Time (ns)    : " 
        << G4std::setw(20) << fStep->GetPreStepPoint()->GetGlobalTime()
        << G4std::setw(20) << fStep->GetPostStepPoint()->GetGlobalTime() << G4endl;
   G4cout << "         Local Time (ns)     : " 
        << G4std::setw(20) << fStep->GetPreStepPoint()->GetLocalTime() 
        << G4std::setw(20) << fStep->GetPostStepPoint()->GetLocalTime() << G4endl;
   G4cout << "         Proper Time (ns)    : " 
        << G4std::setw(20) << fStep->GetPreStepPoint()->GetProperTime()
        << G4std::setw(20) << fStep->GetPostStepPoint()->GetProperTime() << G4endl;
   G4cout << "         Momentum Direct - x : " 
        << G4std::setw(20) << fStep->GetPreStepPoint()->GetMomentumDirection().x()
        << G4std::setw(20) << fStep->GetPostStepPoint()->GetMomentumDirection().x() << G4endl;
   G4cout << "         Momentum Direct - y : " 
        << G4std::setw(20) << fStep->GetPreStepPoint()->GetMomentumDirection().y()
        << G4std::setw(20) << fStep->GetPostStepPoint()->GetMomentumDirection().y() << G4endl;
   G4cout << "         Momentum Direct - z : " 
        << G4std::setw(20) << fStep->GetPreStepPoint()->GetMomentumDirection().z()
        << G4std::setw(20) << fStep->GetPostStepPoint()->GetMomentumDirection().z() << G4endl;
   G4cout << "         Momentum - x (MeV/c): " 
        << G4std::setw(20) << fStep->GetPreStepPoint()->GetMomentum().x()
        << G4std::setw(20) << fStep->GetPostStepPoint()->GetMomentum().x() << G4endl;
   G4cout << "         Momentum - y (MeV/c): " 
        << G4std::setw(20) << fStep->GetPreStepPoint()->GetMomentum().y()
        << G4std::setw(20) << fStep->GetPostStepPoint()->GetMomentum().y() << G4endl;
   G4cout << "         Momentum - z (MeV/c): " 
        << G4std::setw(20) << fStep->GetPreStepPoint()->GetMomentum().z()
        << G4std::setw(20) << fStep->GetPostStepPoint()->GetMomentum().z() << G4endl;
   G4cout << "         Total Energy (MeV)  : " 
        << G4std::setw(20) << fStep->GetPreStepPoint()->GetTotalEnergy()
        << G4std::setw(20) << fStep->GetPostStepPoint()->GetTotalEnergy() << G4endl;
   G4cout << "         Kinetic Energy (MeV): " 
        << G4std::setw(20) << fStep->GetPreStepPoint()->GetKineticEnergy()
        << G4std::setw(20) << fStep->GetPostStepPoint()->GetKineticEnergy() << G4endl;
   G4cout << "         Velocity (mm/ns)    : " 
        << G4std::setw(20) << fStep->GetPreStepPoint()->GetVelocity()
        << G4std::setw(20) << fStep->GetPostStepPoint()->GetVelocity() << G4endl;
   G4cout << "         Volume Name         : " 
        << G4std::setw(20) << fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()
        << G4std::setw(20) << fStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() << G4endl;
   G4cout << "         Safety (mm)         : " 
        << G4std::setw(20) << fStep->GetPreStepPoint()->GetSafety()
        << G4std::setw(20) << fStep->GetPostStepPoint()->GetSafety() << G4endl;
   G4cout << "         Polarization - x    : " 
        << G4std::setw(20) << fStep->GetPreStepPoint()->GetPolarization().x()
        << G4std::setw(20) << fStep->GetPostStepPoint()->GetPolarization().x() << G4endl;
   G4cout << "         Polarization - y    : " 
        << G4std::setw(20) << fStep->GetPreStepPoint()->GetPolarization().y()
        << G4std::setw(20) << fStep->GetPostStepPoint()->GetPolarization().y() << G4endl;
   G4cout << "         Polarization - Z    : " 
        << G4std::setw(20) << fStep->GetPreStepPoint()->GetPolarization().z()
        << G4std::setw(20) << fStep->GetPostStepPoint()->GetPolarization().z() << G4endl;
   G4cout << "         Weight              : " 
        << G4std::setw(20) << fStep->GetPreStepPoint()->GetWeight()
        << G4std::setw(20) << fStep->GetPostStepPoint()->GetWeight() << G4endl;
   G4cout << "         Step Status         : " ;
        G4StepStatus  tStepStatus = fStep->GetPreStepPoint()->GetStepStatus();
        if( tStepStatus == fGeomBoundary ){
           G4cout << G4std::setw(20) << "Geom Limit";
        } else if ( tStepStatus == fAlongStepDoItProc ){
          G4cout << G4std::setw(20) << "AlongStep Proc.";
        } else if ( tStepStatus == fPostStepDoItProc ){
           G4cout << G4std::setw(20) << "PostStep Proc";
        } else if ( tStepStatus == fAtRestDoItProc ){
           G4cout << G4std::setw(20) << "AtRest Proc";
        } else if ( tStepStatus == fUndefined ){
           G4cout << G4std::setw(20) << "Undefined";
        }

        tStepStatus = fStep->GetPostStepPoint()->GetStepStatus();
        if( tStepStatus == fGeomBoundary ){
           G4cout << G4std::setw(20) << "Geom Limit";
        } else if ( tStepStatus == fAlongStepDoItProc ){
           G4cout << G4std::setw(20) << "AlongStep Proc.";
        } else if ( tStepStatus == fPostStepDoItProc ){
           G4cout << G4std::setw(20) << "PostStep Proc";
        } else if ( tStepStatus == fAtRestDoItProc ){
           G4cout << G4std::setw(20) << "AtRest Proc";
        } else if ( tStepStatus == fUndefined ){
           G4cout << G4std::setw(20) << "Undefined";
        }

        G4cout << G4endl;
        G4cout << "         Process defined Step: " ;
        if( fStep->GetPreStepPoint()->GetProcessDefinedStep() == NULL ){
 	   G4cout << G4std::setw(20) << "Undefined";
        } else {
  	   G4cout << G4std::setw(20) << fStep->GetPreStepPoint()->GetProcessDefinedStep()
                                             ->GetProcessName();
        }
        if( fStep->GetPostStepPoint()->GetProcessDefinedStep() == NULL){
  	   G4cout << G4std::setw(20) << "Undefined";
        } else {
 	   G4cout << G4std::setw(20) << fStep->GetPostStepPoint()->GetProcessDefinedStep()
                                              ->GetProcessName(); 
        }

   G4cout << G4endl;
   G4cout << "      -------------------------------------------------------" 
        << "----------------" << G4endl;
}


