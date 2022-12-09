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
// G4SteppingVerbose class implementation
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
// --------------------------------------------------------------------

#include "G4SteppingVerbose.hh"
#include "G4SteppingManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VSensitiveDetector.hh"    // Include from 'hits/digi'
#include "G4StepStatus.hh"            // Include from 'tracking'

// #define G4_USE_G4BESTUNIT_FOR_VERBOSE 1

#ifdef G4_USE_G4BESTUNIT_FOR_VERBOSE
#include "G4UnitsTable.hh"
#else
#define G4BestUnit(a,b) a
#endif

G4int G4SteppingVerbose::useBestUnitPrecision = -1;

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
  if(Silent==1) { return; }

  G4VProcess* ptProcManager;
  CopyState();

  if(verboseLevel >= 3)
  {
    G4int npt=0;
    G4cout << " **List of AtRestDoIt invoked:" << G4endl;
    for(std::size_t np=0; np<MAXofAtRestLoops; ++np)
    {
      std::size_t npGPIL = MAXofAtRestLoops-np-1;
      if( (*fSelectedAtRestDoItVector)[npGPIL] == G4ForceCondition::Forced )
      {
        ++npt;                
        ptProcManager = (*fAtRestDoItVector)[(G4int)np];
        G4cout << "   # " << npt << " : " 
               << ptProcManager->GetProcessName() 
               << " (Forced)" << G4endl;
      }
      else if ( (*fSelectedAtRestDoItVector)[npGPIL] == G4ForceCondition::NotForced )
      {
        ++npt;                
        ptProcManager = (*fAtRestDoItVector)[(G4int)np];
        G4cout << "   # " << npt << " : "  << ptProcManager->GetProcessName()
               << G4endl;
      }
    }
     
    G4cout << "   Generated secondries # : " << fN2ndariesAtRestDoIt << G4endl;
     
    if( fN2ndariesAtRestDoIt > 0 )
    {
      G4cout << "   -- List of secondaries generated : "
             << "(x,y,z,kE,t,PID) --" << G4endl; 
      for( std::size_t lp1=(*fSecondary).size()-fN2ndariesAtRestDoIt;
                       lp1<(*fSecondary).size(); ++lp1)
      {
        G4cout << "      "
               << std::setw( 9)
               << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length")
               << " " << std::setw( 9)
               << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length")
               << " " << std::setw( 9)
               << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length")
               << " " << std::setw( 9)
               << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy")
               << " " << std::setw( 9)
               << G4BestUnit((*fSecondary)[lp1]->GetGlobalTime(),"Time")
               << " " << std::setw(18)
               << (*fSecondary)[lp1]->GetDefinition()->GetParticleName()
               << G4endl;
      }
    }
  }
   
  if( verboseLevel >= 4 )
  { 
    ShowStep();
    G4cout << G4endl;
  }
}

/////////////////////////////////////////////////////
void G4SteppingVerbose::AlongStepDoItAllDone()
/////////////////////////////////////////////////////
{
  if(Silent==1){ return; }

  G4VProcess* ptProcManager;

  CopyState();

  if(verboseLevel >= 3)
  {
    G4cout << G4endl;
    G4cout << " >>AlongStepDoIt (after all invocations):" << G4endl;
    G4cout << "    ++List of invoked processes " << G4endl;

    for(std::size_t ci=0; ci<MAXofAlongStepLoops; ++ci)
    {
      ptProcManager = (*fAlongStepDoItVector)((G4int)ci);
      G4cout << "      " << ci+1 << ") ";
      if(ptProcManager != nullptr)
      {
        G4cout << ptProcManager->GetProcessName() << G4endl;
      }
    }

    ShowStep();
    G4cout << G4endl;
    G4cout << "    ++List of secondaries generated " 
           << "(x,y,z,kE,t,PID):"
           << "  No. of secondaries = " 
           << (*fSecondary).size() << G4endl;

    if((*fSecondary).size()>0)
    {
      for(std::size_t lp1=0; lp1<(*fSecondary).size(); ++lp1)
      {
        G4cout << "      "
               << std::setw( 9)
               << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length")
               << " " << std::setw( 9)
               << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length")
               << " " << std::setw( 9)
               << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length")
               << " " << std::setw( 9)
               << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy")
               << " " << std::setw( 9)
               << G4BestUnit((*fSecondary)[lp1]->GetGlobalTime(),"Time")
               << " " << std::setw(18)
               << (*fSecondary)[lp1]->GetDefinition()->GetParticleName()
               << G4endl;
      }
    }
  }
}

////////////////////////////////////////////////////
void G4SteppingVerbose::PostStepDoItAllDone()
////////////////////////////////////////////////////
{
  if(Silent==1) { return; }

  G4VProcess* ptProcManager;

  CopyState();

  if( (fStepStatus == fPostStepDoItProc) |
      (fCondition  == Forced)            |
      (fCondition  == Conditionally)     |
      (fCondition  == ExclusivelyForced) |
      (fCondition  == StronglyForced) )
  {
    if(verboseLevel >= 3)
    {
      G4int npt=0;
      G4cout << G4endl;
      G4cout << " **PostStepDoIt (after all invocations):" << G4endl;
      G4cout << "    ++List of invoked processes " << G4endl;

      for(std::size_t np=0; np<MAXofPostStepLoops; ++np)
      {
        std::size_t npGPIL = MAXofPostStepLoops-np-1;
        if( (*fSelectedPostStepDoItVector)[npGPIL] == G4ForceCondition::Forced )
        {
          ++npt;
          ptProcManager = (*fPostStepDoItVector)[(G4int)np];
          G4cout << "      " << npt << ") "
                 << ptProcManager->GetProcessName()
                 << " (Forced)" << G4endl;
        }
        else if ( (*fSelectedPostStepDoItVector)[npGPIL] == G4ForceCondition::NotForced )
        {
          ++npt;
          ptProcManager = (*fPostStepDoItVector)[(G4int)np];
          G4cout << "      " << npt << ") " << ptProcManager->GetProcessName()
                 << G4endl;
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

      if((*fSecondary).size()>0)
      {
        for(std::size_t lp1=0; lp1<(*fSecondary).size(); ++lp1)
        {
          G4cout << "      "
                 << std::setw( 9)
                 << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(), "Length")
                 << " " << std::setw( 9)
                 << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(), "Length")
                 << " " << std::setw( 9)
                 << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(), "Length")
                 << " " << std::setw( 9)
                 << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(), "Energy")
                 << " " << std::setw( 9)
                 << G4BestUnit((*fSecondary)[lp1]->GetGlobalTime(), "Time")
                 << " " << std::setw(18)
                 << (*fSecondary)[lp1]->GetDefinition()->GetParticleName()
                 << G4endl;
        }
      }
    }
  }
}

/////////////////////////////////////////
void G4SteppingVerbose::StepInfo()
/////////////////////////////////////////
{
  if(Silent==1) { return; }
  if(SilentStepInfo==1) { return; }

  CopyState();
  G4cout.precision(16); 
  G4long prec = G4cout.precision(3);

  if( verboseLevel >= 1 )
  {
    if( verboseLevel >= 4 ) VerboseTrack();
    if( verboseLevel >= 3 )
    {
      G4cout << G4endl;
#ifdef G4_USE_G4BESTUNIT_FOR_VERBOSE      
      G4cout << std::setw( 5) << "#Step#" << " "
             << std::setw( 8) << "X"      << "     "
             << std::setw( 8) << "Y"      << "     "
             << std::setw( 8) << "Z"      << "     "
             << std::setw( 9) << "KineE"  << "     "
             << std::setw( 8) << "dE"     << "     "
             << std::setw(12) << "StepLeng"   << " "
             << std::setw(12) << "TrackLeng"  << " "
             << std::setw(12) << "NextVolume" << " "
             << std::setw( 8) << "ProcName"   << G4endl;
#else
      G4cout << std::setw( 5) << "#Step#"      << " "
             << std::setw( 8) << "X(mm)"      << " "
             << std::setw( 8) << "Y(mm)"      << " "  
             << std::setw( 8) << "Z(mm)"      << " "
             << std::setw( 9) << "KinE(MeV)"  << " "
             << std::setw( 8) << "dE(MeV)"    << " "  
             << std::setw( 8) << "StepLeng"   << " "
             << std::setw( 9) << "TrackLeng"  << " "
             << std::setw(11) << "NextVolume" << " "
             << std::setw( 8) << "ProcName"   << G4endl;
#endif     
    }
    G4cout << std::setw( 5) << fTrack->GetCurrentStepNumber() << " "
           << std::setw( 8)
           << G4BestUnit(fTrack->GetPosition().x(), "Length") << " "
           << std::setw( 8)
           << G4BestUnit(fTrack->GetPosition().y(), "Length") << " "
           << std::setw( 8)
           << G4BestUnit(fTrack->GetPosition().z(), "Length") << " "
           << std::setw( 9)
           << G4BestUnit(fTrack->GetKineticEnergy(), "Energy") << " "
           << std::setw( 8)
           << G4BestUnit(fStep->GetTotalEnergyDeposit(), "Energy") << " "
           << std::setw( 8)
           << G4BestUnit(fStep->GetStepLength(), "Length") << " "
           << std::setw( 9)
           << G4BestUnit(fTrack->GetTrackLength(), "Length") << " ";

    if( fTrack->GetNextVolume() != 0 )
    { 
      G4cout << std::setw(11) << fTrack->GetNextVolume()->GetName() << " ";
    }
    else
    {
      G4cout << std::setw(11) << "OutOfWorld" << " ";
    }
    if(fStep->GetPostStepPoint()->GetProcessDefinedStep() != 0)
    {
      G4cout << fStep->GetPostStepPoint()
                     ->GetProcessDefinedStep()->GetProcessName();
    }
    else
    {
      G4cout << "User Limit";
    }
    G4cout << G4endl;
    if( verboseLevel == 2 )
    {
      G4int tN2ndariesTot = fN2ndariesAtRestDoIt + fN2ndariesAlongStepDoIt
                          + fN2ndariesPostStepDoIt;
      if(tN2ndariesTot>0)
      {
        G4cout << "    :----- List of 2ndaries - "
               << "#SpawnInStep=" << std::setw(3) << tN2ndariesTot 
               << "(Rest=" << std::setw(2) << fN2ndariesAtRestDoIt
               << ",Along=" << std::setw(2) << fN2ndariesAlongStepDoIt
               << ",Post="  << std::setw(2) << fN2ndariesPostStepDoIt
               << "), "
               << "#SpawnTotal=" << std::setw(3) << (*fSecondary).size()
               << " ---------------"
               << G4endl;

        for(std::size_t lp1=(*fSecondary).size()-tN2ndariesTot;
                        lp1<(*fSecondary).size(); ++lp1)
        {
          G4cout << "    : "
                 << std::setw( 9)
                 << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(), "Length")
                 << " " << std::setw( 9)
                 << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(), "Length")
                 << " " << std::setw( 9)
                 << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(), "Length")
                 << " " << std::setw( 9)
                 << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(), "Energy")
                 << " " << std::setw(18)
                 << (*fSecondary)[lp1]->GetDefinition()->GetParticleName()
                 << G4endl;
        }
        G4cout << "    :-----------------------------"
               << "----------------------------------"
               << "-- EndOf2ndaries Info ---------------"  << G4endl;
      }
    }
  }
  G4cout.precision(prec);
}

////////////////////////////////////////////
void G4SteppingVerbose::DPSLStarted()
////////////////////////////////////////////
{
  if(Silent==1){ return; }
  CopyState();

  if( verboseLevel > 5 )
  {
    G4cout << G4endl
           << " >>DefinePhysicalStepLength (List of proposed StepLengths): "
           << G4endl;
  }
}

//////////////////////////////////////////////
void G4SteppingVerbose::DPSLUserLimit()
//////////////////////////////////////////////
{
  if(Silent==1){ return; }
  CopyState();

  if( verboseLevel > 5 )
  {
    G4cout << G4endl << G4endl;
    G4cout << "=== Defined Physical Step Length (DPSL)" << G4endl;
    G4cout << "    ++ProposedStep(UserLimit) = "
           << std::setw( 9) << physIntLength
           << " : ProcName = User defined maximum allowed Step" << G4endl;
  }
}

/////////////////////////////////////////////
void G4SteppingVerbose::DPSLPostStep()
/////////////////////////////////////////////
{
  if(Silent==1){ return; }
  CopyState();

  if( verboseLevel > 5 )
  {
    G4cout << "    ++ProposedStep(PostStep ) = "
           << std::setw( 9) << physIntLength
           << " : ProcName = " << fCurrentProcess->GetProcessName() << " (";
    if(fCondition==ExclusivelyForced)
    {
      G4cout << "ExclusivelyForced)" << G4endl;
    }
    else if(fCondition==StronglyForced)
    {
      G4cout << "StronglyForced)" << G4endl;
    }
    else if(fCondition==Conditionally)
    {
      G4cout << "Conditionally)" << G4endl;
    }
    else if(fCondition==Forced)
    {
      G4cout << "Forced)" << G4endl;
    }
    else
    {
      G4cout << "No ForceCondition)" << G4endl;
    }
  }
}

/////////////////////////////////////////////
void G4SteppingVerbose::DPSLAlongStep()
/////////////////////////////////////////////
{
  if(Silent==1){ return; }
  CopyState();

  if( verboseLevel > 5 )
  {
    G4cout << "    ++ProposedStep(AlongStep) = " 
           << std::setw( 9) << G4BestUnit(physIntLength , "Length")
           << " : ProcName = "
           << fCurrentProcess->GetProcessName() 
           << " (";
    if(fGPILSelection==CandidateForSelection)
    {
      G4cout << "CandidateForSelection)" << G4endl;
    }
    else if(fGPILSelection==NotCandidateForSelection)
    {
      G4cout << "NotCandidateForSelection)" << G4endl;
    }
    else
    {
      G4cout << "?!?)" << G4endl;
    }
  }
}

////////////////////////////////////////////////
void G4SteppingVerbose::TrackingStarted()
////////////////////////////////////////////////
{
  if(Silent==1){ return; }

  CopyState();

  G4long prec = G4cout.precision(3);
  if( verboseLevel > 0 )
  {
#ifdef G4_USE_G4BESTUNIT_FOR_VERBOSE
    G4cout << std::setw( 5) << "Step#"  << " "
           << std::setw( 8) << "X"      << "     "
           << std::setw( 8) << "Y"      << "     "  
           << std::setw( 8) << "Z"      << "     "
           << std::setw( 9) << "KineE"  << "     "
           << std::setw( 8) << "dE"     << "     "  
           << std::setw(12) << "StepLeng"   << " "  
           << std::setw(12) << "TrackLeng"  << " "
           << std::setw(12) << "NextVolume" << " "
           << std::setw( 8) << "ProcName"   << G4endl;     
#else
    G4cout << std::setw( 5) << "Step#"      << " "
           << std::setw( 8) << "X(mm)"      << " "
           << std::setw( 8) << "Y(mm)"      << " "  
           << std::setw( 8) << "Z(mm)"      << " "
           << std::setw( 9) << "KinE(MeV)"  << " "
           << std::setw( 8) << "dE(MeV)"    << " "  
           << std::setw( 8) << "StepLeng"   << " "  
           << std::setw( 9) << "TrackLeng"  << " "
           << std::setw(11) << "NextVolume" << " "
           << std::setw( 8) << "ProcName"   << G4endl;     
#endif

    G4cout << std::setw( 5) << fTrack->GetCurrentStepNumber() << " "
           << std::setw( 8)
           << G4BestUnit(fTrack->GetPosition().x(),"Length") << " "
           << std::setw( 8)
           << G4BestUnit(fTrack->GetPosition().y(),"Length") << " "
           << std::setw( 8)
           << G4BestUnit(fTrack->GetPosition().z(),"Length") << " "
           << std::setw( 9)
           << G4BestUnit(fTrack->GetKineticEnergy(),"Energy") << " "
           << std::setw( 8)
           << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy") << " "
           << std::setw( 8)
           << G4BestUnit(fStep->GetStepLength(),"Length") << " "
           << std::setw( 9)
           << G4BestUnit(fTrack->GetTrackLength(),"Length") << " ";

    if(fTrack->GetNextVolume())
    {
      G4cout << std::setw(11) << fTrack->GetNextVolume()->GetName() << " ";
    }
    else
    {
      G4cout << std::setw(11) << "OutOfWorld" << " ";
    }
    G4cout << "initStep" << G4endl;
  }
  G4cout.precision(prec);
}

//////////////////////////////////////////////////////
void G4SteppingVerbose::AlongStepDoItOneByOne()
//////////////////////////////////////////////////////
{ 
  if(Silent==1){ return; }

  CopyState();

  if(verboseLevel >= 4)
  { 
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

    if(fN2ndariesAlongStepDoIt>0)
    {
      for(std::size_t lp1=(*fSecondary).size()-fN2ndariesAlongStepDoIt;
                      lp1<(*fSecondary).size(); ++lp1)
      {
         G4cout  << "      "
                 << std::setw( 9)
                 << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(), "Length")
                 << " " << std::setw( 9)
                 << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(), "Length")
                 << " " << std::setw( 9)
                 << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(), "Length")
                 << " " << std::setw( 9)
                 << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(), "Energy")
                 << " " << std::setw( 9)
                 << G4BestUnit((*fSecondary)[lp1]->GetGlobalTime(), "Time")
                 << " " << std::setw(18)
                 << (*fSecondary)[lp1]->GetDefinition()->GetParticleName()
                 << G4endl;
      }
    }
  }
}

//////////////////////////////////////////////////////
void G4SteppingVerbose::PostStepDoItOneByOne()
//////////////////////////////////////////////////////
{
  if(Silent==1) { return; }

  CopyState();

  if(verboseLevel >= 4)
  { 
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

    if(fN2ndariesPostStepDoIt>0)
    {
      for(std::size_t lp1=(*fSecondary).size()-fN2ndariesPostStepDoIt;
                      lp1<(*fSecondary).size(); ++lp1)
      {
        G4cout << "      "
               << std::setw( 9)
               << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(), "Length")
               << " " << std::setw( 9)
               << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(), "Length")
               << " " << std::setw( 9)
               << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(), "Length")
               << " " << std::setw( 9)
               << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(), "Energy")
               << " " << std::setw( 9)
               << G4BestUnit((*fSecondary)[lp1]->GetGlobalTime(), "Time")
               << " " << std::setw(18)
               << (*fSecondary)[lp1]->GetDefinition()->GetParticleName()
               << G4endl;
      }
    }
  }
}

//////////////////////////////////////
void G4SteppingVerbose::VerboseTrack()
//////////////////////////////////////
{
  if(Silent==1){ return; }

  CopyState();

  G4cout << G4endl;
  G4cout << "    ++G4Track Information " << G4endl;
  G4long prec = G4cout.precision(3);


  G4cout << "      -----------------------------------------------" 
         << G4endl;
  G4cout << "        G4Track Information  " << std::setw(20) << G4endl;
  G4cout << "      -----------------------------------------------" 
         << G4endl;

  G4cout << "        Step number         : " 
         << std::setw(20) << fTrack->GetCurrentStepNumber()
         << G4endl; 
#ifdef G4_USE_G4BESTUNIT_FOR_VERBOSE
  G4cout << "        Position - x        : " 
         << std::setw(20) << G4BestUnit(fTrack->GetPosition().x(), "Length")
         << G4endl; 
  G4cout << "        Position - y        : " 
         << std::setw(20) << G4BestUnit(fTrack->GetPosition().y(), "Length")
         << G4endl; 
  G4cout << "        Position - z        : " 
         << std::setw(20) << G4BestUnit(fTrack->GetPosition().z(), "Length")
         << G4endl;
  G4cout << "        Global Time         : " 
         << std::setw(20) << G4BestUnit(fTrack->GetGlobalTime(), "Time")
         << G4endl;
  G4cout << "        Local Time          : " 
         << std::setw(20) << G4BestUnit(fTrack->GetLocalTime(), "Time")
         << G4endl;
#else
  G4cout << "        Position - x (mm)   : " 
         << std::setw(20) << fTrack->GetPosition().x() /mm
         << G4endl; 
  G4cout << "        Position - y (mm)   : " 
         << std::setw(20) << fTrack->GetPosition().y() /mm
         << G4endl; 
  G4cout << "        Position - z (mm)   : " 
         << std::setw(20) << fTrack->GetPosition().z() /mm
         << G4endl;
  G4cout << "        Global Time (ns)    : " 
         << std::setw(20) << fTrack->GetGlobalTime() /ns
         << G4endl;
  G4cout << "        Local Time (ns)     : " 
         << std::setw(20) << fTrack->GetLocalTime() /ns
         << G4endl;
#endif
  G4cout << "        Momentum Direct - x : " 
         << std::setw(20) << fTrack->GetMomentumDirection().x()
         << G4endl;
  G4cout << "        Momentum Direct - y : " 
         << std::setw(20) << fTrack->GetMomentumDirection().y()
         << G4endl;
  G4cout << "        Momentum Direct - z : " 
         << std::setw(20) << fTrack->GetMomentumDirection().z()
         << G4endl;
#ifdef G4_USE_G4BESTUNIT_FOR_VERBOSE
  G4cout << "        Kinetic Energy      : " 
#else
  G4cout << "        Kinetic Energy (MeV): " 
#endif
         << std::setw(20) << G4BestUnit(fTrack->GetKineticEnergy(), "Energy")
         << G4endl;
  G4cout << "        Polarization - x    : " 
         << std::setw(20) << fTrack->GetPolarization().x()
         << G4endl;
  G4cout << "        Polarization - y    : " 
         << std::setw(20) << fTrack->GetPolarization().y()
         << G4endl;
  G4cout << "        Polarization - z    : " 
         << std::setw(20) << fTrack->GetPolarization().z()
         << G4endl;
  G4cout << "        Track Length        : " 
         << std::setw(20) << G4BestUnit(fTrack->GetTrackLength(), "Length")
         << G4endl;
  G4cout << "        Track ID #          : " 
         << std::setw(20) << fTrack->GetTrackID()
         << G4endl;
  G4cout << "        Parent Track ID #   : " 
         << std::setw(20) << fTrack->GetParentID()
         << G4endl;
  G4cout << "        Next Volume         : " 
         << std::setw(20);
  if( fTrack->GetNextVolume() != 0 )
  { 
    G4cout << fTrack->GetNextVolume()->GetName() << " ";
  }
  else
  {
    G4cout << "OutOfWorld" << " ";
  }
  G4cout << G4endl;
  G4cout << "        Track Status        : " 
         << std::setw(20);
  if( fTrack->GetTrackStatus() == fAlive )
  {
    G4cout << " Alive";
  }
  else if( fTrack->GetTrackStatus() == fStopButAlive )
  {
    G4cout << " StopButAlive";
  }
  else if( fTrack->GetTrackStatus() == fStopAndKill )
  {
    G4cout << " StopAndKill";
  }
  else if( fTrack->GetTrackStatus() == fKillTrackAndSecondaries )
  {
    G4cout << " KillTrackAndSecondaries";
  }
  else if( fTrack->GetTrackStatus() == fSuspend )
  {
    G4cout << " Suspend";
  }
  else if( fTrack->GetTrackStatus() == fPostponeToNextEvent )
  {
    G4cout << " PostponeToNextEvent";
  }
  G4cout << G4endl;
#ifdef G4_USE_G4BESTUNIT_FOR_VERBOSE
  G4cout << "        Vertex - x          : " 
         << std::setw(20)
         << G4BestUnit(fTrack->GetVertexPosition().x(),"Length")
         << G4endl; 
  G4cout << "        Vertex - y          : " 
         << std::setw(20)
         << G4BestUnit(fTrack->GetVertexPosition().y(),"Length")
         << G4endl; 
  G4cout << "        Vertex - z          : " 
         << std::setw(20)
         << G4BestUnit(fTrack->GetVertexPosition().z(),"Length")
         << G4endl;
#else
  G4cout << "        Vertex - x (mm)     : " 
         << std::setw(20) << fTrack->GetVertexPosition().x()/mm
         << G4endl; 
  G4cout << "        Vertex - y (mm)     : " 
         << std::setw(20) << fTrack->GetVertexPosition().y()/mm
         << G4endl; 
  G4cout << "        Vertex - z (mm)     : " 
         << std::setw(20) << fTrack->GetVertexPosition().z()/mm
         << G4endl;
#endif
  G4cout << "        Vertex - Px (MomDir): " 
         << std::setw(20) << fTrack->GetVertexMomentumDirection().x()
         << G4endl;
  G4cout << "        Vertex - Py (MomDir): " 
         << std::setw(20) << fTrack->GetVertexMomentumDirection().y()
         << G4endl;
  G4cout << "        Vertex - Pz (MomDir): " 
         << std::setw(20) << fTrack->GetVertexMomentumDirection().z()
         << G4endl;
#ifdef G4_USE_G4BESTUNIT_FOR_VERBOSE
  G4cout << "        Vertex - KineE      : " 
#else
  G4cout << "        Vertex - KineE (MeV): " 
#endif
         << std::setw(20)
         << G4BestUnit(fTrack->GetVertexKineticEnergy(),"Energy")
         << G4endl;
  
  G4cout << "        Creator Process     : " 
         << std::setw(20);
  if( fTrack->GetCreatorProcess() == 0)
  {
    G4cout << " Event Generator" << G4endl;
  }
  else
  {
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
  if(Silent==1) { return; }

  G4cout << G4endl;
  G4cout << "    ++G4ParticleChange Information " << G4endl;
  fParticleChange->DumpInfo();
}

/////////////////////////////////////////
void G4SteppingVerbose::ShowStep() const
////////////////////////////////////////
{
  if(Silent==1){ return; }
   G4String volName;
   G4long oldprc;

   // Show header
   G4cout << G4endl;
   G4cout << "    ++G4Step Information " << G4endl;
   oldprc = G4cout.precision(16);

   // Show G4Step specific information
   G4cout << "      Address of G4Track    : "
          << fStep->GetTrack() << G4endl;
   G4cout << "      Step Length (mm)      : "
          << fStep->GetTrack()->GetStepLength() << G4endl;
   G4cout << "      Energy Deposit (MeV)  : "
          << fStep->GetTotalEnergyDeposit() << G4endl;

   // Show G4StepPoint specific information
   G4cout << "      -------------------------------------------------------" 
          << "----------------" <<  G4endl;
   G4cout << "        StepPoint Information  "
          << std::setw(20) << "PreStep" 
          << std::setw(20) << "PostStep" << G4endl;
   G4cout << "      -------------------------------------------------------" 
          << "----------------" <<  G4endl;
   G4cout << "         Position - x (mm)   : " 
          << std::setw(20) << fStep->GetPreStepPoint()->GetPosition().x() 
          << std::setw(20) << fStep->GetPostStepPoint()->GetPosition().x()
          << G4endl;
   G4cout << "         Position - y (mm)   : " 
          << std::setw(20) << fStep->GetPreStepPoint()->GetPosition().y() 
          << std::setw(20) << fStep->GetPostStepPoint()->GetPosition().y()
          << G4endl;
   G4cout << "         Position - z (mm)   : " 
          << std::setw(20) << fStep->GetPreStepPoint()->GetPosition().z() 
          << std::setw(20) << fStep->GetPostStepPoint()->GetPosition().z()
          << G4endl;
   G4cout << "         Global Time (ns)    : " 
          << std::setw(20) << fStep->GetPreStepPoint()->GetGlobalTime()
          << std::setw(20) << fStep->GetPostStepPoint()->GetGlobalTime()
          << G4endl;
   G4cout << "         Local Time (ns)     : " 
          << std::setw(20) << fStep->GetPreStepPoint()->GetLocalTime() 
          << std::setw(20) << fStep->GetPostStepPoint()->GetLocalTime()
          << G4endl;
   G4cout << "         Proper Time (ns)    : " 
          << std::setw(20) << fStep->GetPreStepPoint()->GetProperTime()
          << std::setw(20) << fStep->GetPostStepPoint()->GetProperTime()
          << G4endl;
   G4cout << "         Momentum Direct - x : " 
       << std::setw(20) << fStep->GetPreStepPoint()->GetMomentumDirection().x()
       << std::setw(20) << fStep->GetPostStepPoint()->GetMomentumDirection().x()
          << G4endl;
   G4cout << "         Momentum Direct - y : " 
       << std::setw(20) << fStep->GetPreStepPoint()->GetMomentumDirection().y()
       << std::setw(20) << fStep->GetPostStepPoint()->GetMomentumDirection().y()
          << G4endl;
   G4cout << "         Momentum Direct - z : " 
       << std::setw(20) << fStep->GetPreStepPoint()->GetMomentumDirection().z()
       << std::setw(20) << fStep->GetPostStepPoint()->GetMomentumDirection().z()
          << G4endl;
   G4cout << "         Momentum - x (MeV/c): " 
          << std::setw(20) << fStep->GetPreStepPoint()->GetMomentum().x()
          << std::setw(20) << fStep->GetPostStepPoint()->GetMomentum().x()
          << G4endl;
   G4cout << "         Momentum - y (MeV/c): " 
          << std::setw(20) << fStep->GetPreStepPoint()->GetMomentum().y()
          << std::setw(20) << fStep->GetPostStepPoint()->GetMomentum().y()
          << G4endl;
   G4cout << "         Momentum - z (MeV/c): " 
          << std::setw(20) << fStep->GetPreStepPoint()->GetMomentum().z()
          << std::setw(20) << fStep->GetPostStepPoint()->GetMomentum().z()
          << G4endl;
   G4cout << "         Total Energy (MeV)  : " 
          << std::setw(20) << fStep->GetPreStepPoint()->GetTotalEnergy()
          << std::setw(20) << fStep->GetPostStepPoint()->GetTotalEnergy()
          << G4endl;
   G4cout << "         Kinetic Energy (MeV): " 
          << std::setw(20) << fStep->GetPreStepPoint()->GetKineticEnergy()
          << std::setw(20) << fStep->GetPostStepPoint()->GetKineticEnergy()
          << G4endl;
   G4cout << "         Velocity (mm/ns)    : " 
          << std::setw(20) << fStep->GetPreStepPoint()->GetVelocity()
          << std::setw(20) << fStep->GetPostStepPoint()->GetVelocity()
          << G4endl;
   G4cout << "         Volume Name         : "
          << std::setw(20)
          << fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
   if (fStep->GetPostStepPoint()->GetPhysicalVolume())
   {
     volName = fStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();
   }
   else
   {
     volName = "OutOfWorld";
   }
   G4cout << std::setw(20) << volName << G4endl;
   G4cout << "         Safety (mm)         : " 
          << std::setw(20) << fStep->GetPreStepPoint()->GetSafety()
          << std::setw(20) << fStep->GetPostStepPoint()->GetSafety()
          << G4endl;
   G4cout << "         Polarization - x    : " 
          << std::setw(20) << fStep->GetPreStepPoint()->GetPolarization().x()
          << std::setw(20) << fStep->GetPostStepPoint()->GetPolarization().x()
          << G4endl;
   G4cout << "         Polarization - y    : " 
          << std::setw(20) << fStep->GetPreStepPoint()->GetPolarization().y()
          << std::setw(20) << fStep->GetPostStepPoint()->GetPolarization().y()
          << G4endl;
   G4cout << "         Polarization - Z    : " 
          << std::setw(20) << fStep->GetPreStepPoint()->GetPolarization().z()
          << std::setw(20) << fStep->GetPostStepPoint()->GetPolarization().z()
          << G4endl;
   G4cout << "         Weight              : " 
          << std::setw(20) << fStep->GetPreStepPoint()->GetWeight()
          << std::setw(20) << fStep->GetPostStepPoint()->GetWeight()
          << G4endl;
   G4cout << "         Step Status         : " ;
   G4StepStatus  tStepStatus = fStep->GetPreStepPoint()->GetStepStatus();
   if( tStepStatus == fGeomBoundary )
   {
     G4cout << std::setw(20) << "Geom Limit";
   }
   else if ( tStepStatus == fAlongStepDoItProc )
   {
     G4cout << std::setw(20) << "AlongStep Proc.";
   }
   else if ( tStepStatus == fPostStepDoItProc )
   {
     G4cout << std::setw(20) << "PostStep Proc";
   }
   else if ( tStepStatus == fAtRestDoItProc )
   {
     G4cout << std::setw(20) << "AtRest Proc";
   }
   else if ( tStepStatus == fUndefined )
   {
     G4cout << std::setw(20) << "Undefined";
   }

   tStepStatus = fStep->GetPostStepPoint()->GetStepStatus();
   if( tStepStatus == fGeomBoundary )
   {
     G4cout << std::setw(20) << "Geom Limit";
   }
   else if ( tStepStatus == fAlongStepDoItProc )
   {
     G4cout << std::setw(20) << "AlongStep Proc.";
   }
   else if ( tStepStatus == fPostStepDoItProc )
   {
     G4cout << std::setw(20) << "PostStep Proc";
   }
   else if ( tStepStatus == fAtRestDoItProc )
   {
     G4cout << std::setw(20) << "AtRest Proc";
   }
   else if ( tStepStatus == fUndefined )
   {
     G4cout << std::setw(20) << "Undefined";
   }

   G4cout << G4endl;
   G4cout << "         Process defined Step: " ;
   if( fStep->GetPreStepPoint()->GetProcessDefinedStep() == 0 )
   {
     G4cout << std::setw(20) << "Undefined";
   }
   else
   {
     G4cout << std::setw(20)
            << fStep->GetPreStepPoint()
                    ->GetProcessDefinedStep()->GetProcessName();
   }
   if( fStep->GetPostStepPoint()->GetProcessDefinedStep() == 0)
   {
     G4cout << std::setw(20) << "Undefined";
   }
   else
   {
     G4cout << std::setw(20)
            << fStep->GetPostStepPoint()
                    ->GetProcessDefinedStep()->GetProcessName(); 
   }
   G4cout.precision(oldprc);

   G4cout << G4endl;
   G4cout << "      -------------------------------------------------------" 
          << "----------------" << G4endl;
}

//////////////////////////////////////////////////
void G4SteppingVerbose::UseBestUnit(G4int prec)
//////////////////////////////////////////////////
{ useBestUnitPrecision = prec; }

//////////////////////////////////////////////////
G4int G4SteppingVerbose::BestUnitPrecision()
//////////////////////////////////////////////////
{ return useBestUnitPrecision; }

