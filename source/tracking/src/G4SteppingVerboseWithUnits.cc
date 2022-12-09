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
/// \file 
/// \brief Implementation of the G4SteppingVerboseWithUnits class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//  Stepping Verbose with units for all the applicable double values
//  This class is ported from TestEm2 extended example
//  Original author : Michel Maire (LAPP)
//  Porting with addition of UI command : Makoto Asai (SLAC) Feb.23.2021
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4SteppingVerboseWithUnits.hh"

#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"
#include "G4GenericMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4SteppingVerboseWithUnits::G4SteppingVerboseWithUnits(G4int prec)
  : G4SteppingVerbose(), fprec(prec)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4SteppingVerboseWithUnits::~G4SteppingVerboseWithUnits()
{
  delete fmessenger;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SteppingVerboseWithUnits::SetManager(G4SteppingManager* const fMan)
{
  fManager = fMan;
  fmessenger = new G4GenericMessenger(this,"/tracking/",
               "precision of verbose output");
  auto& cmd = fmessenger->DeclareProperty("setVerbosePrecision",fprec,
               "set precision of verbose output");
  cmd.SetStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SteppingVerboseWithUnits::TrackingStarted()
{  
 CopyState();
 G4long oldprec = G4cout.precision(fprec);
  
 // Step zero
 //  
 if( verboseLevel > 0 )
 {
   G4cout << std::setw(  5)     << "Step#"      << " "
          << std::setw(fprec+3) << "X"          << "    "
          << std::setw(fprec+3) << "Y"          << "    "  
          << std::setw(fprec+3) << "Z"          << "    "
          << std::setw(fprec+6) << "KineE"      << " "
          << std::setw(fprec+10) << "dEStep"     << " "  
          << std::setw(fprec+7) << "StepLeng"  
          << std::setw(fprec+7) << "TrakLeng"
          << std::setw( 10)     << "Volume"     << "  "
          << std::setw( 10)     << "Process"    << G4endl;             

   G4cout << std::setw(5)  << fTrack->GetCurrentStepNumber() << " "
     << std::setw(fprec+3) << G4BestUnit(fTrack->GetPosition().x(),"Length")
     << std::setw(fprec+3) << G4BestUnit(fTrack->GetPosition().y(),"Length")
     << std::setw(fprec+3) << G4BestUnit(fTrack->GetPosition().z(),"Length")
     << std::setw(fprec+3) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
     << std::setw(fprec+7) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
     << std::setw(fprec+3) << G4BestUnit(fStep->GetStepLength(),"Length")
     << std::setw(fprec+3) << G4BestUnit(fTrack->GetTrackLength(),"Length")
     << std::setw(10)      << fTrack->GetVolume()->GetName()
     << std::setw( 9)      << "   initStep" << G4endl;        
  }
  G4cout.precision(oldprec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SteppingVerboseWithUnits::StepInfo()
{  
 CopyState();
 G4long oldprec = G4cout.precision(fprec);

 if( verboseLevel >= 1 )
 {
   if( verboseLevel >= 4 ) VerboseTrack();
   if( verboseLevel >= 3 )
   {
     G4cout << G4endl;    
     G4cout << std::setw(5)       << "#Step#"     << " "
            << std::setw(fprec+3) << "X"          << "    "
            << std::setw(fprec+3) << "Y"          << "    "  
            << std::setw(fprec+3) << "Z"          << "    "
            << std::setw(fprec+6) << "KineE"      << " "
            << std::setw(fprec+10) << "dEStep"     << " "  
            << std::setw(fprec+7) << "StepLeng"     
            << std::setw(fprec+7) << "TrakLeng" 
            << std::setw(10)      << "Volume"    << "  "
            << std::setw(10)      << "Process"   << G4endl;                  
    }

    G4cout << std::setw(5)  << fTrack->GetCurrentStepNumber() << " "
      << std::setw(fprec+3) << G4BestUnit(fTrack->GetPosition().x(),"Length")
      << std::setw(fprec+3) << G4BestUnit(fTrack->GetPosition().y(),"Length")
      << std::setw(fprec+3) << G4BestUnit(fTrack->GetPosition().z(),"Length")
      << std::setw(fprec+3) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
      << std::setw(fprec+7) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
      << std::setw(fprec+3) << G4BestUnit(fStep->GetStepLength(),"Length")
      << std::setw(fprec+3) << G4BestUnit(fTrack->GetTrackLength(),"Length")
      << std::setw(10)      << fTrack->GetVolume()->GetName();

    const G4VProcess* process 
                      = fStep->GetPostStepPoint()->GetProcessDefinedStep();
    G4String procName = " UserLimit";
    if (process != nullptr) procName = process->GetProcessName();
    if (fStepStatus == fWorldBoundary) procName = "OutOfWorld";
    G4cout << "   " << std::setw( 9) << procName;
    G4cout << G4endl;

    if (verboseLevel == 2)
    {
      const std::vector<const G4Track*>* secondary 
                                    = fStep->GetSecondaryInCurrentStep();
      std::size_t nbtrk = (*secondary).size();
      if (nbtrk)
      {
        G4cout << "\n    :----- List of secondaries ----------------" << G4endl;
        G4cout.precision(4);
        for (std::size_t lp=0; lp<(*secondary).size(); ++lp)
        {
          G4cout << "   "
                 << std::setw(13)                 
                 << (*secondary)[lp]->GetDefinition()->GetParticleName()
                 << ":  energy ="
                 << std::setw(6)
                 << G4BestUnit((*secondary)[lp]->GetKineticEnergy(),"Energy")
                 << "  time ="
                 << std::setw(6)
                 << G4BestUnit((*secondary)[lp]->GetGlobalTime(),"Time");
          G4cout << G4endl;
        }
              
        G4cout << "    :------------------------------------------\n" << G4endl;
      }
    }
  }
  G4cout.precision(oldprec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SteppingVerboseWithUnits::AtRestDoItInvoked()
{
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
     
    G4cout << "   Generated secondaries = " << fN2ndariesAtRestDoIt << G4endl;
     
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SteppingVerboseWithUnits::AlongStepDoItAllDone()
{
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SteppingVerboseWithUnits::PostStepDoItAllDone()
{
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
             << "  No. of secondaries = "
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SteppingVerboseWithUnits::DPSLStarted()
{
  CopyState();

  if( verboseLevel > 5 )
  {
    G4cout << G4endl
           << " >>DefinePhysicalStepLength (List of proposed StepLengths): "
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SteppingVerboseWithUnits::DPSLUserLimit()
{
  CopyState();

  if( verboseLevel > 5 )
  {
    G4cout << G4endl << G4endl;
    G4cout << "=== Defined Physical Step Length (DPSL)" << G4endl;
    G4cout << "    ++ProposedStep(UserLimit) = "
           << std::setw( 9) << G4BestUnit(physIntLength , "Length")
           << " : ProcName = User defined maximum allowed Step" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SteppingVerboseWithUnits::DPSLPostStep()
{
  CopyState();

  if( verboseLevel > 5 )
  {
    G4cout << "    ++ProposedStep(PostStep ) = "
           << std::setw( 9) << G4BestUnit(physIntLength , "Length")
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SteppingVerboseWithUnits::DPSLAlongStep()
{
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SteppingVerboseWithUnits::AlongStepDoItOneByOne()
{ 
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
           << "  No. of secondaries = " 
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SteppingVerboseWithUnits::PostStepDoItOneByOne()
{
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
           << "  No. of secondaries = " 
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SteppingVerboseWithUnits::VerboseTrack()
{
  CopyState();

  G4cout << G4endl;
  G4cout << "    ++G4Track Information " << G4endl;
  G4long oldprec = G4cout.precision(fprec);


  G4cout << "      -----------------------------------------------" << G4endl;
  G4cout << "        G4Track Information  " << std::setw(20)        << G4endl;
  G4cout << "      -----------------------------------------------" << G4endl;

  G4cout << "        Step number         : " 
         << std::setw(20) << fTrack->GetCurrentStepNumber() << G4endl; 
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
  G4cout << "        Momentum Direct - x : " 
         << std::setw(20) << fTrack->GetMomentumDirection().x()
         << G4endl;
  G4cout << "        Momentum Direct - y : " 
         << std::setw(20) << fTrack->GetMomentumDirection().y()
         << G4endl;
  G4cout << "        Momentum Direct - z : " 
         << std::setw(20) << fTrack->GetMomentumDirection().z()
         << G4endl;
  G4cout << "        Kinetic Energy      : " 
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
  if (fTrack->GetNextVolume() != 0)  
         G4cout << fTrack->GetNextVolume()->GetName() << " ";
    else G4cout << "OutOfWorld" << " ";
  G4cout << G4endl;
  G4cout << "        Track Status        : " 
         << std::setw(20);
  if( fTrack->GetTrackStatus() == fAlive )
    G4cout << " Alive";
  else if( fTrack->GetTrackStatus() == fStopButAlive )
    G4cout << " StopButAlive";
  else if( fTrack->GetTrackStatus() == fStopAndKill )
    G4cout << " StopAndKill";
  else if( fTrack->GetTrackStatus() == fKillTrackAndSecondaries )
    G4cout << " KillTrackAndSecondaries";
  else if( fTrack->GetTrackStatus() == fSuspend )
    G4cout << " Suspend";
  else if( fTrack->GetTrackStatus() == fPostponeToNextEvent )
    G4cout << " PostponeToNextEvent";
  G4cout << G4endl;
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
  G4cout << "        Vertex - Px (MomDir): " 
         << std::setw(20) << fTrack->GetVertexMomentumDirection().x()
         << G4endl;
  G4cout << "        Vertex - Py (MomDir): " 
         << std::setw(20) << fTrack->GetVertexMomentumDirection().y()
         << G4endl;
  G4cout << "        Vertex - Pz (MomDir): " 
         << std::setw(20) << fTrack->GetVertexMomentumDirection().z()
         << G4endl;
  G4cout << "        Vertex - KineE      : " 
         << std::setw(20)
         << G4BestUnit(fTrack->GetVertexKineticEnergy(),"Energy")
         << G4endl;
  
  G4cout << "        Creator Process     : " 
         << std::setw(20);
  if (fTrack->GetCreatorProcess() == 0)
    G4cout << " Event Generator" << G4endl;
  else
    G4cout << fTrack->GetCreatorProcess()->GetProcessName() << G4endl;

  G4cout << "      -----------------------------------------------"  << G4endl;
  G4cout.precision(oldprec);      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SteppingVerboseWithUnits::VerboseParticleChange()
{
  G4cout << G4endl;
  G4cout << "    ++G4ParticleChange Information " << G4endl;
  fParticleChange->DumpInfo();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SteppingVerboseWithUnits::ShowStep() const
{
   // Show header
   G4cout << G4endl;
   G4cout << "    ++G4Step Information " << G4endl;
   G4long oldprc = G4cout.precision(fprec);

   // Show G4Step specific information
   G4cout << "      Address of G4Track : "
          << fStep->GetTrack() << G4endl;
   G4cout << "      Step Length        : "
          << G4BestUnit(fStep->GetTrack()->GetStepLength(), "Length") << G4endl;
   G4cout << "      Energy Deposit     : "
          << G4BestUnit(fStep->GetTotalEnergyDeposit(), "Energy") << G4endl;

   // Show G4StepPoint specific information
   G4cout << "      -------------------------------------------------------" 
          << "----------------" <<  G4endl;
   G4cout << "        StepPoint Information  "
          << std::setw(20) << "PreStep" << std::setw(20) << "PostStep" << G4endl;
   G4cout << "      -------------------------------------------------------" 
          << "----------------" <<  G4endl;
   G4cout << "         Position - x        : " 
          << std::setw(17) 
	  << G4BestUnit(fStep->GetPreStepPoint()->GetPosition().x(), "Length") 
          << std::setw(17) 
	  << G4BestUnit(fStep->GetPostStepPoint()->GetPosition().x(), "Length")
          << G4endl;
   G4cout << "         Position - y        : " 
          << std::setw(17) 
	  << G4BestUnit(fStep->GetPreStepPoint()->GetPosition().y(), "Length") 
          << std::setw(17) 
	  << G4BestUnit(fStep->GetPostStepPoint()->GetPosition().y(), "Length")
          << G4endl;
   G4cout << "         Position - z        : " 
          << std::setw(17) 
	  << G4BestUnit(fStep->GetPreStepPoint()->GetPosition().z(), "Length") 
          << std::setw(17) 
	  << G4BestUnit(fStep->GetPostStepPoint()->GetPosition().z(), "Length")
          << G4endl;	  
   G4cout << "         Global Time         : " 
          << std::setw(17)
	  << G4BestUnit(fStep->GetPreStepPoint()->GetGlobalTime(), "Time")
          << std::setw(17) 
	  << G4BestUnit(fStep->GetPostStepPoint()->GetGlobalTime(), "Time")
          << G4endl;
   G4cout << "         Local  Time         : " 
          << std::setw(17)
	  << G4BestUnit(fStep->GetPreStepPoint()->GetLocalTime(), "Time")
          << std::setw(17) 
	  << G4BestUnit(fStep->GetPostStepPoint()->GetLocalTime(), "Time")
          << G4endl;
   G4cout << "         Proper Time         : " 
          << std::setw(17)
	  << G4BestUnit(fStep->GetPreStepPoint()->GetProperTime(), "Time")
          << std::setw(17) 
	  << G4BestUnit(fStep->GetPostStepPoint()->GetProperTime(), "Time")
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
   G4cout << "         Momentum - x        : " 
          << std::setw(14) 
	  << G4BestUnit(fStep->GetPreStepPoint()->GetMomentum().x(),"Momentum")
          << std::setw(14) 
	  << G4BestUnit(fStep->GetPostStepPoint()->GetMomentum().x(),"Momentum")
          << G4endl;
   G4cout << "         Momentum - y        : " 
          << std::setw(14) 
	  << G4BestUnit(fStep->GetPreStepPoint()->GetMomentum().y(),"Momentum")
          << std::setw(14) 
	  << G4BestUnit(fStep->GetPostStepPoint()->GetMomentum().y(),"Momentum")
          << G4endl;
   G4cout << "         Momentum - z        : " 
          << std::setw(14) 
	  << G4BestUnit(fStep->GetPreStepPoint()->GetMomentum().z(),"Momentum")
          << std::setw(14) 
	  << G4BestUnit(fStep->GetPostStepPoint()->GetMomentum().z(),"Momentum")
          << G4endl;	  	    
   G4cout << "         Total Energy        : " 
          << std::setw(16) 
	  << G4BestUnit(fStep->GetPreStepPoint()->GetTotalEnergy(),"Energy")
          << std::setw(16) 
	  << G4BestUnit(fStep->GetPostStepPoint()->GetTotalEnergy(),"Energy")
          << G4endl;
   G4cout << "         Kinetic Energy      : " 
          << std::setw(16) 
	  << G4BestUnit(fStep->GetPreStepPoint()->GetKineticEnergy(),"Energy")
          << std::setw(16) 
	  << G4BestUnit(fStep->GetPostStepPoint()->GetKineticEnergy(),"Energy")
          << G4endl;	  
   G4cout << "         Velocity            : " 
          << std::setw(14) 
	  << G4BestUnit(fStep->GetPreStepPoint()->GetVelocity(),"Velocity")
          << std::setw(14) 
	  << G4BestUnit(fStep->GetPostStepPoint()->GetVelocity(),"Velocity")
          << G4endl;
   G4cout << "         Volume Name         : "
          << std::setw(20)
          << fStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
   G4String volName = "OutOfWorld";
   if (fStep->GetPostStepPoint()->GetPhysicalVolume())
     volName = fStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();
   G4cout << std::setw(20) << volName 
          << G4endl;
   G4cout << "         Safety              : " 
          << std::setw(17) 
	  << G4BestUnit(fStep->GetPreStepPoint()->GetSafety(),"Length")
          << std::setw(17) 
	  << G4BestUnit(fStep->GetPostStepPoint()->GetSafety(),"Length")
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
