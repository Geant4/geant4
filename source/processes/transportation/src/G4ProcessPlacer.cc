//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ProcessPlacer.cc,v 1.9 2002-09-18 13:52:11 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ProcessPlacer.cc
//
// ----------------------------------------------------------------------

#include "G4ProcessPlacer.hh"
#include "G4ProcessManager.hh"
#include "G4VProcess.hh"
#include "G4ParticleTable.hh"

G4ProcessPlacer::G4ProcessPlacer(const G4String &particlename)
 : fParticleName(particlename)
{}

void G4ProcessPlacer::RemoveProcess(G4VProcess *process){
  G4cout << "=== G4ProcessPlacer::RemoveProcess: for: " <<  fParticleName 
	 << G4endl;
  G4cout << "  ProcessName: " << process->GetProcessName() 
	 << ", will be removed!" << G4endl;

  G4cout << "  The initial Vectors: " << G4endl;
  PrintPostStepGPILVec();
  PrintPostStepDoItVec();

  GetProcessManager().RemoveProcess(process);

  G4cout << "  The final Vectors: " << G4endl;
  PrintPostStepGPILVec();
  PrintPostStepDoItVec();

  G4cout << "================================================" << G4endl;
  
}

void G4ProcessPlacer::AddProcessAs(G4VProcess *process, SecondOrLast sol)
{
  G4cout << "  ProcessName: " << process->GetProcessName() << G4endl;

  G4cout << "The initial Vectors: " << G4endl;
  PrintPostStepGPILVec();
  PrintPostStepDoItVec();


  if (sol == eLast) {  
    GetProcessManager().AddProcess(process,
				   ordInActive,
				   ordInActive,
				   ordLast);
  } 
  else if (sol == eSecond) {
    // get transportation process
    G4VProcess *transportation = 
     (* (GetProcessManager().GetProcessList()))[0];

    if (!transportation) {
      G4Exception(" G4ProcessPlacer:: could not get process id=0");
    }
    if (transportation->GetProcessName() != "Transportation") {
      G4cout << transportation->GetProcessName() << G4endl;
      G4Exception(" G4ProcessPlacer:: process id=0 is not Transportation");
    }

    // place the given proces as first for the moment
    GetProcessManager().AddProcess(process);
    GetProcessManager().SetProcessOrderingToFirst(process, 
						  idxPostStep);
    // place transportation first again
    GetProcessManager().SetProcessOrderingToFirst(transportation, 
						  idxPostStep);
  }
  
  // for verification inly
  G4cout << "The final Vectors: " << G4endl;
  PrintPostStepGPILVec();
  PrintPostStepDoItVec();
  
  G4cout << "================================================" << G4endl;
  
}

void G4ProcessPlacer::AddProcessAsSecondDoIt(G4VProcess *process)
{
  G4cout << "=== G4ProcessPlacer::AddProcessAsSecondDoIt: for: " 
	 << fParticleName << G4endl;
  AddProcessAs(process, eSecond);
}

void G4ProcessPlacer::AddProcessAsLastDoIt(G4VProcess *process)
{
  G4cout << "=== G4ProcessPlacer::AddProcessAsLastDoIt: for: " 
	 << fParticleName << G4endl;
  AddProcessAs(process, eLast);
}

G4ProcessManager &G4ProcessPlacer::GetProcessManager()
{ 
  // get particle iterator to add processes ---------------------
  G4ParticleTable* theParticleTable;
  G4ParticleTable::G4PTblDicIterator* theParticleIterator;
  theParticleTable = G4ParticleTable::GetParticleTable();
  theParticleIterator = theParticleTable->GetIterator();
  // -------------------------------------------------------
  G4ProcessManager *processmanager = 0;
  // add parallel transport processe ---------------------------
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    if (particle->GetParticleName() == fParticleName) {
      processmanager =  particle->GetProcessManager();
      break;
    }
  }
  // ---------------------------------------------------------
  if (!processmanager) G4Exception(" G4ProcessPlacer::GetProcessManager: no ProcessSampler");
  return *processmanager;
}


void G4ProcessPlacer::PrintPostStepGPILVec(){
  G4cout << "GPIL Vector: " << G4endl;
  G4ProcessVector* processGPILVec = 
    GetProcessManager().GetPostStepProcessVector(typeGPIL);
  PrintProcVec(processGPILVec);
} 

void G4ProcessPlacer::PrintPostStepDoItVec(){
  G4cout << "DoIt Vector: " << G4endl;
  G4ProcessVector* processDoItVec = 
    GetProcessManager().GetPostStepProcessVector(typeDoIt); 
  PrintProcVec(processDoItVec);
}


void G4ProcessPlacer::PrintProcVec(G4ProcessVector* processVec)
{
  if (!processVec) G4Exception("G4ProcessPlacer::G4ProcessPlacer: no processVec");
  G4int len = processVec->length();
  if (len==0) G4Exception("G4ProcessPlacer::G4ProcessPlacer:processVec len = 0");
  for (int pi=0; pi<len; pi++) {
    G4VProcess *p = (*processVec)[pi];
    if (p) {
      G4cout << "   " << p->GetProcessName() << G4endl;
    } 
    else {
      G4cout << "   " << "no process found for position: " << pi 
	     << ", in vector of length: " << len << G4endl;
    }
  }
}


