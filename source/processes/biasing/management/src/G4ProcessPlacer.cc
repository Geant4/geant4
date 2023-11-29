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
{
}

G4ProcessPlacer::~G4ProcessPlacer()
{
}

void G4ProcessPlacer::RemoveProcess(G4VProcess *process)
{
  G4cout << "=== G4ProcessPlacer::RemoveProcess: for: " <<  fParticleName 
         << G4endl;
  G4cout << "  ProcessName: " << process->GetProcessName() 
         << ", will be removed!" << G4endl;

  G4cout << "  The initial AlongStep Vectors: " << G4endl;
  PrintAlongStepGPILVec();
  PrintAlongStepDoItVec();

  G4cout << "  The initial PostStep Vectors: " << G4endl;
  PrintPostStepGPILVec();
  PrintPostStepDoItVec();

  GetProcessManager()->RemoveProcess(process);

  G4cout << "  The final AlongStep Vectors: " << G4endl;
  PrintAlongStepGPILVec();
  PrintAlongStepDoItVec();

  G4cout << "  The final PostStep Vectors: " << G4endl;
  PrintPostStepGPILVec();
  PrintPostStepDoItVec();

  G4cout << "================================================" << G4endl;
  
}

void G4ProcessPlacer::AddProcessAs(G4VProcess *process, SecondOrLast sol)
{
  G4cout << "  Modifying Process Order for ProcessName: " << process->GetProcessName() << G4endl;

  G4cout << "  The initial AlongStep Vectors: " << G4endl;
  PrintAlongStepGPILVec();
  PrintAlongStepDoItVec();

  G4cout << "The initial PostStep Vectors: " << G4endl;
  PrintPostStepGPILVec();
  PrintPostStepDoItVec();

  if (sol == eLast)
  {  
    GetProcessManager()->AddProcess(process, ordInActive, ordInActive, ordLast);
  } 
  else if (sol == eSecond)
  {
    // get transportation process
    G4VProcess *transportation = 
     (* (GetProcessManager()->GetProcessList()))[0];

    if (!transportation)
    {
      G4Exception("G4ProcessPlacer::AddProcessAs","Bias0001",RunMustBeAborted," could not get process id=0");
    }
    if (transportation->GetProcessName() != "Transportation" && transportation->GetProcessName() != "Transportation8" && transportation->GetProcessName() != "CoupledTransportation")
    {
      //      G4cout << " GOT HERE CoupledTransportation" << G4endl;
      G4cout << transportation->GetProcessName() << G4endl;
      G4Exception("G4ProcessPlacer::AddProcessAs","Bias0002",RunMustBeAborted," process id=0 is not Transportation");
    }

    // place the given proces as first for the moment
    // 31/5/11 previously set to first, then transportation set ahead of it, 
    // which is more conveniently correctly set with placing it second!
    GetProcessManager()->AddProcess(process);
    GetProcessManager()->SetProcessOrderingToSecond(process, 
                                                  idxAlongStep);
    GetProcessManager()->SetProcessOrderingToSecond(process, 
                                                  idxPostStep);
    // xx test
    //     if(process->GetProcessName() == "ImportanceProcess") 
    //bug31/10/07    GetProcessManager()->SetProcessOrdering(process, 
    //bug31/10/07					    idxAlongStep, 1);
    // place transportation first again
//     GetProcessManager()->SetProcessOrderingToFirst(transportation, 
//                                                   idxAlongStep);
//     GetProcessManager()->SetProcessOrderingToFirst(transportation, 
//                                                   idxPostStep);
  }
  
  // for verification inly
  G4cout << "  The final AlongStep Vectors: " << G4endl;
  PrintAlongStepGPILVec();
  PrintAlongStepDoItVec();

  G4cout << "The final PostStep Vectors: " << G4endl;
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

G4ProcessManager *G4ProcessPlacer::GetProcessManager()
{ 
  // get particle iterator to add processes ---------------------
  G4ParticleTable* theParticleTable = 0;
  G4ParticleTable::G4PTblDicIterator* theParticleIterator = 0;
  theParticleTable = G4ParticleTable::GetParticleTable();
  theParticleIterator = theParticleTable->GetIterator();
  // -------------------------------------------------------
  G4ProcessManager *processmanager = 0;
  // find process manager ---------------------------
  theParticleIterator->reset();
  while( (*theParticleIterator)() ) /* while checked for unending loop, 30.05.2016, Marc Verderi */
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    if (particle->GetParticleName() == fParticleName)
    {
      processmanager =  particle->GetProcessManager();
      break;
    }
  }
  // ---------------------------------------------------------
  if (!processmanager)
  {
    G4Exception("G4ProcessPlacer::GetProcessManager()", "InvalidSetup",
                FatalException, "NULL pointer to Process Manager ! Sampler.Configure() must be after PhysicsList instantiation");
  }
  return processmanager;
}

void G4ProcessPlacer::PrintAlongStepGPILVec()
{
  G4cout << "GPIL Vector: " << G4endl;
  G4ProcessVector* processGPILVec = 
    GetProcessManager()->GetAlongStepProcessVector(typeGPIL);
  PrintProcVec(processGPILVec);
} 

void G4ProcessPlacer::PrintAlongStepDoItVec()
{
  G4cout << "DoIt Vector: " << G4endl;
  G4ProcessVector* processDoItVec = 
    GetProcessManager()->GetAlongStepProcessVector(typeDoIt); 
  PrintProcVec(processDoItVec);
}


void G4ProcessPlacer::PrintPostStepGPILVec()
{
  G4cout << "GPIL Vector: " << G4endl;
  G4ProcessVector* processGPILVec = 
    GetProcessManager()->GetPostStepProcessVector(typeGPIL);
  PrintProcVec(processGPILVec);
} 

void G4ProcessPlacer::PrintPostStepDoItVec()
{
  G4cout << "DoIt Vector: " << G4endl;
  G4ProcessVector* processDoItVec = 
    GetProcessManager()->GetPostStepProcessVector(typeDoIt); 
  PrintProcVec(processDoItVec);
}


void G4ProcessPlacer::PrintProcVec(G4ProcessVector* processVec)
{
  if (!processVec)
  {
    G4Exception("G4ProcessPlacer::G4ProcessPlacer()", "InvalidArgument",
                FatalException, "NULL pointer to process-vector !");
  }
  G4int len = (G4int)processVec->length();
  if (len==0)
  {
    G4Exception("G4ProcessPlacer::G4ProcessPlacer()", "InvalidSetup",
                FatalException, "Length of process-vector is zero !");
  }
  for (G4int i=0; i<len; ++i)
  {
    G4VProcess *p = (*processVec)[i];
    if (p)
    {
      G4cout << "   " << p->GetProcessName() << G4endl;
    } 
    else
    {
      G4cout << "   " << "no process found for position: " << i 
             << ", in vector of length: " << len << G4endl;
    }
  }
}
