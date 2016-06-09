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
// $Id: G4NewProcessPlacer.cc,v 1.1 2006/11/20 10:02:16 ahoward Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4NewProcessPlacer.cc
//
// ----------------------------------------------------------------------

#include "G4NewProcessPlacer.hh"
#include "G4ProcessManager.hh"
#include "G4VProcess.hh"
#include "G4ParticleTable.hh"

G4NewProcessPlacer::G4NewProcessPlacer(const G4String &particlename)
  : fParticleName(particlename)
{
}

G4NewProcessPlacer::~G4NewProcessPlacer()
{
}

void G4NewProcessPlacer::RemoveProcess(G4VProcess *process)
{
  G4cout << "=== G4NewProcessPlacer::RemoveProcess: for: " <<  fParticleName 
         << G4endl;
  G4cout << "  ProcessName: " << process->GetProcessName() 
         << ", will be removed!" << G4endl;

  G4cout << "  The initial Vectors: " << G4endl;
  PrintPostStepGPILVec();
  PrintPostStepDoItVec();

  GetProcessManager()->RemoveProcess(process);

  G4cout << "  The final Vectors: " << G4endl;
  PrintPostStepGPILVec();
  PrintPostStepDoItVec();

  G4cout << "================================================" << G4endl;
  
}

void G4NewProcessPlacer::AddProcessAs(G4VProcess *process, SecondOrLast sol)
{
  G4cout << "  ProcessName: " << process->GetProcessName() << G4endl;
  G4cout << "The initial Vectors: " << G4endl;
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
      G4Exception(" G4NewProcessPlacer:: could not get process id=0");
    }
    if (transportation->GetProcessName() != "Transportation" && transportation->GetProcessName() != "CoupledTransportation")
    {
      //      G4cout << " GOT HERE CoupledTransportation" << G4endl;
      G4cout << transportation->GetProcessName() << G4endl;
      G4Exception(" G4NewProcessPlacer:: process id=0 is not Transportation");
    }

    // place the given proces as first for the moment
    GetProcessManager()->AddProcess(process);
    GetProcessManager()->SetProcessOrderingToFirst(process, 
                                                  idxPostStep);
    // xx test
    //     if(process->GetProcessName() == "ImportanceProcess") 
    GetProcessManager()->SetProcessOrdering(process, 
					    idxAlongStep, 1);
    // place transportation first again
    GetProcessManager()->SetProcessOrderingToFirst(transportation, 
                                                  idxPostStep);
  }
  
  // for verification inly
  G4cout << "The final Vectors: " << G4endl;
  PrintPostStepGPILVec();
  PrintPostStepDoItVec();
  
  G4cout << "================================================" << G4endl;
}

void G4NewProcessPlacer::AddProcessAsSecondDoIt(G4VProcess *process)
{
  G4cout << "=== G4NewProcessPlacer::AddProcessAsSecondDoIt: for: " 
         << fParticleName << G4endl;
  AddProcessAs(process, eSecond);
}

void G4NewProcessPlacer::AddProcessAsLastDoIt(G4VProcess *process)
{
  G4cout << "=== G4NewProcessPlacer::AddProcessAsLastDoIt: for: " 
         << fParticleName << G4endl;
  AddProcessAs(process, eLast);
}

G4ProcessManager *G4NewProcessPlacer::GetProcessManager()
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
  while( (*theParticleIterator)() )
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
    G4Exception("G4NewProcessPlacer::GetProcessManager()", "InvalidSetup",
                FatalException, "NULL pointer to Process Manager !");
  }
  return processmanager;
}

void G4NewProcessPlacer::PrintPostStepGPILVec()
{
  G4cout << "GPIL Vector: " << G4endl;
  G4ProcessVector* processGPILVec = 
    GetProcessManager()->GetPostStepProcessVector(typeGPIL);
  PrintProcVec(processGPILVec);
} 

void G4NewProcessPlacer::PrintPostStepDoItVec()
{
  G4cout << "DoIt Vector: " << G4endl;
  G4ProcessVector* processDoItVec = 
    GetProcessManager()->GetPostStepProcessVector(typeDoIt); 
  PrintProcVec(processDoItVec);
}


void G4NewProcessPlacer::PrintProcVec(G4ProcessVector* processVec)
{
  if (!processVec)
  {
    G4Exception("G4NewProcessPlacer::G4NewProcessPlacer()", "InvalidArgument",
                FatalException, "NULL pointer to process-vector !");
  }
  G4int len = processVec->length();
  if (len==0)
  {
    G4Exception("G4NewProcessPlacer::G4NewProcessPlacer()", "InvalidSetup",
                FatalException, "Length of process-vector is zero !");
  }
  for (int pi=0; pi<len; pi++)
  {
    G4VProcess *p = (*processVec)[pi];
    if (p)
    {
      G4cout << "   " << p->GetProcessName() << G4endl;
    } 
    else
    {
      G4cout << "   " << "no process found for position: " << pi 
             << ", in vector of length: " << len << G4endl;
    }
  }
}
