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
// G4ProcessTable class implementation
//
// Author: H.Kurashige, 4 August 1998
// --------------------------------------------------------------------

#include "G4ProcessTableMessenger.hh"
#include "G4ProcessTable.hh"
#include "G4ParticleTable.hh"
#include "G4ios.hh"

// --------------------------------------------------------------------
// Static class variable: ptr to single instance of class in a thread
G4ThreadLocal G4ProcessTable* G4ProcessTable::fProcessTable = nullptr;

// --------------------------------------------------------------------
// Default constructor
//
G4ProcessTable::G4ProcessTable()
{
#ifdef G4VERBOSE
  if (verboseLevel>1)
  {
    G4cout << "--  G4ProcessTable constructor  --" << G4endl;
  }
#endif
  fProcTblVector  = new  G4ProcTableVector();
  fProcNameVector = new  G4ProcNameVector();
  tmpTblVector    = new  G4ProcTableVector();
  fProcTblMessenger = new G4ProcessTableMessenger(this);
}

// --------------------------------------------------------------------
// Destructor
//
G4ProcessTable::~G4ProcessTable()
{
  if ( tmpTblVector != nullptr )
  {
    tmpTblVector ->clear();
    delete tmpTblVector;
    tmpTblVector = nullptr;
  }

  if ( fProcTblVector != nullptr )
  {
    for (auto elem : *fProcTblVector)
    {
      delete elem;
    }  
    fProcTblVector ->clear();
    delete fProcTblVector;
    fProcTblVector = nullptr;
  }

  // delete all process except transportation
  for(auto proc : fListProcesses)
  {
    if ( proc != nullptr )
    {
      G4ProcessType type = proc->GetProcessType();
      if (type != fTransportation && type != fParallel
       && type != fParameterisation)
      {
        delete proc;
      } 
    } 
  }

  fListProcesses.clear();

  if ( fProcNameVector != nullptr )
  {
    fProcNameVector ->clear();
    delete fProcNameVector;
    fProcNameVector = nullptr;
  }
  fProcessTable = nullptr;
  delete fProcTblMessenger;
}

// --------------------------------------------------------------------
//
G4ProcessTable* G4ProcessTable::GetProcessTable()
{
  if(fProcessTable == nullptr)
  {
    static G4ThreadLocalSingleton<G4ProcessTable> inst;
    fProcessTable = inst.Instance();
  }
  return fProcessTable;
}

// --------------------------------------------------------------------
//
G4int G4ProcessTable::Insert(G4VProcess* aProcess, 
                             G4ProcessManager* aProcMgr)
{
  if ( (aProcess == nullptr) || ( aProcMgr == nullptr ) || !fProcTblVector )
  {
#ifdef G4VERBOSE
    if (verboseLevel>0)
    {
      G4cout << "G4ProcessTable::Insert() - arguments are null pointer "
             << aProcess << "," << aProcMgr << G4endl;
    }
#endif
    return -1;
  }
    
#ifdef G4VERBOSE
  if (verboseLevel>1)
  {
    G4cout << "G4ProcessTable::Insert() -";
    G4cout << " Process["  << aProcess->GetProcessName() << "]";
    G4cout << " Particle["  << aProcMgr->GetParticleType()->GetParticleName()
           << "]" << G4endl;
  }
#endif

  G4int idxTbl = 0;
  G4int nidx = (G4int)fProcTblVector->size();
  G4ProcTblElement* anElement = nullptr;
  // loop over all elements
  for (; idxTbl < nidx; ++idxTbl)
  {
    anElement = (*fProcTblVector)[idxTbl];
    if(!anElement) { continue; }
    // check if this process is included
    if (aProcess == anElement->GetProcess())
    {
      // add the process manager into the element 
      // unless  this process manager is included
      if (!anElement->Contains(aProcMgr))
      {
        anElement->Insert(aProcMgr);
#ifdef G4VERBOSE
        if (verboseLevel>2)
        {
          G4cout << " This Process Manager is registered !! " << G4endl;
        }
#endif
      }
      return idxTbl;
    }
  }
  // add this process into the table by creating a new element
  if (verboseLevel>2)
  {
    G4cout << " New element is created !! " << G4endl;
  }
  anElement = new G4ProcTblElement(aProcess);
  anElement->Insert(aProcMgr);
  fProcTblVector->push_back(anElement);
  fProcNameVector->push_back(aProcess->GetProcessName() );
  return nidx;
}

// --------------------------------------------------------------------
//
G4int G4ProcessTable::Remove( G4VProcess* aProcess, 
                              G4ProcessManager* aProcMgr )
{
  if ( (aProcess == nullptr) || ( aProcMgr == nullptr ) || !fProcTblVector )
  {
#ifdef G4VERBOSE
    if (verboseLevel>0)
    {
      G4cout << "G4ProcessTable::Remove() - arguments are null pointer "
             << G4endl;
    }
#endif
    return -1;
  }
    
#ifdef G4VERBOSE
  if (verboseLevel>1)
  {
    G4cout << "G4ProcessTable::Remove() -";
    G4cout << " Process["  << aProcess->GetProcessName() << "]";
    G4cout << " Particle[" << aProcMgr->GetParticleType()->GetParticleName()
           << "]" << G4endl;
  }
#endif

  G4int idxTbl = 0;
  G4int nidx = (G4int)fProcTblVector->size();
  G4ProcTblElement* anElement =nullptr;
  // loop over all elements
  for (; idxTbl < nidx; ++idxTbl)
  {
    anElement = (*fProcTblVector)[idxTbl];
    if(!anElement) { continue; }

    // check if this process is included
    if (aProcess == anElement->GetProcess())
    {
      if(anElement->Contains(aProcMgr))
      {
        // remove the process manager from the element
        anElement->Remove(aProcMgr);
#ifdef G4VERBOSE
        if (verboseLevel>2)
        {
          G4cout << " This Process Manager is removed !! " << G4endl;
        }
#endif
        if(anElement->Length() == 0)
        {
          delete anElement;
          (*fProcTblVector)[idxTbl] = nullptr;
#ifdef G4VERBOSE
          if (verboseLevel>1)
          {
            G4cout << " This Process is removed !! " << G4endl;
          }
#endif
        }
        return idxTbl;
      }
    }
  }
#ifdef G4VERBOSE
  if (verboseLevel>1)
  {
    G4cout << " This Process Manager is not registered to the process!! " 
           << G4endl;
  }
#endif
  return -1;
}

// --------------------------------------------------------------------
//
void G4ProcessTable::RegisterProcess(G4VProcess* ptr)
{
  for(auto proc : fListProcesses)
  {
    if(ptr == proc) { return; } 
  }
  fListProcesses.push_back(ptr);
}

// --------------------------------------------------------------------
//
void G4ProcessTable::DeRegisterProcess(G4VProcess* ptr)
{  
  std::size_t nn = fListProcesses.size();
  for(std::size_t i=0; i<nn; ++i)
  {
    if(ptr == fListProcesses[i])
    { 
      fListProcesses[i] = nullptr;
      return; 
    }
  }
}

// --------------------------------------------------------------------
//
G4VProcess* G4ProcessTable::FindProcess(const G4String& processName, 
                                        const G4String& particleName) const
{
  return FindProcess(processName,
          G4ParticleTable::GetParticleTable()->FindParticle(particleName));
}

// --------------------------------------------------------------------
//
G4VProcess* G4ProcessTable::FindProcess(const G4String& processName, 
                                        const G4ProcessManager* processManager)
                                        const
{
  for (auto anElement : *fProcTblVector)
  {
    // check name and if the processManage is included
    if (anElement && anElement->GetProcessName() == processName
        && anElement->Contains(processManager))
    {
      return anElement->GetProcess();
    }
  }
#ifdef G4VERBOSE
  if (verboseLevel > 1)
  {
    G4cout << " G4ProcessTable::FindProcess() -" ;
    G4cout << " The Process[" << processName << "] is not found  ";
    G4cout << " for [" << processManager->GetParticleType()->GetParticleName()
           << "]" << G4endl;
  }
#endif
  return nullptr;
}

// --------------------------------------------------------------------
//
G4VProcess*
G4ProcessTable::FindProcess(G4ProcessType processType, 
                            const G4ParticleDefinition* particle) const
{
  // find the first process of given type for this particle

  const G4ProcessManager* processManager = particle->GetProcessManager();
  for (auto anElement : *fProcTblVector)
  {
    if (anElement && anElement->GetProcess()->GetProcessType() == processType
        && anElement->Contains(processManager))
    {
      return anElement->GetProcess();
    }
  }
#ifdef G4VERBOSE
  if (verboseLevel > 1)
  {
    G4cout << " G4ProcessTable::FindProcess() -";
    G4cout << " The Process Type " << processType << " is not found  ";
    G4cout << " for [" << particle->GetParticleName() << "]" << G4endl;
  }
#endif
  return nullptr;
}

// --------------------------------------------------------------------
//
G4VProcess*
G4ProcessTable::FindProcess(G4int procSubType, 
                            const G4ParticleDefinition* particle) const
{
  // find the first process of given type for this particle

  const G4ProcessManager* processManager = particle->GetProcessManager();
  for (auto anElement : *fProcTblVector)
  {
    if ( anElement != nullptr
      && anElement->GetProcess()->GetProcessSubType() == procSubType
      && anElement->Contains(processManager) )
    {
      return anElement->GetProcess();
    }
  }
#ifdef G4VERBOSE
  if (verboseLevel > 1)
  {
    G4cout << " G4ProcessTable::FindProcess() -";
    G4cout << " The Process SubType " << procSubType << " is not found  ";
    G4cout << " for [" << particle->GetParticleName() << "]" << G4endl;
  }
#endif
  return nullptr;
}

// --------------------------------------------------------------------
//
G4ProcessTable::G4ProcTableVector*
G4ProcessTable::Find(const G4String& processName)
{
  tmpTblVector->clear();

  G4bool isFound = false;
  G4ProcTblElement* anElement = nullptr;
  for (auto itr=fProcTblVector->cbegin(); itr!=fProcTblVector->cend(); ++itr)
  {
    anElement = (*itr);
    // check name
    if ( anElement != nullptr && anElement->GetProcessName() == processName )
    {
      isFound = true;
      tmpTblVector->push_back(anElement);
    }
  }

  if (!isFound && verboseLevel>0)
  {
#ifdef G4VERBOSE
    G4cout << " G4ProcessTable::Find() -" ;
    G4cout << " The Process[" << processName << "] is not found  " << G4endl;
#endif
  }

  return tmpTblVector;
}

// --------------------------------------------------------------------
//
G4ProcessTable::G4ProcTableVector*
G4ProcessTable::Find(G4ProcessType processType)
{
  tmpTblVector->clear();

  G4bool isFound = false;
  G4ProcTblElement* anElement = nullptr;
  for (auto itr=fProcTblVector->cbegin(); itr!=fProcTblVector->cend(); ++itr)
  {
    anElement = (*itr);
    // check name
    if ( anElement != nullptr && anElement->GetProcess()->GetProcessType() == processType )
    {
      isFound = true;
      tmpTblVector->push_back(anElement);
    }
  }

  if (!isFound && verboseLevel>0)
  {
#ifdef G4VERBOSE
    G4cout << " G4ProcessTable::Find() -" ;
    G4cout << " The ProcessType[" << processType << "] is not found  "
           << G4endl;
#endif
  }

  return tmpTblVector;
}     

// --------------------------------------------------------------------
//
G4ProcessVector*
G4ProcessTable::ExtractProcesses(G4ProcTableVector* procTblVector) const
{
  G4ProcessVector* procList = new G4ProcessVector();
  // loop over all elements
  for (auto itr=procTblVector->cbegin(); itr!=procTblVector->cend(); ++itr)
  {
    G4ProcTblElement* anElement = (*itr);
    if ( anElement != nullptr) procList->insert( anElement->GetProcess() );
  }
  return procList;
}

// --------------------------------------------------------------------
//
void G4ProcessTable::SetProcessActivation(const G4String& processName, 
                                          const G4String& particleName, 
                                                G4bool fActive )
{
  if (particleName == "ALL" )
  {
    SetProcessActivation( processName, fActive); 
  }
  else
  {
    SetProcessActivation(processName, 
        G4ParticleTable::GetParticleTable()->FindParticle(particleName),
        fActive );
  }
}

// --------------------------------------------------------------------
//
void G4ProcessTable::SetProcessActivation(G4ProcessType processType,
                                          const G4String& particleName , 
                                                G4bool fActive)
{
  if ((particleName == "ALL" ) || (particleName == "all" ))
  {
    SetProcessActivation( processType, fActive ); 
  }
  else
  {
    SetProcessActivation(processType,         
        G4ParticleTable::GetParticleTable()->FindParticle(particleName),
        fActive );
  }
}

// --------------------------------------------------------------------
//
void G4ProcessTable::SetProcessActivation( const G4String& processName, 
                                           G4bool          fActive  )
{
#ifdef G4VERBOSE
  if (verboseLevel>1)
  {
    G4cout << " G4ProcessTable::SetProcessActivation() -" ;
    G4cout << " The Process[" << processName << "] "<< G4endl;
  }
#endif

  G4ProcTableVector* pTblVector =  Find(processName);
  G4ProcTblElement* anElement;  
   // loop over all elements
  for (auto itr=pTblVector->cbegin(); itr!=pTblVector->cend(); ++itr)
  {
    anElement = (*itr);
    if ( anElement == nullptr ) continue;
    G4VProcess* process = anElement->GetProcess();
    for (G4int idx = 0 ; idx < anElement->Length(); ++idx)
    {
      G4ProcessManager* manager = anElement->GetProcessManager(idx);
      manager->SetProcessActivation(process, fActive);
#ifdef G4VERBOSE
      if (verboseLevel>1)
      {
        G4cout << "  for " << manager->GetParticleType()->GetParticleName();
        G4cout << "  Index = " << manager->GetProcessIndex(process); 
        G4cout << G4endl;
      }
#endif
    }
  }
}

// --------------------------------------------------------------------
//
void G4ProcessTable::SetProcessActivation(const G4String& processName, 
                                          G4ProcessManager* processManager, 
                                          G4bool fActive  )
{
#ifdef G4VERBOSE
  if (verboseLevel>1)
  {
    G4cout << " G4ProcessTable::SetProcessActivation() -" ;
    G4cout << " The Process[" << processName << "] "<< G4endl;
  }
#endif
  
  G4VProcess* process = FindProcess( processName,  processManager);
  if ( process != nullptr)
  {
    processManager->SetProcessActivation(process, fActive);
#ifdef G4VERBOSE
    if (verboseLevel>1)
    {
      G4cout << "  for "
             << processManager->GetParticleType()->GetParticleName();
      G4cout << "  Index = "
             << processManager->GetProcessIndex(process) << G4endl;
    }
#endif
  } 
}

// --------------------------------------------------------------------
//
void G4ProcessTable::SetProcessActivation( G4ProcessType processType, 
                                           G4bool        fActive  )
{
#ifdef G4VERBOSE
  if (verboseLevel>1)
  {
    G4cout << " G4ProcessTable::SetProcessActivation() -" ;
    G4cout << " The ProcessType[" << G4int(processType) << "] "<< G4endl;
  }
#endif

  G4ProcTableVector* pTblVector = Find(processType);
  G4ProcTblElement* anElement;  
  // loop over all elements
  for (auto itr=pTblVector->cbegin(); itr!=pTblVector->cend(); ++itr)
  {
    anElement = (*itr);
    if ( anElement == nullptr ) continue;
    G4VProcess* process = anElement->GetProcess();
#ifdef G4VERBOSE
    if (verboseLevel>1)
    {
      G4cout << " The Process[" << process->GetProcessName()<< "] "<< G4endl;
    }
#endif
    for (G4int idx = 0 ; idx < anElement->Length(); ++idx)
    {
      G4ProcessManager* manager = anElement->GetProcessManager(idx);
      manager->SetProcessActivation(process, fActive);
#ifdef G4VERBOSE
      if (verboseLevel>1)
      {
        G4cout << "  for " << manager->GetParticleType()->GetParticleName();
        G4cout << "  Index = " << manager->GetProcessIndex(process) << G4endl;
      }
#endif
    }
  }
}

// --------------------------------------------------------------------
//
void G4ProcessTable::SetProcessActivation( G4ProcessType   processType, 
                                           G4ProcessManager* processManager, 
                                           G4bool          fActive  )
{
#ifdef G4VERBOSE
  if (verboseLevel>1)
  {
    G4cout << " G4ProcessTable::SetProcessActivation() -" ;
    G4cout << " The ProcessType[" << G4int(processType) << "] "<< G4endl;
  }
#endif
  
  G4ProcessVector* procList =  processManager->GetProcessList();
  for (G4int idx = 0; idx < (G4int)procList->length(); ++idx)
  {
    G4VProcess* process = (*procList)(idx);
    if ( process->GetProcessType() == processType)
    {
      processManager->SetProcessActivation(process, fActive);
#ifdef G4VERBOSE
      if (verboseLevel>1)
      {
        G4cout << " The Process[" << process->GetProcessName()<< "] "<< G4endl;
        G4cout << "  for "
               << processManager->GetParticleType()->GetParticleName();
        G4cout << "  Index = " << idx << G4endl;
      }
#endif
    }
  }
}

// --------------------------------------------------------------------
//
void G4ProcessTable::DumpInfo(G4VProcess* process, 
                              const G4ParticleDefinition* particle)
{
  G4int idxTbl=0;
  G4ProcTblElement* anElement = nullptr;
  G4bool isFoundInTbl = false;
  G4ProcessManager* manager = nullptr;
  G4int idx;
  // loop over all elements
  for (auto itr=fProcTblVector->cbegin();
            itr!=fProcTblVector->cend(); ++itr, ++idxTbl)
  {
    anElement = (*itr);
    if ( anElement == nullptr ) continue;
    if (process == anElement->GetProcess() )
    {
      if (particle != nullptr)
      {
        for (idx=0; idx<anElement->Length(); ++idx)
        {
          manager = anElement->GetProcessManager(idx);
          if (particle == manager->GetParticleType())
          {
            isFoundInTbl = true;
            break;
          }
        }
      }
      else
      {
        isFoundInTbl = true;
      }
      break;
    }
  }
  if  (!isFoundInTbl ) return;
  
  G4int tmpVerbose = process->GetVerboseLevel();
  process->SetVerboseLevel(verboseLevel);
  process->DumpInfo();
  process->SetVerboseLevel(tmpVerbose);
  if (particle == nullptr)
  {
    for (idx=0; idx<anElement->Length(); ++idx)
    {
      manager = anElement->GetProcessManager(idx);
      G4cout << " for " << manager->GetParticleType()->GetParticleName();
      G4cout << G4endl;
#ifdef G4VERBOSE
      if (verboseLevel >2)
      {
        tmpVerbose = manager->GetVerboseLevel();
        manager->SetVerboseLevel(verboseLevel);
        manager->DumpInfo();
        manager->SetVerboseLevel(tmpVerbose);
      }
#endif
    }
  }
  else
  {
    G4cout << " for " << manager->GetParticleType()->GetParticleName();
    G4cout << G4endl;
#ifdef G4VERBOSE
    if (verboseLevel >2)
    {
      tmpVerbose = manager->GetVerboseLevel();
      manager->SetVerboseLevel(verboseLevel);
      manager->DumpInfo();
      manager->SetVerboseLevel(tmpVerbose);
    }
#endif
  }
}
