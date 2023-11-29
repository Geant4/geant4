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
// G4ProcessManager implementation
//
// Authors:
// - 02.12.1995, G.Cosmo - First implementation, based on object model
// - 06.05.1996, G.Cosmo - Revised; added vector of processes at rest
// - 08.01.1997, H.Kurashige - New Physics scheme
// --------------------------------------------------------------------

#include <iomanip>

#include "G4ios.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessManagerMessenger.hh"
#include "G4ProcessAttribute.hh"
#include "G4StateManager.hh"
#include "G4ProcessTable.hh"

G4ThreadLocal G4ProcessManagerMessenger*
              G4ProcessManager::fProcessManagerMessenger = nullptr;
G4ThreadLocal G4int G4ProcessManager::counterOfObjects = 0;

// --------------------------------------------------------------------
G4ProcessManager::G4ProcessManager(const G4ParticleDefinition* aParticleType)
  : theParticleType(aParticleType)
{
  // create the process List
  theProcessList = new G4ProcessVector();
  if ( theProcessList == nullptr)
  {
    G4Exception( "G4ProcessManager::G4ProcessManager()","ProcMan012",
                 FatalException, "Can not create G4ProcessList ");
  }

  // create process vector
  for (G4int i=0; i<SizeOfProcVectorArray; ++i)
  {
    theProcVector[i] = new G4ProcessVector();
    if ( theProcVector[i] == nullptr)
    {
      G4Exception( "G4ProcessManager::G4ProcessManager()","ProcMan012",
                   FatalException, "Can not create G4ProcessVector ");
    }
  }

  // create Process Attribute vector
  theAttrVector = new G4ProcessAttrVector();

  // create Process Manager Messenger
  if (fProcessManagerMessenger == nullptr)
  {
    fProcessManagerMessenger = new G4ProcessManagerMessenger();
  }

  for (G4int i=0; i<NDoit; ++i)
  {
    isSetOrderingFirstInvoked[i] = false;
    isSetOrderingLastInvoked[i] = false;
  }

  // Increment counter of G4ProcessManager objects
  ++counterOfObjects; 
}

// --------------------------------------------------------------------
G4ProcessManager::G4ProcessManager(G4ProcessManager& right)
  : theParticleType(right.theParticleType),
    verboseLevel(right.verboseLevel)
{
#ifdef G4VERBOSE
   if (GetVerboseLevel() > 2)
   {
     G4cout << "G4ProcessManager::G4ProcessManager() [copy constructor]"
            << G4endl; 
   }
#endif

   // create the process List and ProcessAttr Vector
   theProcessList = new G4ProcessVector();
   theAttrVector = new G4ProcessAttrVector();
   if ( ( theProcessList == nullptr) || (theAttrVector == nullptr) )
   {
     G4Exception( "G4ProcessManager::G4ProcessManager() [copy constructor]",
                  "ProcMan011",FatalException, "Cannot create G4ProcessList");
   }

   for (G4int idx=0; idx < right.numberOfProcesses; ++idx)
   {
     // copy contents in theProcessList
     theProcessList->insert((*right.theProcessList)[idx]);
     // create a G4ProcessAttribute same as source's one
     G4ProcessAttribute* sAttr = (*right.theAttrVector)[idx];
     G4ProcessAttribute* dAttr = new G4ProcessAttribute(*sAttr);
     // adds  a G4ProcessAttribute object
     theAttrVector->push_back(dAttr);
     ++numberOfProcesses;
   }
  
   // fill up theProcVector
   for (G4int i=0; i<SizeOfProcVectorArray; ++i)
   {
     // create i-th ProcessVector in theProcVector
     theProcVector[i] = new G4ProcessVector();
     if ( theProcVector[i] == nullptr )
     {
       G4Exception("G4ProcessManager::G4ProcessManager() [copy constructor]",
                  "ProcMan011",FatalException, "Cannot create G4ProcessVector");
     }

     G4ProcessTable* theProcessTable = G4ProcessTable::GetProcessTable();
     G4ProcessVector* src = right.theProcVector[i];
     for (G4int j=0; j< (G4int)src->entries() ; ++j)
     {
       // copy j-th process in i-th ProcessVector 
       theProcVector[i]->insert((*src)[j]);
       //add aProcess and this ProcessManager into ProcesssTable 
       if (  (*src)[j] != nullptr )
       {
         theProcessTable->Insert((*src)[j], this);
       }
     }
   }

   for (G4int i=0; i<NDoit; ++i)
   {
     isSetOrderingFirstInvoked[i]= right.isSetOrderingFirstInvoked[i];
     isSetOrderingLastInvoked[i] = right.isSetOrderingLastInvoked[i];
   }

   // Increment counter of G4ProcessManager objects
   ++counterOfObjects; 
}

// --------------------------------------------------------------------
G4ProcessManager::~G4ProcessManager()
{
  for (G4int i=0; i<SizeOfProcVectorArray; ++i)
  {
    if (theProcVector[i])
    {
      theProcVector[i]->clear();
      delete theProcVector[i];
    }
  }
  theProcessList->clear();
  delete theProcessList;

  for (auto itr = theAttrVector->cbegin(); itr!= theAttrVector->cend(); ++itr)
  {
    delete (*itr);
  }
  theAttrVector->clear();
  delete theAttrVector;

  --counterOfObjects; 

  // delete messenger if this object is last one
  if ( counterOfObjects == 0 )
  {
    delete fProcessManagerMessenger;
    fProcessManagerMessenger = nullptr;
#ifdef G4VERBOSE
    if (GetVerboseLevel() > 1)
    {
      G4cout << "G4ProcessManagerMessenger is deleted" << G4endl;
    } 
#endif
  }
}

// --------------------------------------------------------------------
G4int G4ProcessManager::GetProcessVectorIndex(
                           G4VProcess* aProcess,
                           G4ProcessVectorDoItIndex idx,
                           G4ProcessVectorTypeIndex typ
                           ) const
{
  G4int idxVect = -1;
  G4int idxProc = GetProcessIndex(aProcess); 
  G4int ivec = GetProcessVectorId(idx, typ);

  if ( ( idxProc >=0) && (ivec >=0) )
  {
    idxVect = GetAttribute(idxProc)->idxProcVector[ivec];
  }
  else
  {
#ifdef G4VERBOSE
    if (verboseLevel>0)
    {
      G4cout << " G4ProcessManager::GetProcessVectorIndex:";
      G4cout << "particle[" << theParticleType->GetParticleName() << "] " ;
      G4cout <<  "process[" << aProcess->GetProcessName() << "]" ;
      G4cout << G4endl;
      if (idxProc <0)
      { 
        G4cout << " is not registered yet ";
      }
      if (ivec <0)
      {
        G4cout << " illegal DoIt Index [= " << G4int(idx) << ","
                                            << G4int(typ) << "]";
      }
      G4cout << G4endl;
    }
#endif
  }
  return idxVect;
}

// --------------------------------------------------------------------
G4ProcessAttribute* G4ProcessManager::GetAttribute(G4int index) const
{
  // check index range
  if ((index<0) || (index>=numberOfProcesses))
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
    {
      G4cout << "G4ProcessManager::GetAttribute():";
      G4cout << " particle[" << theParticleType->GetParticleName() << "]";
      G4cout << G4endl;
      G4cout << "  index out of range " << G4endl;
      G4cout << "  #processes[" << numberOfProcesses << "]"; 
      G4cout << "  index [" << index << "]" << G4endl;
    }
#endif
    return nullptr;
  } 

  // check process pointer is not null
  G4VProcess* aProcess = (*theProcessList)[index];
  if (aProcess == nullptr)
  {
    G4String aErrorMessage("Bad ProcessList: Null Pointer for ");
    aErrorMessage += theParticleType->GetParticleName() ;
    G4Exception("G4ProcessManager::GetAttribute()","ProcMan012",
                FatalException, aErrorMessage);
    return nullptr;
  }    

  // find the process attribute
  if ( ((*theAttrVector)[index])->idxProcessList == index )
  {
    return  (*theAttrVector)[index];
  }
  else
  { 
    // !! Error !!
    // attribute vector index is inconsistent with process List index
#ifdef G4VERBOSE
   if (GetVerboseLevel()>0)
   { 
      G4cout << "G4ProcessManager::GetAttribute():";
      G4cout << " particle[" << theParticleType->GetParticleName() << "]"
             << G4endl;
      G4cout << "Warning: attribute vector index is inconsistent"
             << " with process List index" 
             << G4endl; 
   }
#endif
   // re-ordering attribute vector 
    G4ProcessAttribute* pAttr = nullptr;
    for (auto itr = theAttrVector->cbegin(); itr!= theAttrVector->cend(); ++itr)
    {
      if ( (*itr)->idxProcessList == index)
      {
        pAttr = (*itr);
        break;
      }
    }      
    return pAttr;
  } 
}

// --------------------------------------------------------------------
G4ProcessAttribute * G4ProcessManager::GetAttribute(G4VProcess* aProcess) const
{
  return GetAttribute( GetProcessIndex(aProcess));
}

// --------------------------------------------------------------------
G4int G4ProcessManager::InsertAt(G4int ip, G4VProcess* process, G4int ivec)
{
  G4ProcessVector* pVector = theProcVector[ivec];
  // check position
  if ( (ip<0) || (ip > G4int(pVector->entries())) ) return -1;
  
  // insert in pVector
  pVector->insertAt(ip, process);

  // correct index in ProcessAttributes of processes
  for (G4int iproc=0; iproc<numberOfProcesses; ++iproc)
  {
    G4ProcessAttribute* aAttr = (*theAttrVector)[iproc];
    if (aAttr != nullptr)
    {
      if (aAttr->idxProcVector[ivec] >= ip)
      {
        aAttr->idxProcVector[ivec] += 1;
      }
    }
    else
    {
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0)
      { 
        G4cout << " G4ProcessManager::InsertAt : No Process Attribute "
               << G4endl;
      }
#endif
    }
  }
  return ip;
}

// --------------------------------------------------------------------
G4int G4ProcessManager::RemoveAt(G4int ip, G4VProcess*, G4int ivec)
{
  G4ProcessVector* pVector = theProcVector[ivec];

  // check position
  if ( (ip<0) || (ip >= G4int(pVector->entries())) ) return -1;

  // remove process
  pVector->removeAt(ip);

  // correct index
  for(G4int iproc=0; iproc<numberOfProcesses; ++iproc)
  {
    G4ProcessAttribute* aAttr = (*theAttrVector)[iproc];
    if (aAttr != nullptr)
    {
      if (ip < aAttr->idxProcVector[ivec])
      {
        aAttr->idxProcVector[ivec] -=1;
      }
      else if (ip ==  aAttr->idxProcVector[ivec])
      {
        aAttr->idxProcVector[ivec] = -1;
        aAttr->ordProcVector[ivec] = ordInActive;
      }
    }
    else
    {
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0)
      { 
        G4cout << " G4ProcessManager::RemoveAt(): No Process Attribute "
               << G4endl;
      }
#endif
    }
  } 
  return ip;
}

// --------------------------------------------------------------------
G4int G4ProcessManager::FindInsertPosition(G4int ord, G4int ivec)
{
  G4ProcessVector* pVector = theProcVector[ivec];
  G4int ip = (G4int)pVector->entries();
  G4int tmp = INT_MAX;
  if (ord == ordLast) return ip;

  // find insert position
  for (G4int iproc=0; iproc<numberOfProcesses; ++iproc)
  {
    G4ProcessAttribute* aAttr = (*theAttrVector)[iproc];
    if ( (aAttr->ordProcVector[ivec] > ord )
      && (tmp > aAttr->ordProcVector[ivec]))
    {
      tmp = aAttr->ordProcVector[ivec] ;
      if (ip > aAttr->idxProcVector[ivec]) ip = aAttr->idxProcVector[ivec];
    }
  }
  return ip;
}

// --------------------------------------------------------------------
G4int G4ProcessManager::AddProcess(
                 G4VProcess* aProcess,
                 G4int      ordAtRestDoIt,
                 G4int      ordAlongStepDoIt,
                 G4int      ordPostStepDoIt 
                )
{
  // check the process is applicable to this particle type
  if (  !aProcess->IsApplicable(*theParticleType) )
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1)
    {
      G4cout << "G4ProcessManager::AddProcess()" << G4endl;
      G4cout << "This process is not applicable to this particle" << G4endl;
    }
#endif
    return -1;
  }

#ifdef G4VERBOSE
  if (GetVerboseLevel()>2)
  {
    G4cout << "G4ProcessManager::AddProcess()" << G4endl;
  }
#endif

  // add aProcess and this ProcessManager into ProcesssTable
  G4ProcessTable* theProcessTable = G4ProcessTable::GetProcessTable();
  theProcessTable->Insert(aProcess, this);

  // add aProcess to process List
  theProcessList->insert(aProcess);  
  G4int idx = G4int(theProcessList->entries() - 1);

  // check size of the ProcessVector[0]
  if (numberOfProcesses != idx)
  {
    theProcessList->removeLast();
    G4String anErrorMessage("Inconsistent process List size for ");
    anErrorMessage += "process[" + aProcess->GetProcessName() + "]";
    anErrorMessage += " particle[" + theParticleType->GetParticleName() + "]";
    G4Exception("G4ProcessManager::AddProcess()", "ProcMan012",
                FatalException, anErrorMessage);
    return -1;
  }

  // create ProcessAttribute
  G4ProcessAttribute* pAttr = new G4ProcessAttribute(aProcess);
  pAttr->idxProcessList = idx;

  // check if ordering parameter is non-zero
  if (ordAtRestDoIt==0)    ordAtRestDoIt    = 1;
  if (ordAlongStepDoIt==0) ordAlongStepDoIt = 1;
  if (ordPostStepDoIt==0)  ordPostStepDoIt  = 1;

  // ordering parameter
  pAttr->ordProcVector[0] = ordAtRestDoIt;
  pAttr->ordProcVector[1] = ordAtRestDoIt;
  pAttr->ordProcVector[2] = ordAlongStepDoIt;
  pAttr->ordProcVector[3] = ordAlongStepDoIt;
  pAttr->ordProcVector[4] = ordPostStepDoIt;
  pAttr->ordProcVector[5] = ordPostStepDoIt;
 
  // add aProccess in Process vectors
  for (G4int ivec=1; ivec<SizeOfProcVectorArray; ivec+=2)
  {
    if (pAttr->ordProcVector[ivec] < 0 )
    {
      // DoIt is inactive if ordering parameter is negative
      pAttr->idxProcVector[ivec] = -1;
    }
    else
    {
      // add aProcess in ordering of ordProcVector
      // G4ProcessVector* pVector = theProcVector[ivec];
      // find insert position
      G4int ip = FindInsertPosition(pAttr->ordProcVector[ivec], ivec);
      // insert 
      InsertAt(ip, aProcess, ivec);
      // set index in Process Attribute
      pAttr->idxProcVector[ivec] = ip;

#ifdef G4VERBOSE
      if (verboseLevel>2)
      {
        G4cout << "G4ProcessManager::AddProcess()" << G4endl;
        G4cout << aProcess->GetProcessName() << " is inserted at "<< ip;
        G4cout << " in ProcessVetor[" << ivec<< "]";
        G4cout << " with Ordering parameter = " ;
        G4cout <<  pAttr->ordProcVector[ivec]  << G4endl;
      }
#endif
    }
  }

  // add ProcessAttribute to ProcessAttrVector
  theAttrVector->push_back(pAttr);

  numberOfProcesses += 1;

  // check consistencies between ordering parameters and process 
  CheckOrderingParameters(aProcess);

  CreateGPILvectors();

  // inform process manager pointer to the process 
  aProcess->SetProcessManager(this);

  return idx;
}

// --------------------------------------------------------------------
G4VProcess* G4ProcessManager::RemoveProcess(G4int index)
{
  //find the process attribute
  G4ProcessAttribute* pAttr = GetAttribute(index);
  if (pAttr == nullptr) return nullptr;

  // remove process
  G4VProcess* removedProcess = (*theProcessList)[index];

  if (!(pAttr->isActive)) { ActivateProcess(index);}
  // remove process from vectors if the process is active
  for (G4int ivec=0; ivec<SizeOfProcVectorArray; ++ivec)
  {
    G4ProcessVector* pVector = theProcVector[ivec];
    G4int idx = pAttr->idxProcVector[ivec];
    if ((idx >= 0) && (idx < G4int(pVector->entries())))
    {
      // remove
      if (RemoveAt(idx, removedProcess, ivec) <0)
      {
        G4String anErrorMessage("Bad index in attribute");
        anErrorMessage += "for particle["
                          + theParticleType->GetParticleName() + "] ";
        anErrorMessage += "process["
                          + removedProcess->GetProcessName() + "]  " ;
        G4Exception( "G4ProcessManager::RemoveProcess()","Fatal Error",
                     FatalException,anErrorMessage); 
        return nullptr;
      }    
    }
    else if (idx<0)
    {
      // corresponding DoIt is not active  
    } 
    else
    {
      // idx is out of range
      G4String anErrorMessage("Bad ProcessList: Index is out of range ");
      anErrorMessage += "for particle[" + theParticleType->GetParticleName() + "] ";
      anErrorMessage += "process[" + removedProcess->GetProcessName() + "]  " ;
      G4Exception( "G4ProcessManager::RemoveProcess()","ProcMan012",
                   FatalException,anErrorMessage); 
      return nullptr;
    }
  }
  pAttr->isActive = false;
  // remove from the process List and delete the attribute
  theProcessList->removeAt(index);
  for (auto itr = theAttrVector->cbegin(); itr!= theAttrVector->cend(); ++itr)
  {
    if ( (*itr) == pAttr)
    {
      theAttrVector->erase(itr);
      break;
    }
  }
  delete pAttr;
  --numberOfProcesses;
  
  // correct index
  for(G4int i=0; i<numberOfProcesses; ++i)
  {
    G4ProcessAttribute* aAttr = (*theAttrVector)[i];
    if (index < aAttr->idxProcessList) aAttr->idxProcessList -=1;
  }

  CreateGPILvectors();

  // remove aProcess from ProcesssTable
  G4ProcessTable* theProcessTable = G4ProcessTable::GetProcessTable();
  theProcessTable->Remove(removedProcess, this);

  return removedProcess;
} 
  
// --------------------------------------------------------------------
G4VProcess* G4ProcessManager::RemoveProcess(G4VProcess *aProcess)
{
  return RemoveProcess(GetProcessIndex(aProcess));
} 

// --------------------------------------------------------------------
G4int G4ProcessManager::GetProcessOrdering(
                        G4VProcess *aProcess,
                        G4ProcessVectorDoItIndex idDoIt
                        )
{
  // get Process Vector Id
  G4int ivec = GetProcessVectorId(idDoIt, typeDoIt);
  if (ivec >=0 )
  {
    // get attribute
    G4ProcessAttribute* pAttr = GetAttribute(aProcess); 
    if (pAttr != nullptr)
    { 
      return pAttr->ordProcVector[ivec];
    }
  }
  return -1;
}

// --------------------------------------------------------------------
void G4ProcessManager::SetProcessOrdering(
                        G4VProcess *aProcess,
                        G4ProcessVectorDoItIndex idDoIt,
                        G4int      ordDoIt
                        )
{
  const G4String aErrorMessage("G4ProcessManager::SetProcessOrdering() - ");

#ifdef G4VERBOSE
  if (GetVerboseLevel()>2)
  {
    G4cout << aErrorMessage ;
    G4cout << "particle[" + theParticleType->GetParticleName() +"] " ;
    G4cout <<"process[" + aProcess->GetProcessName() + "]"<<  G4endl;
  }
#endif

  // get Process Vector Id
  G4int ivec = GetProcessVectorId(idDoIt, typeDoIt);
  if (ivec <0 )
  {
#ifdef G4VERBOSE
    if (verboseLevel>0)
    {
      G4cout <<  aErrorMessage << G4endl;
      G4cout << "particle[" << theParticleType->GetParticleName()  << "] " ;
      G4cout << "process[" << aProcess->GetProcessName() << "]"<<  G4endl;
      G4cout << " illegal DoIt Index [= " << G4int(idDoIt) << "]";
      G4cout << G4endl;
    }
#endif
    return;
  }
 
  if (ordDoIt>ordLast) ordDoIt = ordLast;
  // get attribute 
  G4ProcessAttribute* pAttr = GetAttribute(aProcess); 
  if (pAttr == nullptr)
  {
    // can not get process attribute
    return;
  }
  else
  {
    G4int ip = pAttr->idxProcVector[ivec];
    // remove a process from the process vector
    if ( ip >=0 )
    {
      RemoveAt(ip, aProcess, ivec);
    }

    // set ordering parameter to non-zero
    if (ordDoIt == 0) ordDoIt = 1;
    pAttr->ordProcVector[ivec-1] = ordDoIt;
    pAttr->ordProcVector[ivec] = ordDoIt;

    // insert in process vector  if ordDoIt >0
    if (ordDoIt >0)
    {
      // find insert position
      ip = FindInsertPosition(pAttr->ordProcVector[ivec], ivec);
      // insert 
      InsertAt(ip, aProcess, ivec);
      // set index in Process Attribute
      pAttr->idxProcVector[ivec] = ip;
#ifdef G4VERBOSE
      if (verboseLevel>2)
      {
        G4cout << aErrorMessage << G4endl;
        G4cout << "particle[" << theParticleType->GetParticleName() << "] " ;
        G4cout <<"process[" << aProcess->GetProcessName() << "]"<< G4endl;
        G4cout << aProcess->GetProcessName() << " is inserted at "<< ip;
        G4cout << " in ProcessVetor[" << ivec<< "]";
        G4cout << " with Ordering parameter = " << ordDoIt ;
        G4cout << G4endl;
      }
#endif
    }
  }
  // check consistencies between ordering parameters and process 
  CheckOrderingParameters(aProcess);

  // create GPIL vectors 
  CreateGPILvectors();
}
         
// --------------------------------------------------------------------
void G4ProcessManager::SetProcessOrderingToFirst(
                               G4VProcess *aProcess,
                               G4ProcessVectorDoItIndex idDoIt
                               )
{ 
  // get Process Vector Id
  G4int ivec = GetProcessVectorId(idDoIt, typeDoIt);
  if (ivec <0 )
  {
#ifdef G4VERBOSE
    if (verboseLevel>0)
    {
      G4cout << "G4ProcessManager::SetProcessOrderingToFirst(): ";
      G4cout << " illegal DoIt Index [= " << G4int(idDoIt) << "]";
      G4cout << G4endl;
    }
#endif
    return;
  }

    // get attribute
   G4ProcessAttribute* pAttr = GetAttribute(aProcess); 
   if (pAttr == nullptr)
   {
     return;
   }
   else
   {
     G4int ip = pAttr->idxProcVector[ivec];

     // remove a process from the process vector
     if ( ip >=0 )
     {
       RemoveAt(ip, aProcess, ivec);
     }

     // set ordering parameter to zero
     pAttr->ordProcVector[ivec] = 0;
     pAttr->ordProcVector[ivec-1] = 0;

     // insert 
     InsertAt(0, aProcess, ivec);

     // set index in Process Attribute
     pAttr->idxProcVector[ivec] = 0;

#ifdef G4VERBOSE
     if (verboseLevel>2)
     {
       G4cout << "G4ProcessManager::SetProcessOrderingToFirst(): ";
       G4cout << aProcess->GetProcessName() << " is inserted at top ";
       G4cout << " in ProcessVetor[" << ivec<< "]";
       G4cout << G4endl;
     }
#endif
   }
 
   if (isSetOrderingFirstInvoked[idDoIt])
   {
     G4String anErrMsg = "Set Ordering First is invoked twice for ";
     anErrMsg += aProcess->GetProcessName();
     anErrMsg += " to ";
     anErrMsg += theParticleType->GetParticleName();
     G4Exception("G4ProcessManager::SetProcessOrderingToFirst()",
                 "ProcMan113", JustWarning, anErrMsg); 
   }
   isSetOrderingFirstInvoked[idDoIt] = true;

   // check consistencies between ordering parameters and process 
   CheckOrderingParameters(aProcess);

   // create GPIL vectors 
   CreateGPILvectors();
}

// --------------------------------------------------------------------
void G4ProcessManager::SetProcessOrderingToSecond(
                        G4VProcess *aProcess,
                        G4ProcessVectorDoItIndex idDoIt
                        )
{
  const G4String aErrorMessage("G4ProcessManager::SetProcessOrderingToSecond() - ");
 
#ifdef G4VERBOSE
  if (GetVerboseLevel()>2)
  {
    G4cout << aErrorMessage ;
    G4cout << "particle[" << theParticleType->GetParticleName()  << "] " ;
    G4cout <<"process[" << aProcess->GetProcessName() <<  "]"<<  G4endl;
  }
#endif

  // get Process Vector Id
  G4int ivec = GetProcessVectorId(idDoIt, typeDoIt);
  if (ivec <0 )
  {
#ifdef G4VERBOSE
    if (verboseLevel>0)
    {
      G4cout <<  aErrorMessage << G4endl;
      G4cout << "particle[" << theParticleType->GetParticleName() << "] " ;
      G4cout << "process[" << aProcess->GetProcessName() << "]"<<  G4endl;
      G4cout << " illegal DoIt Index [= " << G4int(idDoIt) << "]";
      G4cout << G4endl;
    }
#endif
    return;
  }
 
  // get attribute 
  G4ProcessAttribute* pAttr = GetAttribute(aProcess); 
  if (pAttr == nullptr)
  {
    // can not get process attribute
    return;
  }
  else
  {
    G4int ip = pAttr->idxProcVector[ivec];
    // remove a process from the process vector
    if ( ip >=0 )
    {
      RemoveAt(ip, aProcess, ivec);
    }
  }

  // set ordering parameter
  pAttr->ordProcVector[ivec-1] = 0;
  pAttr->ordProcVector[ivec]   = 0;

  // find insert position
  G4ProcessVector* pVector = theProcVector[ivec];
  G4int ip = (G4int)pVector->entries();
  G4int tmp = INT_MAX;

  // find insert position
  for (G4int iproc=0; iproc<numberOfProcesses; ++iproc)
  {
    G4ProcessAttribute* aAttr = (*theAttrVector)[iproc];
    if ( aAttr->idxProcVector[ivec] >= 0 )
    {
      if ( (aAttr->ordProcVector[ivec] !=0 )  &&
           (tmp >= aAttr->ordProcVector[ivec]) )
      {
        tmp = aAttr->ordProcVector[ivec];
        if ( ip > aAttr->idxProcVector[ivec] )
        {  
          ip = aAttr->idxProcVector[ivec] ;
        }
      }
    }
  }
  
  // insert 
  InsertAt(ip, aProcess, ivec);

  // set index in Process Attribute
  pAttr->idxProcVector[ivec] = ip;
#ifdef G4VERBOSE
  if (verboseLevel>2)
  {
    G4cout <<  aErrorMessage << G4endl;
    G4cout << "particle[" << theParticleType->GetParticleName()  << "] " ;
    G4cout <<"process[" << aProcess->GetProcessName() << "]"<<  G4endl;
    G4cout << aProcess->GetProcessName() << " is inserted at "<< ip;
    G4cout << " in ProcessVetor[" << ivec<< "]";
    G4cout << " with Ordering parameter = 1 ";
    G4cout << G4endl;
  }
#endif
 
  // check consistencies between ordering parameters and process 
  CheckOrderingParameters(aProcess);

  // create GPIL vectors 
  CreateGPILvectors();
}
         
// --------------------------------------------------------------------
void G4ProcessManager::SetProcessOrderingToLast(
                               G4VProcess *aProcess,
                               G4ProcessVectorDoItIndex idDoIt
                               )
{
  SetProcessOrdering(aProcess, idDoIt, ordLast );

  if (isSetOrderingLastInvoked[idDoIt])
  {
    G4String anErrMsg = "Set Ordering Last is invoked twice for ";
    anErrMsg += aProcess->GetProcessName();
    anErrMsg += " to ";
    anErrMsg += theParticleType->GetParticleName();
    G4Exception( "G4ProcessManager::SetProcessOrderingToLast()","ProcMan114",
                 JustWarning,anErrMsg); 
  }
  isSetOrderingLastInvoked[idDoIt] = true;
}

// --------------------------------------------------------------------
G4VProcess* G4ProcessManager::InActivateProcess(G4int index)
{
  G4ApplicationState currentState 
   = G4StateManager::GetStateManager()->GetCurrentState();
  if ( (currentState == G4State_PreInit) || (currentState == G4State_Init) )
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1)
    {
      G4cout << "G4ProcessManager::InActivateProcess is not valid in ";
      if (currentState == G4State_PreInit )
      {
        G4cout << "PreInit ";
      }
      else if  (currentState == G4State_Init )
      {
        G4cout << "Init ";
      } 
      G4cout << "state !" << G4endl;
    }
#endif
   return nullptr;
  }

  // find the process attribute
  G4ProcessAttribute* pAttr = GetAttribute(index);
  if (pAttr == nullptr) return nullptr;

  // remove process
  G4VProcess* pProcess = (*theProcessList)[index];

  const G4String aErrorMessage("G4ProcessManager::InactivateProcess() - ");

  if (pAttr->isActive)
  {
    // remove process from vectors if the process is active
    for (G4int i=0; i<SizeOfProcVectorArray; ++i)
    {
      G4ProcessVector* pVector = theProcVector[i];
      G4int idx = pAttr->idxProcVector[i];

      if (idx<0)
      {
        // corresponding DoIt is not active  
      }
      else if ((idx >= 0) && (idx < G4int(pVector->entries())))
      {
        //check pointer and set to 0
        if ((*pVector)[idx]== pProcess)
        {
          (*pVector)[idx]= nullptr;
        }
        else
        {
          G4String anErrorMessage("Bad ProcessList: Bad index in attribute");
          anErrorMessage += "for particle[" + theParticleType->GetParticleName() + "] ";
          anErrorMessage += "process[" + pProcess->GetProcessName() + "]  " ;
          G4Exception( "G4ProcessManager::InactivateProcess()", "ProcMan012",
                       FatalException, anErrorMessage); 
          return nullptr;
        }    
      }
      else
      {
        // idx is out of range
        G4String anErrorMessage("Bad ProcessList:  Index is out of range");
        anErrorMessage += "for particle[" + theParticleType->GetParticleName() + "] ";
        anErrorMessage += "process[" + pProcess->GetProcessName() + "]  " ;
        G4Exception( "G4ProcessManager::InactivateProcess()", "ProcMan012",
                     FatalException, anErrorMessage); 
        return nullptr;
      }
    } 
    pAttr->isActive = false;
  }
  return pProcess;
} 

// --------------------------------------------------------------------
G4VProcess* G4ProcessManager::ActivateProcess(G4int index)
{
  G4ApplicationState currentState 
   = G4StateManager::GetStateManager()->GetCurrentState();
  if ( (currentState == G4State_PreInit) || (currentState == G4State_Init) )
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1)
    {
      G4cout << "G4ProcessManager::ActivateProcess() is not valid in ";
      if (currentState == G4State_PreInit )
      {
        G4cout << "PreInit ";
      }
      else  if (currentState == G4State_Init )
      {
        G4cout << "Init ";
      } 
      G4cout << "state !" << G4endl;
    }
#endif
   return nullptr;
  }

  //find the process attribute
  G4ProcessAttribute* pAttr = GetAttribute(index);
  if (pAttr == nullptr) return nullptr;

  // remove process
  G4VProcess* pProcess = (*theProcessList)[index];

  if (!pAttr->isActive)
  {
    // remove process from vectors if the process is active
    for (G4int i=0; i<SizeOfProcVectorArray; ++i)
    {
      G4ProcessVector* pVector = theProcVector[i];
      G4int idx = pAttr->idxProcVector[i];
      if (idx<0)
      {
        // corresponding DoIt is not active  
      }
      else if ((idx >= 0) && (idx < G4int(pVector->entries())))
      {
        // check pointer and set
        if ((*pVector)[idx] == nullptr)
        {
          (*pVector)[idx] = pProcess;
        }
        else
        {
          G4String anErrorMessage("Bad ProcessList: Bad index in attribute");
          anErrorMessage += "for particle[" + theParticleType->GetParticleName() + "] ";
          anErrorMessage += "process[" + pProcess->GetProcessName() + "]  " ;
          G4Exception("G4ProcessManager::ActivateProcess()", "ProcMan012",
                      FatalException, anErrorMessage);
          return nullptr;
        }    
      }
      else
      {
        // idx is out of range
        G4String anErrorMessage("bad ProcessList: Index is out of range");
        anErrorMessage += "for particle["
                          + theParticleType->GetParticleName() + "] ";
        anErrorMessage += "process[" + pProcess->GetProcessName() + "]  " ;
        G4Exception("G4ProcessManager::ActivateProcess()", "ProcMan012",
                    FatalException, anErrorMessage); 
         return nullptr;
      }
    } 
    pAttr->isActive = true;
  }
  return pProcess;
} 

// --------------------------------------------------------------------
G4bool G4ProcessManager::operator==(const G4ProcessManager &right) const
{
  return (this == &right);
}

// --------------------------------------------------------------------
G4bool G4ProcessManager::operator!=(const G4ProcessManager &right) const
{
  return (this != &right);
}

// --------------------------------------------------------------------
void G4ProcessManager::DumpInfo()
{ 
  // Dump Information

  // particle type
  G4cout << "G4ProcessManager: particle[" 
         << theParticleType->GetParticleName() << "]"
         << G4endl;

  // loop over all processes
  for (G4int idx=0; idx < (G4int)theProcessList->entries(); ++idx)
  {
    // process name/type
    G4cout << "[" << idx << "]";
    G4cout << "=== process[" << ((*theProcessList)(idx))->GetProcessName()
           << " :"; 
    G4cout << G4VProcess::GetProcessTypeName( ((*theProcessList)(idx))->GetProcessType() )
           << "]";

    // process attribute    
    G4ProcessAttribute* pAttr = (*theAttrVector)[idx];
    // status
    if ( pAttr-> isActive )
    {
      G4cout << " Active ";
    }
    else
    {
      G4cout << " InActive ";
    }
    G4cout << G4endl;

#ifdef G4VERBOSE
    if (verboseLevel>0)
    {
      // order parameter    
      G4cout << "  Ordering::     ";
      G4cout << "        AtRest             AlongStep          PostStep   ";
      G4cout << G4endl;
      G4cout << "                 ";
      G4cout << "   GetPIL/    DoIt    GetPIL/    DoIt    GetPIL/    DoIt ";
      G4cout << G4endl;
      G4cout << "  Ordering::      " << G4endl;
      G4cout << "  index           ";
      for (G4int idx2 = 0; idx2 <6 ; ++idx2)
      {
        G4cout << std::setw(8) << pAttr->idxProcVector[idx2] << ":";
      }
      G4cout << G4endl;
      G4cout << "  parameter       ";
      for (G4int idx3 = 0; idx3 <6 ; ++idx3)
      {
        G4cout << std::setw(8) << pAttr->ordProcVector[idx3] << ":";
      }
      G4cout << G4endl;
    }
#endif 
  }
}

// --------------------------------------------------------------------
void G4ProcessManager::CreateGPILvectors()
{
  // Create GetPhysicalInteractionLength process vectors just as the inverse
  // order of DoIt() process vector

  for(G4int k=0; k<(G4int)theProcessList->entries(); ++k)
  {
    GetAttribute((*theProcessList)[k])->idxProcVector[0]=-1;
    GetAttribute((*theProcessList)[k])->idxProcVector[2]=-1;
    GetAttribute((*theProcessList)[k])->idxProcVector[4]=-1;
  }

  for(G4int i=0; i<SizeOfProcVectorArray; i += 2)
  {
    G4ProcessVector* procGPIL = theProcVector[i];
    G4ProcessVector* procDoIt = theProcVector[i+1];
    G4int nproc = (G4int)procDoIt->entries();
    procGPIL->clear();
    for(G4int j=nproc-1;j>=0;--j)
    {
      G4VProcess* aProc = (*procDoIt)[j];
      procGPIL->insert(aProc);
      GetAttribute(aProc)->idxProcVector[i] = G4int(procGPIL->entries()-1);
    }
  }
}

// --------------------------------------------------------------------
void G4ProcessManager::StartTracking(G4Track* aTrack)
{
  for (G4int idx = 0; idx<(G4int)theProcessList->entries(); ++idx)
  {
    if (GetAttribute(idx)->isActive)
      ((*theProcessList)[idx])->StartTracking(aTrack);
  }
  if(aTrack) duringTracking = true;
}

// --------------------------------------------------------------------
void G4ProcessManager::EndTracking()
{
  for (G4int idx = 0; idx<(G4int)theProcessList->entries(); ++idx)
  {
    if (GetAttribute(idx)->isActive)
      ((*theProcessList)[idx])->EndTracking();
  }
  duringTracking = false;
}

// --------------------------------------------------------------------
G4VProcess* G4ProcessManager::GetProcess(const G4String& processName) const
{
  for (G4int k=0; k<numberOfProcesses; ++k)
  {
    G4VProcess* process = (*theProcessList)[k];
    if (process->GetProcessName() == processName) return process;
  }
  return nullptr;
}

// --------------------------------------------------------------------
G4VProcess* G4ProcessManager::SetProcessActivation(G4VProcess* aProcess, 
                                                   G4bool      fActive  )
{
  return SetProcessActivation(GetProcessIndex(aProcess), fActive);
} 

// --------------------------------------------------------------------
G4VProcess* G4ProcessManager::SetProcessActivation(G4int index, G4bool fActive)
{
  if (fActive) return ActivateProcess(index);
  else         return InActivateProcess(index);
}

// --------------------------------------------------------------------
 G4bool G4ProcessManager::GetProcessActivation(G4VProcess* aProcess) const
{
  return GetProcessActivation(GetProcessIndex(aProcess));
} 

// --------------------------------------------------------------------
G4bool G4ProcessManager::GetProcessActivation(G4int index) const
{
  if (index <0)
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
    {
      G4cout << "G4ProcessManager::GetProcessActivation  ";
      G4cout << " process (or its index) not found ";
    }
#endif
    return false;
  }
  // process attribute    
  G4ProcessAttribute* pAttr = (*theAttrVector)[index];
  // status
  return pAttr->isActive;
}

// --------------------------------------------------------------------
void G4ProcessManager::CheckOrderingParameters(G4VProcess* aProcess) const
{
  if (aProcess == nullptr) return;
  G4ProcessAttribute* pAttr = GetAttribute(aProcess);
  if (pAttr == nullptr)
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
    {
      G4cout << "G4ProcessManager::CheckOrderingParameters()" << G4endl;
      G4cout << " process " << aProcess->GetProcessName() 
             << " has no attribute" << G4endl;
    }
#endif
    return;
  }

  // check consistencies between ordering parameters and 
  // validity of DoIt of the Process  
  G4bool isOK =true;
  if ( (pAttr->ordProcVector[0]>=0) && (!aProcess->isAtRestDoItIsEnabled()) )
  {
 #ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
    {
      G4cerr << "G4ProcessManager::CheckOrderingParameters()" << G4endl;
      G4cerr << "You cannot set ordering parameter ["
             << pAttr->ordProcVector[0]
             << "] for AtRest DoIt  to the process "
             << aProcess->GetProcessName() << G4endl;
    }
#endif
    isOK = false;
  }

  if ((pAttr->ordProcVector[2]>=0) && (!aProcess->isAlongStepDoItIsEnabled()))
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
    {
      G4cerr << "G4ProcessManager::CheckOrderingParameters()" << G4endl;
      G4cerr << "You cannot set ordering parameter ["
             <<  pAttr->ordProcVector[2]
             << "] for AlongStep DoIt to the process "
             << aProcess->GetProcessName() << G4endl;

    }
#endif
    isOK = false;
  } 

  if ((pAttr->ordProcVector[4]>=0) && (!aProcess->isPostStepDoItIsEnabled()))
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0)
    {
      G4cerr << "G4ProcessManager::CheckOrderingParameters()" << G4endl;
      G4cerr << "You cannot set ordering parameter [" 
             << pAttr->ordProcVector[4] 
             << "] for PostStep DoIt to the process"
             << aProcess->GetProcessName() << G4endl;
    }
#endif
    isOK = false;
  } 

  if (!isOK)
  {
    G4String msg;
    msg = "Invalid ordering parameters are set for  ";
    msg +=  aProcess->GetProcessName();
    G4Exception( "G4ProcessManager::CheckOrderingParameters()",
                 "ProcMan013", FatalException, msg);
  }
  
  return;
}
