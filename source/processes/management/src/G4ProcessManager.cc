// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ProcessManager.cc,v 1.2 1999-01-08 11:23:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	For information related to this code contact:  public:
//	CERN, CN Division, ASD Group
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
// ------------------------------------------------------------
//   New Physics scheme           8 Jan. 1997  H.Kurahige
//   remove sprintf              14 Nov 1997  H.Kurahige
//   fixed bugs in FindInsertPosition
//                               18 July 1998 H.Kurashige
// ------------------------------------------------------------

#include "G4ProcessManagerMessenger.hh"
#include "G4ProcessManager.hh"
#include <iomanip.h>
#include "G4ProcessTable.hh"


// ---------------------------------
//  function members implementation
// ---------------------------------
G4ProcessManagerMessenger* G4ProcessManager::fProcessManagerMessenger = NULL;
G4int  G4ProcessManager::counterOfObjects = 0;
// ///////////////////////////////////////
G4ProcessManager::G4ProcessManager(const G4ParticleDefinition* aParticleType):
                theParticleType(aParticleType),
                numberOfProcesses(0),
                duringTracking(false),
		verboseLevel(1)
{
  // create the process List
  theProcessList = new G4ProcessVector();
  if ( theProcessList == NULL) {
    const G4String aErrorMessage("G4ProcessManager::G4ProcessManager():");
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cerr << aErrorMessage;
      G4cerr << " can not create G4ProcessVector " << endl;
    }
#endif
    G4Exception((const char*)aErrorMessage);
  }
  //create process vector
  for (G4int i=0; i<SizeOfProcVectorArray; ++i) {
    theProcVector[i] = new G4ProcessVector();
    if ( theProcVector[i] == NULL) {
      G4String aErrorMessage("G4ProcessManager::G4ProcessManager():");
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0) {
	G4cerr << aErrorMessage;
        G4cerr << " can  not create G4ProcessVector " << endl;
      }
#endif
      G4Exception((const char*)aErrorMessage);
    }
  }
  theAttrVector = new G4ProcessAttrVector();
  if (fProcessManagerMessenger == NULL){
    fProcessManagerMessenger = new G4ProcessManagerMessenger();
#ifdef G4VERBOSE
    if (GetVerboseLevel() > 1) {
      G4cerr << "G4ProcessManagerMessenger is created" << endl;
    } 
#endif
  }
  // Increment counter of G4ProcessManager objects
  counterOfObjects+=1; 
}

// ///////////////////////////////////////
G4ProcessManager::G4ProcessManager(G4ProcessManager &right)
               :duringTracking(false)
{
   verboseLevel = right.verboseLevel;
#ifdef G4VERBOSE
   if (GetVerboseLevel() > 2) {
     G4cerr <<  "G4ProcessManagerMessenger:: copy constructor " <<endl; 
    }
#endif

   theParticleType    = right.theParticleType;
   numberOfProcesses  = right.numberOfProcesses;
 
   // create the process List
   theProcessList = new G4ProcessVector();
   theAttrVector = new G4ProcessAttrVector();
   if ( ( theProcessList == NULL) || (theAttrVector == NULL) ){
     G4String aErrorMessage("G4ProcessManager::G4ProcessManager():");
     G4Exception((const char*)aErrorMessage);
   }

   for (G4int idx=0; idx < right.numberOfProcesses; idx++) {
     // copy contents in theProcessList
     theProcessList->insert((*right.theProcessList)[idx]);
    // create a G4ProcessAttribute same as source's one
     G4ProcessAttribute* sAttr = (*right.theAttrVector)[idx];
     G4ProcessAttribute* dAttr = new G4ProcessAttribute(*sAttr);
     // adds  a G4ProcessAttribute object
     theAttrVector->insert(dAttr);
     numberOfProcesses +=1;
   }
  

  //create theProcVector
  for (G4int i=0; i<SizeOfProcVectorArray; ++i) {
    // create i-th ProcessVector in theProcVector
    theProcVector[i] = new G4ProcessVector();
    if ( theProcVector[i] == NULL) {
      G4String aErrorMessage("G4ProcessManager::G4ProcessManager():");
      G4Exception((const char*)aErrorMessage);
    }
    G4ProcessVector* src = right.theProcVector[i];
    for (G4int j=0; j< src->entries() ; j++){
      // copy j-th process in i-th ProcessVector 
      theProcVector[i]->insert((*src)[j]);
    }
  }
  // Increment counter of G4ProcessManager objects
  counterOfObjects+=1; 
}

// ///////////////////////////////////////
G4ProcessManager::G4ProcessManager():
                theParticleType(NULL),
                numberOfProcesses(0)
{
  if (GetVerboseLevel()>0) {
    G4cerr << "G4ProcessManager: default constructor is called !!" <<endl;
  }
}

// ///////////////////////////////////////
G4ProcessManager & G4ProcessManager::operator=(G4ProcessManager &)
{
  G4Exception("G4ProcessManager: assignment operator is called !!");
  return *this;
}

// ///////////////////////////////////////
G4ProcessManager::~G4ProcessManager()
{
  for (G4int i=0; i<SizeOfProcVectorArray; i++) {
    if (theProcVector[i]) delete theProcVector[i];
  }
  delete theProcessList;

  theAttrVector->clearAndDestroy();
  delete  theAttrVector;

  counterOfObjects-=1; 

  if ( counterOfObjects == 0 ){
    if (fProcessManagerMessenger != NULL){
      delete fProcessManagerMessenger;
      fProcessManagerMessenger = NULL;
#ifdef G4VERBOSE
      if (GetVerboseLevel() > 1) {
	G4cerr << "G4ProcessManagerMessenger is deleted" << endl;
      } 
#endif
    }
  }
}
////////////////////////////////////////////////////////////////
G4int G4ProcessManager::GetProcessVectorIndex(
                           G4VProcess* aProcess,
			   G4ProcessVectorDoItIndex idx,
			   G4ProcessVectorTypeIndex typ
			   ) const
{
  G4int idxVect =  -1;
  G4int idxProc = GetProcessIndex(aProcess); 
  G4int ivec = GetProcessVectorId(idx, typ);
  if ( ( idxProc >=0) && (ivec >=0) ){
    idxVect =  GetAttribute(idxProc)->idxProcVector[ivec];
  } else {
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cerr << " G4ProcessManager::GetProcessVectorIndex:";
      G4cerr << "particle[" << theParticleType->GetParticleName() << "] " ;
      G4cerr <<  "process[" << aProcess->GetProcessName() << "]" ;
      G4cerr << endl;
      if (idxProc <0) { 
	G4cerr << " is not registered yet ";
      }
      if (ivec <0) {
	G4cerr << " illegal DoIt Index [= " << idx << "," << typ << "]";
      }
      G4cerr << endl;
    }
#endif
  }
  return idxVect;
}

/////////////////////////////////////////
G4ProcessAttribute* G4ProcessManager::GetAttribute(G4int index) const
{
  // check index range
  if ((index<0) || (index>=numberOfProcesses)) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cerr << "G4ProcessManager::GetAttribute():";
      G4cerr << " particle[" << theParticleType->GetParticleName() << "]";
      G4cerr << endl;
      G4cerr << "  index out of range " << endl;
      G4cerr << "  #processes[" << numberOfProcesses << "]"; 
      G4cerr << "  index [" << index << "]" << endl;
    }
#endif
    return NULL;
  } 

  // check process pointer is not NULL
  G4VProcess* aProcess = (*theProcessList)[index];
  if (aProcess == NULL) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cerr << "G4ProcessManager::GetAttribute():";
      G4cerr << " particle[" << theParticleType->GetParticleName() << "]";
      G4cerr << endl;
      G4cerr << "process is not defined at " << index << endl;
    }
#endif
    G4String aErrorMessage("G4ProcessManager::GetAttribute():");
    aErrorMessage += " particle[" + theParticleType->GetParticleName() + "]";
    G4Exception((const char*)aErrorMessage); 
    return NULL;
  }    

  //find the process attribute
  if ( ((*theAttrVector)[index])->idxProcessList == index ){
    return  (*theAttrVector)[index];
  } else { 
#ifdef G4VERBOSE
   if (GetVerboseLevel()>0) { 
      G4cerr << "G4ProcessManager::GetAttribute():";
      G4cerr << " particle[" << theParticleType->GetParticleName() << "]";
      G4cerr << endl;
      G4cerr << "Warning:: attribute vector index is inconsistent with process List index" << endl; 
    }
#endif
    G4ProcessAttribute *pAttr = NULL;
    for (G4int i=0; i<theAttrVector->length(); i++){
      if ( ((*theAttrVector)[i])->idxProcessList == index) {
	pAttr = (*theAttrVector)[i];
	break;
      }
    }      
    return pAttr;
  } 
}

// ///////////////////////////////////////
G4ProcessAttribute * G4ProcessManager::GetAttribute(G4VProcess *aProcess) const
{
  return GetAttribute( GetProcessIndex(aProcess));
}

// ///////////////////////////////////////
G4int G4ProcessManager::InsertAt(G4int ip, G4VProcess* process, G4int ivec)
{
  G4ProcessVector* pVector = theProcVector[ivec];
  // check position
  if ( (ip<0) || (ip > pVector->length()) ) return -1;
  // insert in pVector
  pVector->insertAt(ip, process);
  //correct index in ProcessAttributes of processes
  for (G4int iproc=0; iproc<numberOfProcesses; iproc++) {
    G4ProcessAttribute* aAttr = (*theAttrVector)[iproc];
    if (aAttr != NULL) {
      if (aAttr->idxProcVector[ivec] >= ip){
	aAttr->idxProcVector[ivec] += 1;
      }
    } else {
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0) { 
	G4cerr << " G4ProcessManager::InsertAt : No Process Attribute " << endl;
      }
#endif
    }
  }
  return ip;
}

// ///////////////////////////////////////
G4int G4ProcessManager::RemoveAt(G4int ip, G4VProcess* process, G4int ivec)
{
  G4ProcessVector* pVector = theProcVector[ivec];
  // check position
  if ( (ip<0) || (ip >= pVector->length()) ) return -1;
  //check pointer and remove
  if ((*pVector)[ip]== process) {
    pVector->removeAt(ip);
  } else {
    return -1;
  }    
  // correct index
  for(G4int iproc=0; iproc<numberOfProcesses; iproc++) {
    G4ProcessAttribute* aAttr = (*theAttrVector)[iproc];
    if (aAttr != NULL) {
      if (ip < aAttr->idxProcVector[ivec]) {
	aAttr->idxProcVector[ivec] -=1;
      } else if (ip ==  aAttr->idxProcVector[ivec]) {
	aAttr->idxProcVector[ivec] = -1;
	aAttr->ordProcVector[ivec] = ordInActive;
      }
    }else {
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0) { 
	G4cerr << " G4ProcessManager::RemoveAt : No Process Attribute " << endl;
      }
#endif
    }
  } 
  return ip;
}

// ///////////////////////////////////////
G4int G4ProcessManager::FindInsertPosition(G4int ord, G4int ivec)
{
  G4ProcessVector* pVector = theProcVector[ivec];
  G4int ip =  pVector->length();
  G4int tmp = INT_MAX;
  // find insert position
  for (G4int iproc=0; iproc<numberOfProcesses; iproc++) {
    G4ProcessAttribute* aAttr = (*theAttrVector)[iproc];
    if ( (aAttr->ordProcVector[ivec] > ord ) && (tmp > aAttr->ordProcVector[ivec])){
      tmp = aAttr->ordProcVector[ivec] ;
      if (ip > aAttr->idxProcVector[ivec]) ip = aAttr->idxProcVector[ivec];
    }
  }
  return ip;
}

// ///////////////////////////////////////
G4int G4ProcessManager::AddProcess(
		 G4VProcess *aProcess,
		 G4int      ordAtRestDoIt,
		 G4int      ordAlongStepDoIt,
		 G4int      ordPostStepDoIt 
		)
{ 
  G4String aErrorMessage("G4ProcessManager::AddProcess():");
  aErrorMessage += "process[" + aProcess->GetProcessName() + "]";
  aErrorMessage += " particle[" + theParticleType->GetParticleName() + "]";
  
  //check the process is applicable to this particle type
  if ( ( !aProcess->IsApplicable(*theParticleType) ) ||
       (theParticleType->IsShortLived() )                ){
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cerr << aErrorMessage << endl;
      G4cerr << "This process is not applicable to this particle";
    }
#endif
    // --comment out for alpha version   ----
    //  G4Exception((const char*)aErrorMessage); 
    //  return -1;
   }

#ifdef G4VERBOSE
  if (GetVerboseLevel()>2) {
    G4cerr << aErrorMessage << endl;
  }
#endif

  //add aProcess and this ProcessManager into ProcesssTable
  G4ProcessTable* theProcessTable = G4ProcessTable::GetProcessTable();
  theProcessTable->Insert(aProcess, this);

  //add aProcess to process List
  theProcessList->insert(aProcess);  
  G4int idx = (theProcessList->length()) - 1;

  // check size of the ProcessVector[0]
  if (numberOfProcesses != idx){
    theProcessList->removeLast();
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cerr << aErrorMessage << endl;
      G4cerr << "inconsistent process List size" << numberOfProcesses << endl;
    }
#endif
    G4Exception((const char*)aErrorMessage); 
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
  for (G4int ivec=1; ivec<SizeOfProcVectorArray; ivec+=2) {
    if (pAttr->ordProcVector[ivec] < 0 ) {
      // DoIt is inactive if ordering parameter is negative
      pAttr->idxProcVector[ivec] = -1;

    } else {
      //add aProcess in ordering of ordProcVector
      G4ProcessVector* pVector = theProcVector[ivec];
      // find insert position
      G4int ip = FindInsertPosition(pAttr->ordProcVector[ivec], ivec);
      // insert 
      InsertAt(ip, aProcess, ivec);
      // set index in Process Attribute
      pAttr->idxProcVector[ivec] = ip;

#ifdef G4VERBOSE
      if (verboseLevel>2) {
	G4cerr <<  aErrorMessage << endl;
	G4cerr << aProcess->GetProcessName() << " is inserted at "<< ip;
	G4cerr << " in ProcessVetor[" << ivec<< "]";
	G4cerr << " with Ordering parameter = " <<  pAttr->ordProcVector[ivec] ;
	G4cerr << endl;
      }
#endif
     }
  }

  //add ProcessAttribute to ProcessAttrVector
  theAttrVector->insert(pAttr);

  numberOfProcesses += 1;

  CreateGPILvectors();

  return idx;
}


// ///////////////////////////////////////
G4VProcess* G4ProcessManager::RemoveProcess(G4int index)
{
  //find the process attribute
  G4ProcessAttribute* pAttr = GetAttribute(index);
  if (pAttr == NULL) return NULL;

  // remove process
  G4VProcess* removedProcess = (*theProcessList)[index];

  const G4String aErrorMessage(" G4ProcessManager::RemoveProcess():");

  if (pAttr->isActive) {
    // remove process from vectors if the process is active
    for (G4int ivec=0; ivec<SizeOfProcVectorArray; ivec++) {
      G4ProcessVector* pVector = theProcVector[ivec];
      G4int idx = pAttr->idxProcVector[ivec];
      if ((idx >= 0)  && (idx < pVector->length())) {
        //remove
	if (RemoveAt(idx, removedProcess, ivec) <0) {
#ifdef G4VERBOSE
	  if (GetVerboseLevel()>0) {
	    G4cerr << aErrorMessage; 
            G4cerr << "particle["<<theParticleType->GetParticleName()<<"] " ;
	    G4cerr << "process["<<removedProcess->GetProcessName()<< "]  " ;
	    G4cerr << "bad index in attribute ";
	    G4cerr << idx << "->" << pVector->index(removedProcess) << endl;
	  }
#endif
          G4Exception((const char*)aErrorMessage); 
	  return NULL;
	}    
      } else if (idx<0) {
        // corresponding DoIt is not active  
      } else {
        // idx is out of range
#ifdef G4VERBOSE
	if (GetVerboseLevel()>0) {
	  G4cerr << aErrorMessage;
	  G4cerr << "particle["<<theParticleType->GetParticleName()<<"] " ;
	  G4cerr << "process["<<removedProcess->GetProcessName()<< "]  " ;
	  G4cerr << "index(=" << idx << ") of process vector";
	  G4cerr << "[size:" << pVector->length() << "] out of range" <<endl;
	}
#endif
	G4Exception((const char*)aErrorMessage); 
	return NULL;
      }
    }
    pAttr->isActive = false;
  } else { 
    // the process is inactive
  }
  // remove from the process List and delete the attribute
  theProcessList->removeAt(index);
  theAttrVector->remove(pAttr);
  delete pAttr;
  numberOfProcesses -= 1;
  
  // correct index
  for(G4int i=0; i<numberOfProcesses; i++) {
    G4ProcessAttribute* aAttr = (*theAttrVector)[i];
    if (index < aAttr->idxProcessList) aAttr->idxProcessList -=1;
  }

  CreateGPILvectors();

  //remove aProcess from ProcesssTable
  G4ProcessTable* theProcessTable = G4ProcessTable::GetProcessTable();
  theProcessTable->Remove(removedProcess, this);

  return removedProcess;
} 
  
// ///////////////////////////////////////
G4VProcess* G4ProcessManager::RemoveProcess(G4VProcess *aProcess)
{
  return RemoveProcess(GetProcessIndex(aProcess));
} 

/////////////////////////////////////////
G4int G4ProcessManager::GetProcessOrdering(
			G4VProcess *aProcess,
			G4ProcessVectorDoItIndex idDoIt
			)
{
  G4ProcessAttribute* pAttr = GetAttribute(aProcess); 
  G4int ivec = GetProcessVectorId(idDoIt, typeDoIt);
  if (ivec <0 ) {
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cerr << "G4ProcessManager::SetProcessOrdering: ";
      G4cerr << " illegal DoIt Index [= " << idDoIt << "]";
      G4cerr << endl;
      return -1;
    }
#endif
  } else if (pAttr != NULL) { 
    return pAttr->ordProcVector[ivec];
  } else {
    return -1;
  }
  return -1;
}


// ///////////////////////////////////////
void G4ProcessManager::SetProcessOrdering(
			G4VProcess *aProcess,
			G4ProcessVectorDoItIndex idDoIt,
			G4int      ordDoIt
			)
{
  const G4String aErrorMessage(" G4ProcessManager::GetProcessVectorIndex:");

#ifdef G4VERBOSE
  if (GetVerboseLevel()>2) {
    G4cerr << aErrorMessage ;
    G4cerr << "particle[" + theParticleType->GetParticleName() +"] " ;
    G4cerr <<"process[" + aProcess->GetProcessName() + "]"<<  endl;
  }
#endif

  G4ProcessAttribute* pAttr = GetAttribute(aProcess); 
  G4int ivec = GetProcessVectorId(idDoIt, typeDoIt);

  if (ivec <0 ) {
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cerr <<  aErrorMessage << endl;
      G4cerr << "particle[" + theParticleType->GetParticleName() +"] " ;
      G4cerr << "process[" + aProcess->GetProcessName() + "]"<<  endl;
      G4cerr << " illegal DoIt Index [= " << idDoIt << "]";
      G4cerr << endl;
    }
#endif
  } else if (pAttr != NULL) {
    G4int ip = pAttr->idxProcVector[ivec];
    // remove a process from the process vector
    if ( ip >=0 ) {
      RemoveAt(ip, aProcess, ivec);
    }
    // set ordering parameter to non-zero
    if (ordDoIt == 0) ordDoIt = 1;
    pAttr->ordProcVector[ivec-1] = ordDoIt;
    pAttr->ordProcVector[ivec] = ordDoIt;
    // insert in process vector  if ordDoIt >0
    if (ordDoIt >0) {
      // find insert position
      ip = FindInsertPosition(pAttr->ordProcVector[ivec], ivec);
      // insert 
      InsertAt(ip, aProcess, ivec);
      // set index in Process Attribute
      pAttr->idxProcVector[ivec] = ip;
#ifdef G4VERBOSE
      if (verboseLevel>2) {
	G4cerr <<  aErrorMessage << endl;
	G4cerr << "particle[" + theParticleType->GetParticleName() +"] " ;
	G4cerr <<"process[" + aProcess->GetProcessName() + "]"<<  endl;
	G4cerr << aProcess->GetProcessName() << " is inserted at "<< ip;
	G4cerr << " in ProcessVetor[" << ivec<< "]";
	G4cerr << " with Ordering parameter = " <<  ordDoIt ;
	G4cerr << endl;
      }
#endif
    }
  }
  // create GPIL vectors 
  CreateGPILvectors();
}
         

// ///////////////////////////////////////
void G4ProcessManager::SetProcessOrderingToFirst(
			       G4VProcess *aProcess,
			       G4ProcessVectorDoItIndex idDoIt
			       )
{ 
  G4ProcessAttribute* pAttr = GetAttribute(aProcess); 
  G4int ivec = GetProcessVectorId(idDoIt, typeDoIt);
  if (ivec <0 ) {
#ifdef G4VERBOSE
    if (verboseLevel>0) {
      G4cerr << "G4ProcessManager::SetProcessOrdering: ";
      G4cerr << " illegal DoIt Index [= " << idDoIt << "]";
      G4cerr << endl;
    }
#endif
  } else if (pAttr != NULL) {
    G4int ip = pAttr->idxProcVector[ivec];
    // remove a process from the process vector
    if ( ip >=0 ) {
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
    if (verboseLevel>2) {
      G4cerr << "G4ProcessManager::SetProcessOrdering: ";
      G4cerr << aProcess->GetProcessName() << " is inserted at top ";
      G4cerr << " in ProcessVetor[" << ivec<< "]";
      G4cerr << endl;
    }
#endif
  }
  // create GPIL vectors 
  CreateGPILvectors();

}

// ///////////////////////////////////////
void G4ProcessManager::SetProcessOrderingToLast(
			       G4VProcess *aProcess,
			       G4ProcessVectorDoItIndex idDoIt
			       )
{
  SetProcessOrdering(aProcess, idDoIt, ordLast );
}

// ///////////////////////////////////////
G4VProcess* G4ProcessManager::InActivateProcess(G4int index)
{
  //find the process attribute
  G4ProcessAttribute* pAttr = GetAttribute(index);
  if (pAttr == NULL) return NULL;

  // remove process
  G4VProcess* pProcess = (*theProcessList)[index];

  const G4String aErrorMessage(" G4ProcessManager::InactivateProcess():");

  if (pAttr->isActive) {
    // remove process from vectors if the process is active
    for (G4int i=0; i<SizeOfProcVectorArray; i++) {
      G4ProcessVector* pVector = theProcVector[i];
      G4int idx = pAttr->idxProcVector[i];
      if (idx<0) {
        // corresponding DoIt is not active  
      } else if ((idx >= 0)  && (idx < pVector->length())) {
        //check pointer and set to NULL
        if ((*pVector)[idx]== pProcess) {
	  (*pVector)[idx]= NULL;
	} else {
#ifdef G4VERBOSE
	  if (GetVerboseLevel()>0) {
	    G4cerr << aErrorMessage;
	    G4cerr << "particle["<<theParticleType->GetParticleName()<<"] " ;
	    G4cerr << "process["<<pProcess->GetProcessName()<<"]   " ;
	    G4cerr << "bad index in attribute ";
	    G4cerr << idx << "->" << pVector->index(pProcess) << endl;
	  }
#endif
          G4Exception((const char*)aErrorMessage); 
	  return NULL;
	}    
      } else {
        // idx is out of range
#ifdef G4VERBOSE
	if (GetVerboseLevel()>0) {
	  G4cerr << aErrorMessage;
	  G4cerr << "particle["<<theParticleType->GetParticleName()<<"] " ;
	  G4cerr << "process["<<pProcess->GetProcessName()<<"]   " ;
	  G4cerr << "index(=" << idx << ") of process vector";
	  G4cerr << "[size:" << pVector->length() << "] out of range" << endl;
	}
#endif
	G4Exception((const char*)aErrorMessage); 
	return NULL;
      }
    } 
    pAttr->isActive = false;
  } else { 
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cerr << aErrorMessage;
      G4cerr << "particle["<<theParticleType->GetParticleName()<<"] " ;
      G4cerr << "process["<<pProcess->GetProcessName()<<"]   " ;
      G4cerr << "The process is already inactive" << endl;
    }
#endif
  }
  return pProcess;
} 

// ///////////////////////////////////////
G4VProcess* G4ProcessManager::ActivateProcess(G4int index)
{
  //find the process attribute
  G4ProcessAttribute* pAttr = GetAttribute(index);
  if (pAttr == NULL) return NULL;

  // remove process
  G4VProcess* pProcess = (*theProcessList)[index];

  const G4String aErrorMessage(" G4ProcessManager::ActivateProcess():");

  if (pAttr->isActive) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cerr << aErrorMessage << endl;
      G4cerr << "The process is already active" << endl;
    }
#endif
  } else { 
    // remove process from vectors if the process is active
    for (G4int i=0; i<SizeOfProcVectorArray; i++) {
      G4ProcessVector* pVector = theProcVector[i];
      G4int idx = pAttr->idxProcVector[i];
       if (idx<0) {
        // corresponding DoIt is not active  
       } else if ((idx >= 0)  && (idx < pVector->length())) {
        //check pointer and set
	if ((*pVector)[idx]== NULL) {
	  (*pVector)[idx] = pProcess;
	} else {
#ifdef G4VERBOSE
	  if (GetVerboseLevel()>0) {
	    G4cerr << aErrorMessage;
	    G4cerr << "particle["<<theParticleType->GetParticleName()<<"] " ;
	    G4cerr << "process["<<pProcess->GetProcessName()<<"]   " ;
	    G4cerr << "bad index in attribute " << idx << "->";
	    G4cerr << pVector->index(pProcess) << endl;
	  }
#endif
          G4Exception((const char*)aErrorMessage); 
	  return NULL;
	}    
      } else {
        // idx is out of range
#ifdef G4VERBOSE
	if (GetVerboseLevel()>0) {
	  G4cerr << aErrorMessage;
	  G4cerr << "particle["<<theParticleType->GetParticleName()<<"] " ;
	  G4cerr << "process["<<pProcess->GetProcessName()<<"]   " ;
	  G4cerr << "index(=" << idx <<")of process vector";
	  G4cerr << "[size:" << pVector->length() << "] out of range" << endl;
	}
#endif
	G4Exception((const char*)aErrorMessage); 
	return NULL;
      }
    } 
    pAttr->isActive = true;
  }
  return pProcess;
} 

// ///////////////////////////////////////
G4int G4ProcessManager::operator==(const G4ProcessManager &right) const
{
  return (this == &right);
}

// ///////////////////////////////////////
G4int G4ProcessManager::operator!=(const G4ProcessManager &right) const
{
  return (this != &right);
}

// ///////////////////////////////////////
void G4ProcessManager::DumpInfo()
{ 
  // Dump Information

  // particle type
  G4String* aMessage = new G4String("G4ProcessManager:");
  *aMessage += "particle[" + theParticleType->GetParticleName() + "]";
  G4cout << *aMessage << endl;
  delete aMessage;

  // loop over all processes
  for (G4int idx=0; idx <theProcessList->length(); idx++){
    // process name/type
    G4cout << "[" << idx << "]";
    G4cout << "=== process[" << ((*theProcessList)(idx))->GetProcessName()<< " :"; 
    G4cout << G4VProcess::GetProcessTypeName( ((*theProcessList)(idx))->GetProcessType() )<< "]";

    // process attribute    
    G4ProcessAttribute* pAttr = (*theAttrVector)(idx);
    // staus
    if ( pAttr-> isActive ) {
      G4cout << " Active ";
    } else {
      G4cout << " InActive ";
    }
    G4cout << endl;

#ifdef G4VERBOSE
    if (verboseLevel>0) {
      // order parameter    
      G4cout << "  Ordering::     ";
      G4cout << "        AtRest             AlongStep          PostStep   ";
      G4cout << endl;
      G4cout << "                 ";
      G4cout << "   GetPIL/    DoIt    GetPIL/    DoIt    GetPIL/    DoIt ";
      G4cout << endl;
      G4cout << "  Ordering::      " << endl;
      G4cout << "  index           ";
      for (G4int idx2 = 0; idx2 <6 ; idx2++) {
	G4cout << setw(8) << pAttr->idxProcVector[idx2] <<":";
      }
      G4cout << endl;
      G4cout << "  parameter       ";
      for (G4int idx3 = 0; idx3 <6 ; idx3++) {
	G4cout << setw(8) << pAttr->ordProcVector[idx3] <<":";
      }
      G4cout << endl;
    }
#endif 
  }
}

void G4ProcessManager::CreateGPILvectors()
{
//-- create GetPhysicalInteractionLength process vectors just as the inverse
//-- order of DoIt process vector
  for(G4int k=0; k<theProcessList->length(); k++) {
    GetAttribute((*theProcessList)(k))->idxProcVector[0]=-1;
    GetAttribute((*theProcessList)(k))->idxProcVector[2]=-1;
    GetAttribute((*theProcessList)(k))->idxProcVector[4]=-1;
  }
  for(G4int i=0; i<SizeOfProcVectorArray; i += 2) {
    G4ProcessVector* procGPIL = theProcVector[i];
    G4ProcessVector* procDoIt = theProcVector[i+1];
    G4int nproc = procDoIt->entries();
    procGPIL->clear();
    for(G4int j=nproc-1;j>=0;j--) {
      G4VProcess* aProc = (*procDoIt)[j];
      procGPIL->insert(aProc);
      GetAttribute(aProc)->idxProcVector[i] = procGPIL->entries()-1;
    }
  }
}






//////////////////////////////////////////
void G4ProcessManager::StartTracking()
{
  for (G4int idx = 0; idx<theProcessList->length(); idx++){
    if (GetAttribute(idx)->isActive)
      ((*theProcessList)[idx])->StartTracking();
  }
  duringTracking = true;
}

/////////////////////////////////////////////
void G4ProcessManager::EndTracking()
{
  for (G4int idx = 0; idx<theProcessList->length(); idx++){
    if (GetAttribute(idx)->isActive)
      ((*theProcessList)[idx])->EndTracking();
  }
  duringTracking = false;
}


/////////////////////////////////////////////
 G4VProcess* G4ProcessManager::SetProcessActivation(G4VProcess *aProcess, 
						   G4bool      fActive  )
{
  return SetProcessActivation(GetProcessIndex(aProcess), fActive);
} 


/////////////////////////////////////////////
 G4VProcess* G4ProcessManager::SetProcessActivation(G4int index, G4bool fActive)
{
  if (fActive) return ActivateProcess(index);
  else         return InActivateProcess(index);
}
