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
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	History: first implementation, based on object model of
//	4th Aug 1998, H.Kurashige
// ------------------------------------------------------------
//  History:
//   Use STL vector instead of RW vector    1. Mar 00 H.Kurashige
//

#include "G4ProcessTableMessenger.hh"
#include "G4ProcessTable.hh"

//  constructor //////////////////////////
G4ProcessTable::G4ProcessTable()
  :verboseLevel(1)
{
#ifdef G4VERBOSE
  if (verboseLevel>1){
    G4cout << "--  G4ProcessTable constructor  --" << G4endl;
  }
#endif
  fProcTblVector  = new  G4ProcTableVector();
  fProcNameVector = new  G4ProcNameVector();
  tmpTblVector    = new  G4ProcTableVector();
  fProcTblMessenger = nullptr;
}

// copy constructor //////////////////////////
G4ProcessTable::G4ProcessTable(const G4ProcessTable &)
  :verboseLevel(1)
{
  fProcTblVector  = nullptr;
  fProcNameVector = nullptr;
  tmpTblVector    = nullptr;
  fProcTblMessenger = nullptr;
#ifdef G4VERBOSE
  G4cout << "--  G4ProcessTable copy constructor  --" << G4endl;
#endif
}

// destructor //////////////////////////
G4ProcessTable::~G4ProcessTable()
{
#ifdef G4VERBOSE
  if (verboseLevel>1){
    G4cout << "--  G4ProcessTable destructor  --" << G4endl;
  }
#endif

  if ( tmpTblVector != nullptr) {
    tmpTblVector ->clear();
    delete tmpTblVector;
    tmpTblVector = nullptr;
  }

  if ( fProcTblVector != nullptr) {
    G4ProcTableVector::iterator idx;
    
    for (idx=fProcTblVector->begin(); idx!=fProcTblVector->end(); ++idx) {
      delete (*idx);
    }  
    fProcTblVector ->clear();
    delete fProcTblVector;
    fProcTblVector = nullptr;
  }

  if ( fProcNameVector != nullptr) {
    fProcNameVector ->clear();
    delete fProcNameVector;
    fProcNameVector = nullptr;
  }
  fProcessTable = nullptr;
}

/////////////////////////
G4UImessenger* G4ProcessTable::CreateMessenger()
{
  if (fProcTblMessenger == nullptr) {
    fProcTblMessenger = new G4ProcessTableMessenger(this);
  }
  return     fProcTblMessenger;
}

/////////////////////////
void  G4ProcessTable::DeleteMessenger()
{
  if (fProcTblMessenger != nullptr) {
    delete fProcTblMessenger;
    fProcTblMessenger = nullptr;
  }
}


//////////////////////////
G4ProcessTable & G4ProcessTable::operator=(const G4ProcessTable &)
{
#ifdef G4VERBOSE
  if (verboseLevel>0){
    G4cout << "--  G4ProcessTable assignment operator  --" << G4endl;
  }
#endif
  return *this;
}

//////////////////////////
G4bool G4ProcessTable::operator==(const G4ProcessTable &right) const
{
  return (this == &right);
}

//////////////////////////
G4bool G4ProcessTable::operator!=(const G4ProcessTable &right) const
{
  return (this != &right);
}

// Static class variable: ptr to single instance of class
G4ThreadLocal G4ProcessTable* G4ProcessTable::fProcessTable = nullptr;


//////////////////////////
G4ProcessTable* G4ProcessTable::GetProcessTable()
{
    static G4ThreadLocal G4ProcessTable *theProcessTable_G4MT_TLS_ = nullptr ; 
    if (theProcessTable_G4MT_TLS_ == nullptr) {
      theProcessTable_G4MT_TLS_ = new  G4ProcessTable  ;
    }  
    G4ProcessTable &theProcessTable = *theProcessTable_G4MT_TLS_;
    if (fProcessTable == nullptr){
      fProcessTable =  &theProcessTable;
    }
    return fProcessTable;
}

//////////////////////////
G4int   G4ProcessTable::Insert(G4VProcess* aProcess, 
			       G4ProcessManager* aProcMgr)
{
  if ( (aProcess == nullptr) || ( aProcMgr == nullptr ) ){
#ifdef G4VERBOSE
    if (verboseLevel>0){
      G4cout << "G4ProcessTable::Insert : arguments are 0 pointer "
	     <<aProcess <<","<<  aProcMgr << G4endl;
    }
#endif
    return -1;
  }
    
#ifdef G4VERBOSE
  if (verboseLevel>1){
    G4cout << "G4ProcessTable::Insert ";
    G4cout << " Process["  << aProcess->GetProcessName() << "]";
    G4cout << " Particle["  << aProcMgr->GetParticleType()->GetParticleName() << "]";
    G4cout << G4endl;
  }
#endif

  G4ProcTableVector::iterator itr; 
  G4int idxTbl=0;
  G4ProcTblElement* anElement = nullptr;
  G4bool isFoundInTbl = false;
  // loop over all elements
  for (itr=fProcTblVector->begin(); itr!=fProcTblVector->end(); ++itr, ++idxTbl) {
    anElement = (*itr);
    // check if this process is included
    if (aProcess == anElement->GetProcess()) {
      isFoundInTbl = true;

      // add the process manager into the element 
      //  unless  this process manager is included
      if (!anElement->Contains(aProcMgr)) {
	anElement->Insert(aProcMgr);
#ifdef G4VERBOSE
	if (verboseLevel>2){
	  G4cout << " This Process Manager is registered !! " << G4endl;
	}
#endif
      }
      break;
    }
  }
  // add this process into the table by creating a new element
  if (!isFoundInTbl) {
    G4ProcTblElement* newElement = new G4ProcTblElement(aProcess);
    newElement->Insert(aProcMgr);
    fProcTblVector->push_back(newElement);
    // add into name vector     
    G4ProcNameVector::iterator ip;
    G4bool isFound = false;
    for (ip=fProcNameVector->begin(); ip!=fProcNameVector->end(); ++ip) {
      isFound |= (aProcess->GetProcessName() == (*ip));
    }
    if (!isFound) {
      fProcNameVector->push_back(aProcess->GetProcessName() );
#ifdef G4VERBOSE
      if (verboseLevel>2){
	G4cout << " This Process is registered !! " << G4endl;
      }
#endif
    }
  }
  return idxTbl;
}

//////////////////////////
G4int  G4ProcessTable::Remove( G4VProcess* aProcess, 
			       G4ProcessManager* aProcMgr)
{
  if ( (aProcess == nullptr) || ( aProcMgr == nullptr ) ){
#ifdef G4VERBOSE
    if (verboseLevel>0){
      G4cout << "G4ProcessTable::Remove : arguments are 0 pointer "<< G4endl;
    }
#endif
    return -1;
  }
    
#ifdef G4VERBOSE
  if (verboseLevel>1){
    G4cout << "G4ProcessTable::Remove ";
    G4cout << " Process["  << aProcess->GetProcessName() << "]";
    G4cout << " Particle[" << aProcMgr->GetParticleType()->GetParticleName() << "]" << G4endl;
  }
#endif

  G4ProcTableVector::iterator itr; 
  G4int idxTbl=0;
  G4ProcTblElement* anElement =nullptr;
  G4bool isFound = false;
  // loop over all elements
  for (itr=fProcTblVector->begin(); itr!=fProcTblVector->end(); ++itr, ++idxTbl) {
    anElement = (*itr);

    // check if this process is included
    if (aProcess == anElement->GetProcess()) {
      isFound = anElement->Contains(aProcMgr);
      // remove the process manager from the element
      anElement->Remove(aProcMgr);
#ifdef G4VERBOSE
      if (verboseLevel>2){
	G4cout << " This Process Manager is removed !! " << G4endl;
      }
#endif
      break;
    }
  }
  // 
  if (!isFound) {
#ifdef G4VERBOSE
    if (verboseLevel>0){
      G4cout << " This Process Manager is not registered !! " << G4endl;
    }
#endif
    return -1;
  }
  // remove the element if it has no entry
  if (anElement->Length() == 0){
    fProcTblVector->erase(itr);
    delete anElement;
    // check other prcesses with same name exist or not
    G4bool isSameName = false;
    for (itr=fProcTblVector->begin(); itr!=fProcTblVector->end(); ++itr) {
      anElement = (*itr);
      if (anElement->GetProcessName() == aProcess->GetProcessName()) {
	isSameName = true;
	break;
      }
    }
    // remove from name vector
    if (!isSameName ) {
      G4ProcNameVector::iterator i; 
      for (i=fProcNameVector->begin(); i!=fProcNameVector->end(); ++i) {
	if ( *i == aProcess->GetProcessName() ) {
	  fProcNameVector->erase(i);
	  break;
	}
      }
    }
#ifdef G4VERBOSE
    if (verboseLevel>1){
      G4cout << " This Process is removed !! " << G4endl;
    }
#endif
  }
  return idxTbl;
}

//////////////////////////
G4VProcess* G4ProcessTable::FindProcess(const G4String& processName, 
			                const G4ProcessManager* processManager)
                                        const
{
  G4ProcTableVector::iterator itr;   
  G4int idxTbl = 0;
  G4bool isFound = false;
  G4ProcTblElement* anElement =nullptr;
  for (itr=fProcTblVector->begin(); itr!=fProcTblVector->end(); ++itr, ++idxTbl) {
    anElement = (*itr);
    // check name
    if ( anElement->GetProcessName() == processName ) {
      // check if the processManage is included
      if ( anElement->Contains(processManager) ) {
	isFound = true;
	break;
      }
    }
  }
#ifdef G4VERBOSE
  if (!isFound && verboseLevel>1){
    G4cout << " G4ProcessTable::FindProcess :" ;
    G4cout << " The Process[" << processName << "] is not found  ";
    G4cout << " for " << processManager->GetParticleType()->GetParticleName() << G4endl;
  }
#endif
  
  if (isFound) return anElement->GetProcess();
  else         return nullptr;
}

///////////////
G4ProcessTable::G4ProcTableVector* G4ProcessTable::Find( 
					 G4ProcTableVector*,
					 const G4String& processName )
{
  tmpTblVector->clear();

  G4ProcTableVector::iterator itr; 
  G4bool isFound = false;
  G4ProcTblElement* anElement =nullptr;
  for (itr=fProcTblVector->begin(); itr!=fProcTblVector->end(); ++itr) {
    anElement = (*itr);
    // check name
    if ( anElement->GetProcessName() == processName ) {
      isFound = true;
      tmpTblVector->push_back(anElement);
    }
  }

  if (!isFound && verboseLevel>0){
#ifdef G4VERBOSE
    G4cout << " G4ProcessTable::Find :" ;
    G4cout << " The Process[" << processName << "] is not found  " << G4endl;
#endif
  }
  
  return tmpTblVector;

}     
///////////////
G4ProcessTable::G4ProcTableVector* G4ProcessTable::Find( 
					 G4ProcTableVector*,
					 G4ProcessType   processType )
{
  tmpTblVector->clear();

  G4ProcTableVector::iterator itr; 
  G4bool isFound = false;
  G4ProcTblElement* anElement =nullptr;
  for (itr=fProcTblVector->begin(); itr!=fProcTblVector->end(); ++itr) {
    anElement = (*itr);
    // check name
    if ( anElement->GetProcess()->GetProcessType() == processType ) {
      isFound = true;
      tmpTblVector->push_back(anElement);
    }
  }

  if (!isFound && verboseLevel>0){
#ifdef G4VERBOSE
    G4cout << " G4ProcessTable::Find :" ;
    G4cout << " The ProcessType[" << processType << "] is not found  " << G4endl;
#endif
  }
  
  return tmpTblVector;

}     

///////////////
G4ProcessVector* G4ProcessTable::ExtractProcesses( G4ProcTableVector* procTblVector)
{
  G4ProcessVector* procList = new G4ProcessVector();
  G4ProcTableVector::iterator itr; 
  // loop over all elements
  for (itr=procTblVector->begin(); itr!=procTblVector->end(); ++itr) {
    G4ProcTblElement* anElement = (*itr);
    procList->insert( anElement->GetProcess() );
  }
  return procList;
}

///////////////
G4ProcessVector* G4ProcessTable::FindProcesses()
{
  return ExtractProcesses(fProcTblVector);
}

///////////////
G4ProcessVector* G4ProcessTable::FindProcesses( const G4ProcessManager* pManager )
{
  G4ProcessVector* procList = pManager->GetProcessList();
  return new G4ProcessVector(*procList);
}

///////////////
G4ProcessVector* G4ProcessTable::FindProcesses( const G4String& processName )
{
  G4ProcTableVector* pTblVector =  Find(fProcTblVector, processName);
  return ExtractProcesses(pTblVector);
}

///////////////
G4ProcessVector* G4ProcessTable::FindProcesses( G4ProcessType processType)
{
  G4ProcTableVector* pTblVector =  Find(fProcTblVector, processType);
  return ExtractProcesses(pTblVector);
}

///////////////
void G4ProcessTable::SetProcessActivation( const G4String& processName, 
			                   G4bool          fActive  )
{
#ifdef G4VERBOSE
  if (verboseLevel>1){
    G4cout << " G4ProcessTable::SetProcessActivation:" ;
    G4cout << " The Process[" << processName << "] "<< G4endl;
  }
#endif

  G4ProcTableVector* pTblVector =  Find(fProcTblVector, processName);
  G4ProcTableVector::iterator itr; 
  G4ProcTblElement* anElement;  
   // loop over all elements
  for (itr=pTblVector->begin(); itr!=pTblVector->end(); ++itr) {
    anElement = (*itr);
    G4VProcess* process = anElement->GetProcess();
    for (G4int idx = 0 ; idx < anElement->Length(); idx++) {
      G4ProcessManager* manager = anElement->GetProcessManager(idx);
      manager->SetProcessActivation(process, fActive);
#ifdef G4VERBOSE
      if (verboseLevel>1){
        G4cout << "  for " << manager->GetParticleType()->GetParticleName();
	G4cout << "  Index = " << manager->GetProcessIndex(process); 
        G4cout << G4endl;
      }
#endif
    }
  }
}

///////////////
void G4ProcessTable::SetProcessActivation( 
			   const G4String& processName, 
		           G4ProcessManager* processManager, 
			   G4bool          fActive  )
{
#ifdef G4VERBOSE
  if (verboseLevel>1){
    G4cout << " G4ProcessTable::SetProcessActivation:" ;
    G4cout << " The Process[" << processName << "] "<< G4endl;
  }
#endif
  
  G4VProcess* process = FindProcess( processName,  processManager);
  if ( process != nullptr) {
    processManager->SetProcessActivation(process, fActive);
#ifdef G4VERBOSE
    if (verboseLevel>1){
      G4cout << "  for " << processManager->GetParticleType()->GetParticleName();
      G4cout << "  Index = " << processManager->GetProcessIndex(process) << G4endl;
    }
#endif
  } 
}


///////////////
void G4ProcessTable::SetProcessActivation( G4ProcessType   processType, 
			                   G4bool          fActive  )
{
#ifdef G4VERBOSE
  if (verboseLevel>1){
    G4cout << " G4ProcessTable::SetProcessActivation:" ;
    G4cout << " The ProcessType[" << G4int(processType) << "] "<< G4endl;
  }
#endif

  G4ProcTableVector* pTblVector =  Find(fProcTblVector, processType);
  G4ProcTableVector::iterator itr; 
  G4ProcTblElement* anElement;  
  // loop over all elements
  for (itr=pTblVector->begin(); itr!=pTblVector->end(); ++itr) {
    anElement = (*itr);
    G4VProcess* process = anElement->GetProcess();
#ifdef G4VERBOSE
    if (verboseLevel>1){
      G4cout << " The Process[" << process->GetProcessName()<< "] "<< G4endl;
    }
#endif
    for (G4int idx = 0 ; idx < anElement->Length(); idx++) {
      G4ProcessManager* manager = anElement->GetProcessManager(idx);
      manager->SetProcessActivation(process, fActive);
#ifdef G4VERBOSE
      if (verboseLevel>1){
        G4cout << "  for " << manager->GetParticleType()->GetParticleName();
	G4cout << "  Index = " << manager->GetProcessIndex(process) << G4endl;
      }
#endif
    }
  }
}

///////////////
void G4ProcessTable::SetProcessActivation( 
			   G4ProcessType   processType, 
		           G4ProcessManager* processManager, 
			   G4bool          fActive  )
{
#ifdef G4VERBOSE
  if (verboseLevel>1){
    G4cout << " G4ProcessTable::SetProcessActivation:" ;
    G4cout << " The ProcessType[" << G4int(processType) << "] "<< G4endl;
  }
#endif
  
  G4ProcessVector* procList =  processManager->GetProcessList();
  for (G4int idx = 0; idx < procList->length(); idx++) {
    G4VProcess* process = (*procList)(idx);
    if ( process->GetProcessType() == processType) {
      processManager->SetProcessActivation(process, fActive);
#ifdef G4VERBOSE
      if (verboseLevel>1){
        G4cout << " The Process[" << process->GetProcessName()<< "] "<< G4endl;
	G4cout << "  for " << processManager->GetParticleType()->GetParticleName();
	G4cout << "  Index = " << idx << G4endl;
      }
#endif
    }
  }
}


/////////////
void G4ProcessTable::DumpInfo(G4VProcess* process, 
			      G4ParticleDefinition* particle)
{
  G4ProcTableVector::iterator itr; 
  G4int idxTbl=0;
  G4ProcTblElement* anElement =nullptr;
  G4bool isFoundInTbl = false;
  G4ProcessManager* manager =nullptr;
  G4int idx;
  // loop over all elements
  for (itr=fProcTblVector->begin(); itr!=fProcTblVector->end(); ++itr, ++idxTbl) {
    anElement = (*itr);
    if (process == anElement->GetProcess() ){
      if (particle != nullptr) {
	for (idx=0; idx<anElement->Length(); idx++){
	  manager = anElement->GetProcessManager(idx);
	  if (particle == manager->GetParticleType()) {
	    isFoundInTbl = true;
	    break;
	  }
	}
      } else {
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
  if (particle == nullptr) {
    for (idx=0; idx<anElement->Length(); idx++){
      manager = anElement->GetProcessManager(idx);
      G4cout << " for " << manager->GetParticleType()->GetParticleName();
      G4cout << G4endl;
#ifdef G4VERBOSE
      if (verboseLevel >2){
	tmpVerbose = manager->GetVerboseLevel();
	manager->SetVerboseLevel(verboseLevel);
	manager->DumpInfo();
	manager->SetVerboseLevel(tmpVerbose);
      }
#endif
    }
  } else {
    G4cout << " for " << manager->GetParticleType()->GetParticleName();
    G4cout << G4endl;
#ifdef G4VERBOSE
    if (verboseLevel >2){
      tmpVerbose = manager->GetVerboseLevel();
      manager->SetVerboseLevel(verboseLevel);
      manager->DumpInfo();
      manager->SetVerboseLevel(tmpVerbose);
    }
#endif
  }
}






