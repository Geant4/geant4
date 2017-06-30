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
// $Id: G4VModularPhysicsList.cc 103803 2017-04-27 14:03:05Z gcosmo $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class implementation file 
// ------------------------------------------------------------
// - Add  ReplacePhysics             14 Mar 2011 by H.Kurashige
// - Add  RemovePhysics               2 May 2011 by H.Kurashige
//  
#include "G4VModularPhysicsList.hh"
#include "G4StateManager.hh"
#include <algorithm>

// This macros change the references to fields that are now encapsulated
// in the class G4VMPLData.
#define G4MT_physicsVector ((G4VMPLsubInstanceManager.offset[g4vmplInstanceID]).physicsVector)

G4VMPLManager G4VModularPhysicsList::G4VMPLsubInstanceManager;

void G4VMPLData::initialize()
{
    physicsVector = new G4PhysConstVectorData();
}

//G4ThreadLocal G4VModularPhysicsList::G4PhysConstVector* G4VModularPhysicsList::physicsVector = 0;

G4VModularPhysicsList::G4VModularPhysicsList()
  : G4VUserPhysicsList(),
    verboseLevel(0)
{
    g4vmplInstanceID = G4VMPLsubInstanceManager.CreateSubInstance();
}

G4VModularPhysicsList::~G4VModularPhysicsList()
{
  for (auto itr = G4MT_physicsVector->begin(); itr!= G4MT_physicsVector->end(); ++itr) {
    delete (*itr);
  }
  G4MT_physicsVector->clear();
  delete G4MT_physicsVector;
}

G4VModularPhysicsList::G4VModularPhysicsList(const G4VModularPhysicsList& right)  : G4VUserPhysicsList(right),
    verboseLevel(0)
{
    g4vmplInstanceID = G4VMPLsubInstanceManager.CreateSubInstance();
}

G4VModularPhysicsList & G4VModularPhysicsList::operator=(const G4VModularPhysicsList& right)
{
   if (this != &right) {
    defaultCutValue = right.defaultCutValue;
    isSetDefaultCutValue = right.isSetDefaultCutValue;
    fRetrievePhysicsTable = right.fRetrievePhysicsTable;
    fStoredInAscii = right.fStoredInAscii;
    fIsCheckedForRetrievePhysicsTable = right.fIsCheckedForRetrievePhysicsTable;
    fIsRestoredCutValues = right.fIsRestoredCutValues;
    directoryPhysicsTable = right.directoryPhysicsTable;
    //fDisplayThreshold = static_cast<const G4VUserPhysicsList&>(right).GetSubInstanceManager().offset[right.GetInstanceID()]._fDisplayThreshold;
    (this->subInstanceManager.offset[this->g4vuplInstanceID])._fDisplayThreshold=
    		static_cast<const G4VUserPhysicsList&>(right).GetSubInstanceManager().offset[right.GetInstanceID()]._fDisplayThreshold;
    //fIsPhysicsTableBuilt = static_cast<const G4VUserPhysicsList&>(right).GetSubInstanceManager().offset[right.GetInstanceID()]._fIsPhysicsTableBuilt;
    (this->subInstanceManager.offset[this->g4vuplInstanceID])._fDisplayThreshold=
    		static_cast<const G4VUserPhysicsList&>(right).GetSubInstanceManager().offset[right.GetInstanceID()]._fIsPhysicsTableBuilt;
    //fDisplayThreshold = right.fDisplayThreshold;
    fDisableCheckParticleList = right.fDisableCheckParticleList;
    verboseLevel = right.verboseLevel;
    
    if(G4MT_physicsVector !=0) {
      for (auto itr = G4MT_physicsVector->begin(); itr!= G4MT_physicsVector->end(); ++itr) {
	delete (*itr);
      }
      G4MT_physicsVector->clear();
      delete G4MT_physicsVector;
    }
    g4vmplInstanceID = G4VMPLsubInstanceManager.CreateSubInstance();
  }
  return *this;
}

void G4VModularPhysicsList::ConstructParticle()
{
  // create particles
  for (auto itr = G4MT_physicsVector->begin(); itr!= G4MT_physicsVector->end(); ++itr) {
    (*itr)->ConstructParticle();;
  }
}


//Andrea Dotti: May 6 2013
//Current limitation being debugged: Construction of physics processes
//needs to be sequential (there is at least one HAD processes creating problems)
//This is not yet understood and needs to be debugged since we do not want
//this part to be sequential (imagine when one has 100 threads)
//TODO: Remove this lock
#include "G4AutoLock.hh"
namespace {
    G4Mutex constructProcessMutex = G4MUTEX_INITIALIZER;
}

void G4VModularPhysicsList::ConstructProcess()
{
    G4AutoLock l(&constructProcessMutex); //Protection to be removed (A.Dotti)
 AddTransportation();
    
 for (auto itr = G4MT_physicsVector->begin(); itr!= G4MT_physicsVector->end(); ++itr) {
    (*itr)->ConstructProcess();
 }
}



void G4VModularPhysicsList::RegisterPhysics(G4VPhysicsConstructor* fPhysics)
{
  G4StateManager*    stateManager = G4StateManager::GetStateManager();
  G4ApplicationState currentState = stateManager->GetCurrentState();
  if(!(currentState==G4State_PreInit)){
    G4Exception("G4VModularPhysicsList::RegisterPhysics",
		"Run0201", JustWarning,
                "Geant4 kernel is not PreInit state : Method ignored.");
    return;
  }

  G4String pName = fPhysics->GetPhysicsName();  
  G4int pType = fPhysics->GetPhysicsType(); 
  // If physics_type is equal to 0, 
  // following duplication check is omitted 
  // This is TEMPORAL treatment.
  if (pType == 0) {
    G4MT_physicsVector->push_back(fPhysics);
#ifdef G4VERBOSE
    if (verboseLevel >1){
      G4cout << "G4VModularPhysicsList::RegisterPhysics: "
	     <<  pName  << "with type : " << pType 
	     << "  is added" 
	     << G4endl;
    }
#endif
    return;
  }
   
  // Check if physics with the physics_type same as one of given physics 
  auto itr = G4MT_physicsVector->begin();
  for (; itr!= G4MT_physicsVector->end(); ++itr) {
    if ( pType == (*itr)->GetPhysicsType()) break;
  }
  if (itr!= G4MT_physicsVector->end()) {
#ifdef G4VERBOSE
    if (verboseLevel >0){
      G4cout << "G4VModularPhysicsList::RegisterPhysics: "
	     << "a physics with given type already exists " 
	     << G4endl;
      G4cout << " Type = " << pType << " : " 
	     << "  existing physics is " << (*itr)->GetPhysicsName()  
	     << G4endl;
      G4cout << pName << " can not be registered "<<G4endl;
    }
#endif
    G4String comment ="Duplicate type for ";
    comment += pName;
    G4Exception("G4VModularPhysicsList::RegisterPhysics",
		"Run0202", JustWarning, comment); 
    return;
  }

  // register 
  G4MT_physicsVector->push_back(fPhysics);

}    

void G4VModularPhysicsList::ReplacePhysics(G4VPhysicsConstructor* fPhysics)
{
  G4StateManager*    stateManager = G4StateManager::GetStateManager();
  G4ApplicationState currentState = stateManager->GetCurrentState();
  if(!(currentState==G4State_PreInit)){
    G4Exception("G4VModularPhysicsList::ReplacePhysics",
		"Run0203", JustWarning,
                "Geant4 kernel is not PreInit state : Method ignored.");
     return;
  }

  G4String pName = fPhysics->GetPhysicsName();  
  G4int pType = fPhysics->GetPhysicsType(); 
  // If physics_type is equal to 0, 
  // duplication check is omitted and just added. 
  // This is TEMPORAL treatment.
  if (pType == 0) {
    // register 
    G4MT_physicsVector->push_back(fPhysics);
#ifdef G4VERBOSE
    if (verboseLevel >0){
      G4cout << "G4VModularPhysicsList::ReplacePhysics: "
	     <<  pName  << "with type : " << pType 
	     << "  is added" 
	     << G4endl;
    }
#endif
    return;
  }

  // Check if physics with the physics_type same as one of given physics 
  auto itr= G4MT_physicsVector->begin();
  for (itr = G4MT_physicsVector->begin(); itr!= G4MT_physicsVector->end(); ++itr) {
    if ( pType == (*itr)->GetPhysicsType()) break;
  }
  if (itr == G4MT_physicsVector->end()) {
    // register 
    G4MT_physicsVector->push_back(fPhysics);
  } else {
#ifdef G4VERBOSE
    if (verboseLevel >0){
      G4cout << "G4VModularPhysicsList::ReplacePhysics: "
	     << (*itr)->GetPhysicsName()  << "with type : " << pType 
	     << "  is replaces with " << pName
	     << G4endl;
    }
#endif

    //  delete exsiting one   
    delete (*itr);
    // replace with given one
    (*itr) = fPhysics;

  }

  return;  
}    

void  G4VModularPhysicsList::RemovePhysics(G4int pType)
{
  G4StateManager*    stateManager = G4StateManager::GetStateManager();
  G4ApplicationState currentState = stateManager->GetCurrentState();
  if(!(currentState==G4State_PreInit)){
    G4Exception("G4VModularPhysicsList::RemovePhysics",
		"Run0204", JustWarning,
                "Geant4 kernel is not PreInit state : Method ignored.");
    return;
  }

  for (auto itr = G4MT_physicsVector->begin();
      itr!= G4MT_physicsVector->end();) {
      if ( pType  == (*itr)->GetPhysicsType()) {
	  G4String pName = (*itr)->GetPhysicsName();
#ifdef G4VERBOSE
	  if (verboseLevel > 0){
	      G4cout << "G4VModularPhysicsList::RemovePhysics: "
		     <<  pName  << "  is removed"
		     << G4endl;
	  }
#endif
	  G4MT_physicsVector->erase(itr);
	  break;
      } else {
	  itr++;
      }
    }
}

void G4VModularPhysicsList::RemovePhysics(G4VPhysicsConstructor* fPhysics)
{
  G4StateManager*    stateManager = G4StateManager::GetStateManager();
  G4ApplicationState currentState = stateManager->GetCurrentState();
  if(!(currentState==G4State_PreInit)){
    G4Exception("G4VModularPhysicsList::RemovePhysics",
		"Run0205", JustWarning,
                "Geant4 kernel is not PreInit state : Method ignored.");
    return;
  }

  for (auto itr = G4MT_physicsVector->begin();
                                     itr!= G4MT_physicsVector->end();) {
    if ( fPhysics == (*itr)) {
      G4String pName = (*itr)->GetPhysicsName();  
#ifdef G4VERBOSE
    if (verboseLevel > 0 ){
      G4cout << "G4VModularPhysicsList::RemovePhysics: "
	     <<  pName  << "  is removed" 
	     << G4endl;
    }
#endif
      G4MT_physicsVector->erase(itr);
      break;
    } else {
      itr++;
    }
  }
}
void G4VModularPhysicsList::RemovePhysics(const G4String& name)
{
  G4StateManager*    stateManager = G4StateManager::GetStateManager();
  G4ApplicationState currentState = stateManager->GetCurrentState();
  if(!(currentState==G4State_PreInit)){
    G4Exception("G4VModularPhysicsList::RemovePhysics",
		"Run0206", JustWarning,
                "Geant4 kernel is not PreInit state : Method ignored.");
    return;
  }

  for (auto itr = G4MT_physicsVector->begin();
                                     itr!= G4MT_physicsVector->end();) {
    G4String pName = (*itr)->GetPhysicsName();  
    if ( name == pName) {
#ifdef G4VERBOSE
    if (verboseLevel > 0){
      G4cout << "G4VModularPhysicsList::RemovePhysics: "
	     <<  pName  << "  is removed" 
	     << G4endl;
    }
#endif
      G4MT_physicsVector->erase(itr);
      break;
    } else {
      itr++;
    }
  }
}

const G4VPhysicsConstructor* G4VModularPhysicsList::GetPhysics(G4int idx) const
{
  G4int i;
  auto itr= G4MT_physicsVector->begin();
  for (i=0; i<idx && itr!= G4MT_physicsVector->end() ; ++i) ++itr;
  if (itr!= G4MT_physicsVector->end()) return (*itr);
  else return 0;
}

const G4VPhysicsConstructor* G4VModularPhysicsList::GetPhysics(const G4String& name) const
{
  auto itr = G4MT_physicsVector->begin();
  for (; itr!= G4MT_physicsVector->end(); ++itr) {
    if ( name == (*itr)->GetPhysicsName()) break;
  }
  if (itr!= G4MT_physicsVector->end()) return (*itr);
  else return 0;
}

const G4VPhysicsConstructor* G4VModularPhysicsList::GetPhysicsWithType(G4int pType) const
{
  auto itr = G4MT_physicsVector->begin();
  for (; itr!= G4MT_physicsVector->end(); ++itr) {
    if ( pType == (*itr)->GetPhysicsType()) break;
  }
  if (itr!= G4MT_physicsVector->end()) return (*itr);
  else return 0;
}


void G4VModularPhysicsList::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
  // Loop over constructors
  for (auto itr = G4MT_physicsVector->begin(); itr!= G4MT_physicsVector->end(); ++itr) {
    (*itr)->SetVerboseLevel(verboseLevel);
  }

}

void G4VModularPhysicsList::TerminateWorker()
{
  //See https://jira-geant4.kek.jp/browse/DEV-284
  std::for_each( G4MT_physicsVector->begin() , G4MT_physicsVector->end() ,
		 [](G4PhysConstVector::value_type el) { el->TerminateWorker();});
  G4VUserPhysicsList::TerminateWorker();
}
