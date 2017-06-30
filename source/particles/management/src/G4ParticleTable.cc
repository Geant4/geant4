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
// $Id: G4ParticleTable.cc 103108 2017-03-16 13:00:35Z gcosmo $
//
// class G4ParticleTable
//
// Implementation
//
// History:
//      modified                                Apr., 97 H.Kurashige
//      added fParticleMessenger             14 Nov., 97 H.Kurashige
//      added GetParticle()                  13 Dec., 97 H.Kurashige
//      added IonTable and ShortLivedTable   27 June, 98 H.Kurashige 
//      modified FindIon                     02 Aug., 98 H.Kurashige
//      added dictionary for encoding    24 Sep., 98 H.Kurashige
//      fixed bugs in destruction of IonTable 08 Nov.,98 H.Kurashige
//      commented out G4cout/G4cout in the constructor 10 Nov.,98 H.Kurashige
//         --------------------------------
//      modified destructor for STL interface 18 May 1999
//      fixed  some improper codings     08 Apr., 99 H.Kurashige
//      modified FindIon/GetIon methods  17 AUg., 99 H.Kurashige
//      implement new version for using STL map instaed of 
//      RW PtrHashedDictionary           28 Oct., 99  H.Kurashige
//      remove G4ShortLivedTable         25 July, 13 H.Kurashige
//      remove FindIon/GetIon            25 Sep. 14 H.Kurashige
// 

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleTable.hh"
#include "G4UImessenger.hh"
#include "G4ParticleMessenger.hh"
#include "G4IonTable.hh"
#include "G4StateManager.hh"

// These fields should be thread local or thread private. For a singleton
// class, we can change any member field as static without any problem
// because there is only one instance. Then we are allowed to add 
// "G4ThreadLocal".
//
G4ThreadLocal G4ParticleMessenger* G4ParticleTable::fParticleMessenger = 0;
G4ThreadLocal G4ParticleTable::G4PTblDictionary*  G4ParticleTable::fDictionary = 0;
G4ThreadLocal G4ParticleTable::G4PTblDicIterator* G4ParticleTable::fIterator = 0;
G4ThreadLocal G4ParticleTable::G4PTblEncodingDictionary* G4ParticleTable::fEncodingDictionary = 0;

// This field should be thread private. However, we have to keep one copy
// of the ion table pointer. So we change all important fields of G4IonTable
// to the thread local variable.
//
G4IonTable*            G4ParticleTable::fIonTable = 0;


// These shadow pointers are used by each worker thread to copy the content
// from the master thread. 
//
G4ParticleMessenger* G4ParticleTable::fParticleMessengerShadow = 0;
G4ParticleTable::G4PTblDictionary*  G4ParticleTable::fDictionaryShadow = 0;
G4ParticleTable::G4PTblDicIterator* G4ParticleTable::fIteratorShadow = 0;
G4ParticleTable::G4PTblEncodingDictionary* G4ParticleTable::fEncodingDictionaryShadow = 0;

// Static class variable: ptr to single instance of class
G4ParticleTable* G4ParticleTable::fgParticleTable =0;

#ifdef G4MULTITHREADED
// Lock for particle table accesses.
//
G4Mutex G4ParticleTable::particleTableMutex = G4MUTEX_INITIALIZER;
G4int G4ParticleTable::lockCount = 0;
#endif 

////////////////////
G4ParticleTable* G4ParticleTable::GetParticleTable()
{
    static G4ParticleTable theParticleTable;
    if (!fgParticleTable){
      fgParticleTable =  &theParticleTable;
    }

    // Here we initialize all thread private data members.
    //
    if (fDictionary == 0) fgParticleTable->WorkerG4ParticleTable();

    return fgParticleTable;
}

////////////////////
G4ParticleTable::G4ParticleTable()
     :verboseLevel(1),
      noName(" "),
      readyToUse(false),
      genericIon(NULL)
{
  fDictionary = new G4PTblDictionary();

  // Set up the shadow pointer used by worker threads.
  //
  if (fDictionaryShadow == 0)
  {
    fDictionaryShadow = fDictionary;
  }

  fIterator   = new G4PTblDicIterator( *fDictionary );

  // Set up the shadow pointer used by worker threads.
  //
  if (fIteratorShadow == 0)
  {
    fIteratorShadow = fIterator;
  }
 
  fEncodingDictionary = new G4PTblEncodingDictionary();
  // Set up the shadow pointer used by worker threads.
  //
  if (fEncodingDictionaryShadow == 0)
  {
    fEncodingDictionaryShadow = fEncodingDictionary;
  }


   // Ion Table
  fIonTable = new G4IonTable();

}

// This method is similar to the constructor. It is used by each worker
// thread to achieve the partial effect as that of the master thread.
// Here we initialize all thread private data members.
//
void G4ParticleTable::SlaveG4ParticleTable()
{
  G4Exception("G4ParticleTable::SlaveG4ParticleTable()","G4MT0000",FatalException,"Obsolete");
}

void G4ParticleTable::WorkerG4ParticleTable()
{
  // The iterator for the shadow particle table is not sharable.
  //
#ifdef G4MULTITHREADED
  G4MUTEXLOCK(&G4ParticleTable::particleTableMutex);
  G4ParticleTable::lockCount++;
#endif
  if(fDictionary == 0) { 
    fDictionary = new G4PTblDictionary();   
  } else { 
    fDictionary->clear(); 
  }

  if(fEncodingDictionary == 0){
    fEncodingDictionary = new G4PTblEncodingDictionary(); 
  } else { 
    fEncodingDictionary->clear(); 
  }

  fIteratorShadow->reset(false);
  while( (*fIteratorShadow)() ) { // Loop checking, 09.08.2015, K.Kurashige
    G4ParticleDefinition* particle = fIteratorShadow->value();
    fDictionary->insert( std::pair<G4String, G4ParticleDefinition*>(GetKey(particle), particle) );
    G4int code = particle->GetPDGEncoding();
    if (code !=0 ) {
      fEncodingDictionary->insert( std::pair<G4int, G4ParticleDefinition*>(code ,particle) );
    }
  }       
  fIterator =  new G4PTblDicIterator( *fDictionary);

#ifdef G4MULTITHREADED
  G4MUTEXUNLOCK(&G4ParticleTable::particleTableMutex);
#endif

  fIonTable->WorkerG4IonTable();

}

////////////////////
G4ParticleTable::~G4ParticleTable()
{
   readyToUse = false;
   
   // remove all items from G4ParticleTable
   RemoveAllParticles();

  //delete Ion Table 
  if (fIonTable!=0) delete fIonTable;
  fIonTable =0;

  // delete dictionary for encoding
  if (fEncodingDictionary!=0){
    fEncodingDictionary -> clear();
    delete fEncodingDictionary;
    fEncodingDictionary =0;
  }

  if(fDictionary){
    if (fIterator!=0 )delete fIterator;
    fIterator =0;

    fDictionary->clear();
    delete fDictionary;
    fDictionary =0;
  }

  if (fParticleMessenger!=0) delete fParticleMessenger;  
  fParticleMessenger =0;

  fgParticleTable =0;

  G4ParticleDefinition::Clean();  // Delete sub-instance static data
}

////////////////////
void G4ParticleTable::DestroyWorkerG4ParticleTable()
{
  //delete Ion Table in worker thread
  if (fIonTable!=0) fIonTable->DestroyWorkerG4IonTable();

  // delete dictionary for encoding
  if (fEncodingDictionary!=0){
    fEncodingDictionary -> clear();
    delete fEncodingDictionary;
    fEncodingDictionary =0;
  }

  if(fDictionary){
    if (fIterator!=0 )delete fIterator;
    fIterator =0;

    fDictionary->clear();
    delete fDictionary;
    fDictionary =0;
  }

  if (fParticleMessenger!=0) delete fParticleMessenger;  
  fParticleMessenger =0;
}

////////////////////
G4ParticleTable::G4ParticleTable(const G4ParticleTable &right)
  :verboseLevel(1),
   noName(" "),
   readyToUse(false)
{
  fParticleMessenger = 0 ;

  G4Exception("G4ParticleTable::G4ParticleTable()",
	      "PART001", FatalException,
	      "Illegal call of copy constructor for G4ParticleTable");    
  fDictionary = new G4PTblDictionary(*(right.fDictionary));
  fIterator   = new G4PTblDicIterator(*fDictionary);
}

////////////////////
G4ParticleTable & G4ParticleTable::operator=(const G4ParticleTable & right)
{
  if (this != &right) {
    G4Exception("G4ParticleTable::G4ParticleTable()",
		"PART001", FatalException,
		"Illegal call of assignment operator for G4ParticleTable");    
    fDictionary = new G4PTblDictionary(*(right.fDictionary));
    fIterator   = new G4PTblDicIterator(*fDictionary);
  }
  return *this;
}

////////////////////
G4UImessenger* G4ParticleTable::CreateMessenger()
{
  if (fParticleMessenger== 0) {
    //UI messenger
    fParticleMessenger = new G4ParticleMessenger(this);
  }
  return fParticleMessenger;
}

////////////////////
void G4ParticleTable::DeleteMessenger()
{
  if (fParticleMessenger!= 0) {
    //UI messenger
    delete fParticleMessenger;
    fParticleMessenger= 0;
  }

}

////////////////////
void G4ParticleTable::DeleteAllParticles()
{
  //set readyToUse false  
  readyToUse = false;

#ifdef G4VERBOSE
  if (verboseLevel>1){
    G4cout << "G4ParticleTable::DeleteAllParticles() " << G4endl;
  }
#endif

  // delete all particles 
  G4PTblDicIterator *piter = fIterator; 
  piter -> reset(false);
  while( (*piter)() ){// Loop checking, 09.08.2015, K.Kurashige
#ifdef G4VERBOSE
    if (verboseLevel>2){
      G4cout << "Delete " << (piter->value())->GetParticleName() 
	     << " " << (piter->value()) << G4endl;
    }
#endif
    delete (piter->value());
  }
  RemoveAllParticles();
}

////////////////////
void G4ParticleTable::RemoveAllParticles()
{
  if (readyToUse) {
    G4Exception("G4ParticleTable::RemoveAllParticle()",
		"PART115", JustWarning,
		"No effects because readyToUse is true.");    
    return;
  }
  
#ifdef G4VERBOSE
  if (verboseLevel>1){
    G4cout << "G4ParticleTable::RemoveAllParticles() " << G4endl;
  }
#endif

  //remove all contnts in Ion Table
  if (fIonTable!=0) {
    fIonTable->clear();
  }

  // clear dictionary
  if (fDictionary) {
    fDictionary->clear();
  }
}

////////////////////
G4ParticleDefinition* G4ParticleTable::Insert(G4ParticleDefinition *particle)
{

  // check particle name
  if ((particle == 0) || (GetKey(particle).isNull())) {
    G4Exception("G4ParticleTable::Insert()",
		"PART121", FatalException,
		"Particle witnout name can not be registered.");    
#ifdef G4VERBOSE
    if (verboseLevel>1){
      G4cout << "The particle[Addr:" << particle << "] has no name "<< G4endl;
    }
#endif
    return 0;

  }else {  

    if (contains(particle)) {  
#ifdef G4VERBOSE
      if (verboseLevel>2){
        FindParticle(particle) -> DumpTable();
      }
#endif
      G4String msg = "The particle ";
      msg += particle->GetParticleName();
      msg += "  has already been registered in the Particle Table ";
      G4Exception("G4ParticleTable::Insert()",
		  "PART122", FatalException,msg);	
////////////////////      return  FindParticle(particle);
      return particle;

    } else {
      G4PTblDictionary *pdic =  fDictionaryShadow;

      // insert into Dictionary
      pdic->insert( std::pair<G4String, G4ParticleDefinition*>(GetKey(particle), particle) );
#ifdef G4MULTITHREADED
      if(G4Threading::IsWorkerThread())
      { fDictionary->insert( std::pair<G4String, G4ParticleDefinition*>(GetKey(particle), particle) ); }
#endif

      G4PTblEncodingDictionary *pedic =  fEncodingDictionaryShadow;  
      // insert into EncodingDictionary
      G4int code = particle->GetPDGEncoding();
      if (code !=0 ) {
        pedic->insert( std::pair<G4int, G4ParticleDefinition*>(code ,particle) );
#ifdef G4MULTITHREADED
        if(G4Threading::IsWorkerThread())
        { fEncodingDictionary->insert( std::pair<G4int, G4ParticleDefinition*>(code ,particle) ); }
#endif
      }       

      // insert it in IonTable if "nucleus"
      if (fIonTable->IsIon(particle) ){
        fIonTable->Insert(particle);
      }

      // set Verbose Level same as ParticleTable
      particle->SetVerboseLevel(verboseLevel);

#ifdef G4VERBOSE
      if (verboseLevel>3){
        G4cout << "The particle "<< particle->GetParticleName() 
	       << " is inserted in the ParticleTable " << G4endl;
      }
#endif

      return particle;
    }
  }
}

////////////////////
G4ParticleDefinition* G4ParticleTable::Remove(G4ParticleDefinition* particle)
{
  if(!particle) return 0;
#ifdef G4MULTITHREADED
  if(G4Threading::IsWorkerThread()) {
    G4ExceptionDescription ed;
    ed << "Request of removing " << particle->GetParticleName()
       << " is ignored as it is invoked from a worker thread.";
    G4Exception("G4ParticleTable::Remove()","PART10117",JustWarning,ed);
    return 0;
  }
#endif
  if (readyToUse) {
    G4StateManager* pStateManager = G4StateManager::GetStateManager();
    G4ApplicationState currentState = pStateManager->GetCurrentState();
    if (currentState != G4State_PreInit) {
      G4String msg = "Request of removing ";
      msg += particle->GetParticleName();  
      msg += " has No effects other than Pre_Init";
      G4Exception("G4ParticleTable::Remove()",
		"PART117", JustWarning, msg);
      return 0;
    } else {
#ifdef G4VERBOSE
      if (verboseLevel>0){
	G4cout << particle->GetParticleName()
	       << " will be removed from the ParticleTable " << G4endl;
      }
#endif
    }
  }
  
  G4PTblDictionary::iterator it =  fDictionaryShadow->find(GetKey(particle));
  if (it != fDictionaryShadow->end()) {
    fDictionaryShadow->erase(it);
    // remove from EncodingDictionary
    G4int code = particle->GetPDGEncoding();
    if (code !=0 ) {
      fEncodingDictionaryShadow->erase(fEncodingDictionaryShadow->find(code)); 
    }
  } else {
    return 0;
  }

  // remove it from IonTable if "nucleus"
  if (fIonTable->IsIon(particle) ){
    fIonTable->Remove(particle);
  }
  
#ifdef G4VERBOSE
  if (verboseLevel>3){
    G4cout << "The particle "<< particle->GetParticleName()
           << " is removed from the ParticleTable " << G4endl;
  }
#endif

  return particle;
}


////////////////////
G4ParticleDefinition* G4ParticleTable::GetParticle(G4int index) const
{
   CheckReadiness();
  if ( (index >=0) && (index < entries()) ) {
    G4PTblDicIterator *piter = fIterator; 
    piter -> reset(false);
    G4int counter = 0;
    while( (*piter)() ){ // Loop checking, 09.08.2015, K.Kurashige
      if ( counter == index ) return piter->value();
      counter++;
    }
  } 
#ifdef G4VERBOSE
  if (verboseLevel>1){
    G4cout << " G4ParticleTable::GetParticle"
           << " invalid index (=" << index << ")" << G4endl;
  }
#endif
  return 0;
}

////////////////////
const G4String& G4ParticleTable::GetParticleName(G4int index) const
{
  G4ParticleDefinition* aParticle =GetParticle(index);
  if (aParticle != 0) {
    return aParticle->GetParticleName();
  } else {
    return noName;
  }
}

////////////////////
G4ParticleDefinition* G4ParticleTable::FindParticle(const G4String &particle_name) 
{
  G4PTblDictionary::iterator it =  fDictionary->find(particle_name);
  if (it != fDictionary->end()) {
    return (*it).second;
  } else {
#ifdef G4MULTITHREADED
    G4ParticleDefinition* ptcl = 0;
    if(G4Threading::IsWorkerThread())
    {
      G4MUTEXLOCK(&G4ParticleTable::particleTableMutex);
      G4PTblDictionary::iterator its = fDictionaryShadow->find(particle_name);
      if(its != fDictionaryShadow->end())
      {
        fDictionary->insert(*its);
        ptcl = (*its).second;
        G4int code = ptcl->GetPDGEncoding();
        if(code!=0) fEncodingDictionary->insert(std::pair<G4int, G4ParticleDefinition*>(code,ptcl) );
      }
      G4MUTEXUNLOCK(&G4ParticleTable::particleTableMutex);
    }
    return ptcl;
#else
    return 0;
#endif
  }
}

////////////////////
G4ParticleDefinition* G4ParticleTable::FindParticle(const G4ParticleDefinition *particle )
{
  CheckReadiness();
  G4String key = GetKey(particle);
  return FindParticle(key);
}

////////////////////
G4ParticleDefinition* G4ParticleTable::FindParticle(G4int aPDGEncoding ) 
{
   CheckReadiness();
    // check aPDGEncoding is valid
    if (aPDGEncoding == 0){ 
#ifdef G4VERBOSE
      if (verboseLevel>1){
        G4cout << "PDGEncoding  [" <<  aPDGEncoding << "] is not valid " << G4endl;
      }
#endif
      return 0;
    }
    
    G4PTblEncodingDictionary *pedic =  fEncodingDictionary;
    G4ParticleDefinition* particle =0;  

    G4PTblEncodingDictionary::iterator it =  pedic->find(aPDGEncoding );
    if (it != pedic->end()) {
      particle = (*it).second;
    }

#ifdef G4MULTITHREADED
    if(particle == 0 && G4Threading::IsWorkerThread())
    {
      G4MUTEXLOCK(&G4ParticleTable::particleTableMutex);
      G4PTblEncodingDictionary::iterator its = fEncodingDictionaryShadow->find(aPDGEncoding);
      if(its!=fEncodingDictionaryShadow->end())
      {
        particle = (*its).second;
        fEncodingDictionary->insert(*its);
        G4String key = GetKey(particle);
        fDictionary->insert( std::pair<G4String, G4ParticleDefinition*>(key,particle) );
      }
      G4MUTEXUNLOCK(&G4ParticleTable::particleTableMutex);  
    }
#endif

#ifdef G4VERBOSE
    if ((particle == 0) && (verboseLevel>1) ){
      G4cout << "CODE:" << aPDGEncoding << " does not exist in ParticleTable " << G4endl;
    }
#endif
    return particle;
}

////////////////////
void G4ParticleTable::DumpTable(const G4String &particle_name) 
{
  CheckReadiness();
  if (( particle_name == "ALL" ) || (particle_name == "all")){
    // dump all particles 
    G4PTblDicIterator *piter = fIterator; 
    piter -> reset();
    while( (*piter)() ){// Loop checking, 09.08.2015, K.Kurashige
      (piter->value())->DumpTable();
    }
  } else {
    // dump only particle with name of  particle_name
    G4ParticleDefinition *ptr;
    ptr = FindParticle(particle_name);
    if ( ptr != 0) {
      ptr->DumpTable();
    } else { 
#ifdef G4VERBOSE
      if (verboseLevel>1) {
	G4cout << " G4ParticleTable::DumpTable : " 
	       << particle_name << " does not exist in ParticleTable " <<G4endl;
      }
#endif
    }
  }
}

void G4ParticleTable::CheckReadiness() const
{
  if(!readyToUse) {
   G4String msg;
   msg = "Illegal use of G4ParticleTable : ";
   msg += " Access to G4ParticleTable for finding a particle or equivalent\n";
   msg += "operation occurs before G4VUserPhysicsList is instantiated and\n";
   msg += "assigned to G4RunManager. Such an access is prohibited by\n";
   msg += "Geant4 version 8.0. To fix this problem, please make sure that\n";
   msg += "your main() instantiates G4VUserPhysicsList and set it to\n";
   msg += "G4RunManager before instantiating other user classes such as\n";
   msg += "G4VUserPrimaryParticleGeneratorAction.";
   G4Exception("G4ParticleTable::CheckReadiness()",
              "PART002",FatalException,msg);
  }
}

 

G4IonTable*  G4ParticleTable::GetIonTable() const
{
  return fIonTable;
}
 
const G4ParticleTable::G4PTblDictionary* G4ParticleTable::GetDictionary() const
{
  return fDictionary;
}

G4ParticleTable::G4PTblDicIterator* G4ParticleTable::GetIterator() const
{
  return fIterator;
}

const G4ParticleTable::G4PTblEncodingDictionary* G4ParticleTable::GetEncodingDictionary() const
{
  return fEncodingDictionary;
}

G4bool  G4ParticleTable::contains(const G4String& particle_name) const
{
  G4PTblDictionary::iterator it =  fDictionaryShadow->find(particle_name);
  return (it != fDictionaryShadow->end());
}

G4int G4ParticleTable::entries() const
{
  return fDictionary->size();
}

G4int G4ParticleTable::size() const
{
  return fDictionary->size();
}


