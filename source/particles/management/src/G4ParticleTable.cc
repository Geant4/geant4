// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParticleTable.cc,v 1.1 1999-01-07 16:10:35 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
//      commented out G4cout/G4cerr in the constructor 10 Nov.,98 H.Kurashige

#include "G4ios.hh"
#include "globals.hh"
#include "G4ParticleTable.hh"
#include "G4UImessenger.hh"
#include "G4ParticleMessenger.hh"
#include "G4IonTable.hh"
#include "G4ShortLivedTable.hh"

const G4int G4ParticleTableDefaultBucket = 32;
const G4int G4ParticleTableMaxInABucket = 64;

G4ParticleTable::G4ParticleTable():verboseLevel(0),fParticleMessenger(NULL)
{
  DictionaryBucketSize = G4ParticleTableDefaultBucket;
  fDictionary = new G4PTblDictionary(G4ParticleTable::HashFun,DictionaryBucketSize);
  fIterator   = new G4PTblDicIterator( *fDictionary );
  fEncodingDictionary = new G4PTblEncodingDictionary(G4ParticleTable::EncodingHashFun,DictionaryBucketSize);

 // Ion Table
  fIonTable = new G4IonTable();

  // short lived table
  fShortLivedTable = new G4ShortLivedTable();
}


G4ParticleTable::~G4ParticleTable()
{
  // delete G4String objects for key
  G4String* pS;
  if (fDictionary->entries()>0){
    fIterator -> reset();
    while( (*fIterator)() ){
      if (pS =  fIterator -> key()) delete pS;
    }
    fDictionary -> clear();
  }
  delete fIterator;
  delete fDictionary;
  if (fParticleMessenger!=NULL) delete fParticleMessenger;  

  // delete dictionary for encoding
  if (fEncodingDictionary!=NULL){
    fEncodingDictionary -> clear();
    delete fEncodingDictionary;
  }

  //delete Ion Table and contents
  if (fIonTable!=NULL) delete fIonTable;
 
  // delete Short Lived table and contents
  if (fShortLivedTable!=NULL) delete fShortLivedTable;
}

G4ParticleTable::G4ParticleTable(const G4ParticleTable &right)
{
  G4Exception("you call copy constructor of G4ParticleTable");    
  fDictionary = new G4PTblDictionary(*(right.fDictionary));
  fIterator   = new G4PTblDicIterator(*fDictionary);
}

// Static class variable: ptr to single instance of class
G4ParticleTable* G4ParticleTable::fgParticleTable =0;

G4ParticleTable* G4ParticleTable::GetParticleTable()
{
    static G4ParticleTable theParticleTable;
    if (!fgParticleTable){
      fgParticleTable =  &theParticleTable;
    }
    return fgParticleTable;
}

G4UImessenger* G4ParticleTable::CreateMessenger()
{
  if (fParticleMessenger== NULL) {
    //UI messenger
    fParticleMessenger = new G4ParticleMessenger(this);
  }
  return fParticleMessenger;
}

void G4ParticleTable::DeleteMessenger()
{
  if (fParticleMessenger!= NULL) {
    //UI messenger
    delete fParticleMessenger;
    fParticleMessenger= NULL;
    // remove all items from G4ParticleTable
    RemoveAllParticles();
  }
}

void G4ParticleTable::RemoveAllParticles()
{

#ifdef G4VERBOSE
  if (verboseLevel>1){
    G4cerr << "G4ParticleTable::RemoveAllParticles() " << endl;
  }
#endif

  //delete Ion Table and contents
  if (fIonTable!=NULL) {
    delete fIonTable;
    fIonTable = NULL;
  }

  // delete Short Lived table and contents
  if (fShortLivedTable!=NULL) {
    delete fShortLivedTable;
    fShortLivedTable = NULL;
  }

  // delete dictionary for encoding
  if (fEncodingDictionary!=NULL){
   G4int* pCode;
   G4PTblEncodingDicIterator* fEncodeIterator;
   fEncodeIterator = new  G4PTblEncodingDicIterator(*fEncodingDictionary);
   fEncodeIterator-> reset();
   while( (*fEncodeIterator)() ){
     if (pCode = fEncodeIterator-> key()) delete pCode;
   }
   fEncodingDictionary -> clear();
   delete fEncodeIterator;
   delete fEncodingDictionary;
   fEncodingDictionary = NULL;
  }

  // delete G4String objects for key
  G4String* pS;
  fIterator -> reset();
  while( (*fIterator)() ){
    if (pS =  fIterator -> key()) delete pS;
  }
  fDictionary -> clear();

}

G4ParticleDefinition* G4ParticleTable::Insert(G4ParticleDefinition *particle)
{

  // check particle name
  if ((particle == NULL) || (GetKey(particle).isNull())) {
    //#ifdef G4VERBOSE
    //if (verboseLevel>0){
    //   G4cerr << "The particle[Addr:" << particle << "] has no name "<< endl;
    //}
    //#endif
    return NULL;
  }else {  
    if (contains(particle)) {
      //#ifdef G4VERBOSE
      //if (verboseLevel>1){
      //  G4cerr << "The particle[Addr:" << particle; 
      //  G4cerr << "] has same name as "<< endl;
      //  FindParticle(particle) -> DumpTable();
      //}
      //#endif
      return  FindParticle(particle);
    } else {
      G4PTblDictionary *pdic = GetDictionary();
      pdic -> insertKeyAndValue(new G4String(GetKey(particle)), particle);
      // if HashDictionary get too big, rehash all of the key
      if (pdic -> entries() > G4ParticleTableMaxInABucket*DictionaryBucketSize) {
	DictionaryBucketSize *= 2;
	pdic -> resize(DictionaryBucketSize);
	G4PTblDicIterator *piter = GetIterator(); 
	piter -> reset();
      }

      // insert into EncodingDictionary
      G4int code = particle->GetPDGEncoding();
      if (code !=0 ) {
	G4PTblEncodingDictionary *pedic = GetEncodingDictionary();  
	pedic -> insertKeyAndValue(new G4int(code), particle);
      }       
 
      // insert it in IonTable if "nucleus"
      if (fIonTable->IsIon(particle) ){
	fIonTable->Insert(particle);
      }

      // insert it in SHortLivedTable if "shortlived"
      if (particle->IsShortLived() ){
	fShortLivedTable->Insert(particle);
      }

      return particle;
    }
  }
}

G4ParticleDefinition* G4ParticleTable::Remove(G4ParticleDefinition* particle)
{
  if (! contains(particle) ) return NULL;

  G4String particle_name = GetKey(particle);
  fDictionary->remove(&particle_name );

  // remove from EncodingDictionary
  G4int code = particle->GetPDGEncoding();
  if (code !=0 ) {
    G4PTblEncodingDictionary *pedic = GetEncodingDictionary();  
    pedic -> remove(&code);
  }       
  
  // remove it from IonTable if "nucleus"
  if (fIonTable->IsIon(particle) ){
    fIonTable->Remove(particle);
  }
  
  // Remove it from ShortLivedTable if "shortlived"
  if (particle->IsShortLived() ){
    fShortLivedTable->Remove(particle);
  }
  return particle;
}

G4ParticleDefinition* G4ParticleTable::FindIon(G4int Z, G4int A, G4int J, G4int Q)
{
   if (Z<=0) return NULL;
   if (A<Z) return NULL;
   if (J<0) return NULL;
   if (Q<0) Q=Z;
   return fIonTable->GetIon(Z, A, J, Q);
}

G4ParticleDefinition* G4ParticleTable::GetParticle(G4int index)
{
  if ( (index >=0) && (index < entries()) ) {
    G4PTblDicIterator *piter = GetIterator(); 
    piter -> reset();
    G4int counter = 0;
    while( (*piter)() ){
      if ( counter == index ) return piter->value();
      counter++;
    }
  }
#ifdef G4VERBOSE
  if (verboseLevel>0){
    G4cerr << " G4ParticleTable::GetParticle";
    G4cerr << " invalid index (=" << index << ")" << endl;
  }
#endif
  return NULL;
}

G4ParticleDefinition* G4ParticleTable::FindParticle(G4int aPDGEncoding )  const
{
    // check aPDGEncoding is valid
    if (aPDGEncoding == 0){ 
#ifdef G4VERBOSE
      if (verboseLevel>0){
        G4cerr << "PDGEncoding  [" <<  aPDGEncoding << "] is not valid " << endl;
      }
#endif
      return NULL;
    }

    G4PTblEncodingDictionary *pedic = GetEncodingDictionary();  
    G4ParticleDefinition* particle = pedic -> findValue( &aPDGEncoding );

#ifdef G4VERBOSE
    if ((particle == NULL) && (verboseLevel>0) ){
      G4cerr << "CODE:" << aPDGEncoding << " does not exist in ParticleTable " << endl;
    }
#endif
    return particle;
}

void G4ParticleTable::DumpTable(const G4String &particle_name) const
{
  if (( particle_name == "ALL" ) || (particle_name == "all")){
    // dump all particles 
    G4PTblDicIterator *piter = GetIterator(); 
    piter -> reset();
    while( (*piter)() ){
      (piter->value())->DumpTable();
    }
  } else {
    // dump only particle with name of  particle_name
    G4ParticleDefinition *ptr;
    if ( (ptr = FindParticle(particle_name)) != NULL) {
      ptr->DumpTable();
    } else {
#ifdef G4VERBOSE
      if (verboseLevel>0){
        G4cout << particle_name << " does not exist in ParticleTable " <<endl;
      }
#endif
    }
  }
}








