// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParticleTable.cc,v 1.17 2000-10-20 11:35:57 kurasige Exp $
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
//      commented out G4cout/G4cout in the constructor 10 Nov.,98 H.Kurashige
//         --------------------------------
//      modified destructor for STL interface 18 May 1999
//      fixed  some improper codings     08 Apr., 99 H.Kurashige
//      modified FindIon/GetIon methods  17 AUg., 99 H.Kurashige
//      implement new version for using STL map instaed of RW PtrHashedDictionary
//                                       28 ct., 99  H.Kurashige


#include "G4ios.hh"
#include "globals.hh"
#include "G4ParticleTable.hh"
#include "G4UImessenger.hh"
#include "G4ParticleMessenger.hh"
#include "G4IonTable.hh"
#include "G4ShortLivedTable.hh"


////////////////////
G4ParticleTable::G4ParticleTable():verboseLevel(0),fParticleMessenger(0),noName(" ")
{
  fDictionary = new G4PTblDictionary();
  fIterator   = new G4PTblDicIterator( *fDictionary );
  fEncodingDictionary = new G4PTblEncodingDictionary();

 // Ion Table
  fIonTable = new G4IonTable();

  // short lived table
  fShortLivedTable = new G4ShortLivedTable();
}


////////////////////
G4ParticleTable::~G4ParticleTable()
{
  // delete Short Lived table and contents
  if (fShortLivedTable!=0) delete fShortLivedTable;
  fShortLivedTable =0;


  //delete Ion Table and contents
  if (fIonTable!=0) delete fIonTable;
  fIonTable =0;

  // delete dictionary for encoding
  if (fEncodingDictionary!=0){
    fEncodingDictionary -> clear();
    delete fEncodingDictionary;
    fEncodingDictionary =0;
  }

  if(fDictionary){
    fDictionary->clear();
    delete fDictionary;
    fDictionary =0;
    if (fIterator!=0 )delete fIterator;
    fIterator =0;
  }

  if (fParticleMessenger!=0) delete fParticleMessenger;  
  fParticleMessenger =0;
}

////////////////////
G4ParticleTable::G4ParticleTable(const G4ParticleTable &right)
{
  G4Exception("you call copy constructor of G4ParticleTable");    
  fDictionary = new G4PTblDictionary(*(right.fDictionary));
  fIterator   = new G4PTblDicIterator(*fDictionary);
}



// Static class variable: ptr to single instance of class
G4ParticleTable* G4ParticleTable::fgParticleTable =0;

////////////////////
G4ParticleTable* G4ParticleTable::GetParticleTable()
{
    static G4ParticleTable theParticleTable;
    if (!fgParticleTable){
      fgParticleTable =  &theParticleTable;
    }
    return fgParticleTable;
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
    // remove all items from G4ParticleTable
    RemoveAllParticles();
  }
}

////////////////////
void G4ParticleTable::RemoveAllParticles()
{

#ifdef G4VERBOSE
  if (verboseLevel>1){
    G4cout << "G4ParticleTable::RemoveAllParticles() " << G4endl;
  }
#endif

  //delete Ion Table and contents
  if (fIonTable!=0) {
    delete fIonTable;
    fIonTable = 0;
  }

  // delete Short Lived table and contents
  if (fShortLivedTable!=0) {
    delete fShortLivedTable;
    fShortLivedTable = 0;
  }

  // delete dictionary for encoding
  if (fEncodingDictionary)
    {
      fEncodingDictionary->clear();
      delete fEncodingDictionary;
      fEncodingDictionary = 0;
    }

  // delete dictionary
  if (fDictionary)
    {
      if (fIterator!=0 )delete fIterator;
      fIterator =0;
      fDictionary->clear();
      delete fDictionary;
      fDictionary = 0;
    }
}

////////////////////
G4ParticleDefinition* G4ParticleTable::Insert(G4ParticleDefinition *particle)
{

  // check particle name
  if ((particle == 0) || (GetKey(particle).isNull())) {
#ifdef G4VERBOSE
    if (verboseLevel>0){
      G4cout << "The particle[Addr:" << particle << "] has no name "<< G4endl;
    }
#endif
    return 0;

  }else {  

    if (contains(particle)) {
#ifdef G4VERBOSE
      if (verboseLevel>0){
	G4cout << "The particle has same name "<< G4endl;
      }
      if (verboseLevel>1){
        FindParticle(particle) -> DumpTable();
      }
#endif
      return  FindParticle(particle);

    } else {
      G4PTblDictionary *pdic =  fDictionary;
      G4PTblEncodingDictionary *pedic =  fEncodingDictionary;  

      (*pdic)[GetKey(particle)] = particle;
      // insert into EncodingDictionary
      G4int code = particle->GetPDGEncoding();
      if (code !=0 ) {
       (*pedic)[code] = particle;
      }       

      // insert it in IonTable if "nucleus"
      if (fIonTable->IsIon(particle) ){
        fIonTable->Insert(particle);
      }

      // insert it in ShortLivedTable if "shortlived"
      if (particle->IsShortLived() ){
	fShortLivedTable->Insert(particle);
      }

      return particle;
    }
  }
}


////////////////////
G4ParticleDefinition* G4ParticleTable::Remove(G4ParticleDefinition* particle)
{
  G4PTblDictionary::iterator it =  fDictionary->find(GetKey(particle));
  if (it != fDictionary->end()) {
    fDictionary->erase(it);
    // remove from EncodingDictionary
    G4int code = particle->GetPDGEncoding();
    if (code !=0 ) {
      fEncodingDictionary->erase(fEncodingDictionary->find(code)); 
    }
  } else {
    return 0;
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

////////////////////
G4ParticleDefinition* G4ParticleTable::FindIon(G4int Z, G4int A, G4int J, G4int Q)
{
   if (Z<=0) return 0;
   if (A<Z) return 0;
   return fIonTable->GetIon(Z, A);
}

////////////////////
G4ParticleDefinition* G4ParticleTable::GetIon(G4int Z, G4int A, G4double E)
{
   if (Z<=0) return 0;
   if (A<Z) return 0;
   if (E<0.) return 0;
   return fIonTable->GetIon(Z, A, E);
}

////////////////////
G4ParticleDefinition* G4ParticleTable::FindIon(G4int Z, G4int A, G4double E)
{
   if (Z<=0) return 0;
   if (A<Z) return 0;
   if (E<0.) return 0;
   return fIonTable->FindIon(Z, A, E);
}

////////////////////
G4ParticleDefinition* G4ParticleTable::GetParticle(G4int index)
{
  if ( (index >=0) && (index < entries()) ) {
    G4PTblDicIterator *piter = fIterator; 
    piter -> reset();
    G4int counter = 0;
    while( (*piter)() ){
      if ( counter == index ) return piter->value();
      counter++;
    }
  }
#ifdef G4VERBOSE
  if (verboseLevel>0){
    G4cout << " G4ParticleTable::GetParticle";
    G4cout << " invalid index (=" << index << ")" << G4endl;
  }
#endif
  return 0;
}

////////////////////
G4ParticleDefinition* G4ParticleTable::FindParticle(const G4ParticleDefinition *particle)
{
  G4String key = GetKey(particle);
  return FindParticle(key);
}

////////////////////
G4ParticleDefinition* G4ParticleTable::FindParticle(G4int aPDGEncoding )
{
    // check aPDGEncoding is valid
    if (aPDGEncoding == 0){ 
#ifdef G4VERBOSE
      if (verboseLevel>0){
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

#ifdef G4VERBOSE
    if ((particle == 0) && (verboseLevel>0) ){
      G4cout << "CODE:" << aPDGEncoding << " does not exist in ParticleTable " << G4endl;
    }
#endif
    return particle;
}

////////////////////
void G4ParticleTable::DumpTable(const G4String &particle_name)  
{
  if (( particle_name == "ALL" ) || (particle_name == "all")){
    // dump all particles 
    G4PTblDicIterator *piter = fIterator; 
    piter -> reset();
    while( (*piter)() ){
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
      if (verboseLevel>0){
        G4cout << particle_name << " does not exist in ParticleTable " <<G4endl;
      }
#endif
    }
  }
}








