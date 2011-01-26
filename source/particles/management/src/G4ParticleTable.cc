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
// $Id: G4ParticleTable.cc,v 1.38 2010-12-22 07:07:59 kurasige Exp $
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
//      implement new version for using STL map instaed of 
//      RW PtrHashedDictionary           28 ct., 99  H.Kurashige


#include "G4ios.hh"
#include "globals.hh"
#include "G4ParticleTable.hh"
#include "G4UImessenger.hh"
#include "G4ParticleMessenger.hh"
#include "G4IonTable.hh"
#include "G4ShortLivedTable.hh"

////////////////////
G4ParticleTable::G4ParticleTable()
     :verboseLevel(1),fParticleMessenger(0),
      noName(" "),
      readyToUse(false)
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
  
   // remove all items from G4ParticleTable
   RemoveAllParticles();

  // delete Short Lived table 
  if (fShortLivedTable!=0) delete fShortLivedTable;
  fShortLivedTable =0;


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

}

////////////////////
G4ParticleTable::G4ParticleTable(const G4ParticleTable &right)
  :verboseLevel(1),fParticleMessenger(0),
   noName(" "),
   readyToUse(false)
{
  G4Exception("G4ParticleTable::G4ParticleTable()",
	      "PART001", FatalException,
	      "Illegal call of copy constructor for G4ParticleTable");    
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
  }

}

////////////////////
void G4ParticleTable::DeleteAllParticles()
{

#ifdef G4VERBOSE
  if (verboseLevel>1){
    G4cout << "G4ParticleTable::DeleteAllParticles() " << G4endl;
  }
#endif

  // delete all particles 
  G4PTblDicIterator *piter = fIterator; 
  piter -> reset();
  while( (*piter)() ){
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

#ifdef G4VERBOSE
  if (verboseLevel>1){
    G4cout << "G4ParticleTable::RemoveAllParticles() " << G4endl;
  }
#endif

  //remove all contnts in Ion Table
  if (fIonTable!=0) {
    fIonTable->clear();
  }

  // remomve all contents in hort Lived table 
  if (fShortLivedTable!=0) {
    fShortLivedTable->clear();
  }

  // clear dictionary for encoding
  if (fEncodingDictionary) {
      fEncodingDictionary->clear();
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
#ifdef G4VERBOSE
    if (verboseLevel>0){
      G4cerr << "The particle[Addr:" << particle << "] has no name "<< G4endl;
    }
#endif
    return 0;

  }else {  

    if (contains(particle)) {
#ifdef G4VERBOSE
      if (verboseLevel>0){
	G4cerr << "The particle " << particle->GetParticleName() 
	       << "has been already registered in the Particle Table "<< G4endl;
      }
      if (verboseLevel>1){
        FindParticle(particle) -> DumpTable();
      }
#endif
      return  FindParticle(particle);

    } else {
      G4PTblDictionary *pdic =  fDictionary;
      G4PTblEncodingDictionary *pedic =  fEncodingDictionary;  

      // insert into Dictionary
      pdic->insert( std::pair<G4String, G4ParticleDefinition*>(GetKey(particle), particle) );

      // insert into EncodingDictionary
      G4int code = particle->GetPDGEncoding();
      if (code !=0 ) {
        pedic->insert( std::pair<G4int, G4ParticleDefinition*>(code ,particle) );
      }       

      // insert it in IonTable if "nucleus"
      if (fIonTable->IsIon(particle) ){
        fIonTable->Insert(particle);
      }

      // insert it in ShortLivedTable if "shortlived"
      if (particle->IsShortLived() ){
	fShortLivedTable->Insert(particle);
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

#ifdef G4VERBOSE
  if (verboseLevel>3){
    G4cout << "The particle "<< particle->GetParticleName()
           << " is removed from the ParticleTable " << G4endl;
  }
#endif

  return particle;
}

////////////////////
G4ParticleDefinition* G4ParticleTable::FindIon(G4int Z, G4int A, G4int , G4int )
{
   CheckReadiness();
   if (Z<=0) return 0;
   if (A<Z) return 0;
   return fIonTable->GetIon(Z, A);
}

////////////////////
G4ParticleDefinition* G4ParticleTable::GetIon(G4int Z, G4int A, G4double E)
{
   CheckReadiness();
   if (Z<=0) return 0;
   if (A<Z) return 0;
   if (E<0.) return 0;
   return fIonTable->GetIon(Z, A, E);
}

////////////////////
G4ParticleDefinition* G4ParticleTable::GetIon(G4int Z, G4int A, G4int L, G4double E)
{
   CheckReadiness();
   if (Z<=0) return 0;
   if (A-L<Z) return 0;
   if (L<0) return 0; 
   if (E<0.) return 0;
   return fIonTable->GetIon(Z, A, L, E);
}

////////////////////
G4ParticleDefinition* G4ParticleTable::FindIon(G4int Z, G4int A, G4double E)
{
   CheckReadiness();
   if (Z<=0) return 0;
   if (A<Z) return 0;
   if (E<0.) return 0;
   return fIonTable->FindIon(Z, A, E);
}

////////////////////
G4ParticleDefinition* G4ParticleTable::FindIon(G4int Z, G4int A, G4int L, G4double E)
{
   CheckReadiness();
   if (Z<=0) return 0;
   if (A-L<Z) return 0;
   if (L<0) return 0;
   if (E<0.) return 0;
   return fIonTable->FindIon(Z, A, L, E);
}

////////////////////
G4ParticleDefinition* G4ParticleTable::GetParticle(G4int index)
{
   CheckReadiness();
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
  if (verboseLevel>1){
    G4cerr << " G4ParticleTable::GetParticle"
           << " invalid index (=" << index << ")" << G4endl;
  }
#endif
  return 0;
}

////////////////////
G4ParticleDefinition* G4ParticleTable::FindParticle(const G4ParticleDefinition *particle)
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
        G4cerr << "PDGEncoding  [" <<  aPDGEncoding << "] is not valid " << G4endl;
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
    if ((particle == 0) && (verboseLevel>1) ){
      G4cerr << "CODE:" << aPDGEncoding << " does not exist in ParticleTable " << G4endl;
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
      G4cerr << " G4ParticleTable::DumpTable : " 
	     << particle_name << " does not exist in ParticleTable " <<G4endl;
    }
  }
}

void G4ParticleTable::CheckReadiness()
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







