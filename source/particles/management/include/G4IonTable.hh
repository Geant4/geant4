// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IonTable.hh,v 1.17 2000-10-20 11:34:45 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, IT Division, ASD group
//	History: first implementation, 
//      based on object model of June 27, 98 H.Kurashige
// ------------------------------------------------------------
//      modified GetIon                 02 Aug., 98 H.Kurashige
//      added Remove()                  06 Nov.,98 H.Kurashige
//      add GetNucleusMass              15 Mar. 99  H.Kurashige
//          -----
//      Modified GetIon methods                  17 Aug. 99 H.Kurashige
//      New design using G4VIsotopeTable          5 Oct. 99 H.Kurashige

#ifndef G4IonTable_h
#define G4IonTable_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"

#include "g4std/vector"

class G4ParticleTable;
class G4VIsotopeTable; 
class G4IsotopeProperty;

class G4IonTable
{
 // Class Description
 //   G4IonTable is the table of pointer to G4ParticleDefinition
 //   In G4IonTable, each G4ParticleDefinition pointer is stored
 //

 public:
   // Use STL Vector as list of ions
   typedef G4std::vector<G4ParticleDefinition*> G4IonList;

 public:
  // constructor
   G4IonTable();

 protected:
   // hide copy construictor as protected 
   G4IonTable(const  G4IonTable &right);
   G4IonTable & operator = (const G4IonTable &) {return *this;}

 public:
  // destructor
   virtual ~G4IonTable();

 public: // With Description
   G4int GetNumberOfElements() const;
   // Get number of elements defined in the IonTable

   // Register Isotope table
   void RegisterIsotopeTable(G4VIsotopeTable* table);
   G4VIsotopeTable* GetIsotopeTable() const;
   // G4IonTable asks properties of isotopes to this G4VIsotopeTable 
   // by using FindIsotope(G4IsotopeProperty* property) method.
   
   // ---------------------------  
   // FindIon/GetIon
   //   FindIon methods return pointer of ion if it exists       
   //   GetIon methods also return pointer of ion. In GetIon 
   //   methods the designated ion will be created if it does not exist.
   //
   // !! PDGCharge inG4ParticleDefinition of ions is           !!
   // !! electric charge of nucleus (i.e. fully ionized ions)  !!
   // -----------------------------

   // Find/Get "ground state" 
   G4ParticleDefinition* GetIon(G4int Z, G4int A, G4int J=0);
   // The ion is assumed to be ground state (i.e Excited energy = 0) 
   //   Z: Atomic Number
   //   A: Atomic Mass
   //   J: Total Angular momentum (in unit of 1/2)

   // Find/Get "excited state" 
   G4ParticleDefinition* FindIon(G4int Z, G4int A, G4double E, G4int J=0);
   G4ParticleDefinition* GetIon(G4int Z, G4int A, G4double E, G4int J=0);
   //   Z: Atomic Number
   //   A: Atomic Mass
   //   J: Total Angular momentum (in unit of 1/2)
   //   E: Excitaion energy

   G4ParticleDefinition* GetIon(G4int Z, G4int A, G4int J, G4int Q);
   // This method is provided for compatibilties 
   // The third and last arguments gives no effect

   static G4bool        IsIon(G4ParticleDefinition*);
   // return true if the particle is ion

   G4String             GetIonName(G4int Z, G4int A, G4double E) const;
   // get ion name

   G4double             GetIonMass(G4int Z, G4int A) const;
   G4double             GetNucleusMass(G4int Z, G4int A) const;
   // These two methods returns Nucleus (i.e. full ionized atom) mass 
   // ,where Z is Atomic Number (number of protons) and
   //  A is Atomic Number (number of nucleons)

  
   G4int                 Entries() const;
   // Return number of ions in the table

   G4ParticleDefinition* GetParticle(G4int index) const;
   // Return the pointer of index-th ion in the table

   G4bool                Contains(const G4ParticleDefinition *particle) const;
   // Return 'true' if the ion exists

   void                  Insert(G4ParticleDefinition* particle);
   void                  Remove(G4ParticleDefinition* particle);
   // Insert/Remove an ion in the table

    void DumpTable(const G4String &particle_name = "ALL") const;
   // dump information of particles specified by name 


 protected:
   G4ParticleDefinition* CreateIon(G4int Z, G4int A, G4double E, G4int J);
   // Create Ion 
   
   G4IsotopeProperty* FindIsotope(G4int Z, G4int A, G4double E, G4int J);
   // Ask properties of isotopes to this G4VIsotopeTable 
   
   G4ParticleDefinition* GetLightIon(G4int Z, G4int A) const;
   
   
   G4bool                IsLightIon(G4ParticleDefinition*) const;
   // return true if the particle is pre-defined ion
 
   void                  AddProcessManager(const G4String& ionName);
   // Add process manager to ions with name of 'ionName'

   void                  SetCuts(G4ParticleDefinition* ion);
   // Set cut value same as "GenericIon" and build physics tables

   G4int                GetVerboseLevel() const;
   // get Verbose Level defined in G4ParticleTable

 private:
   G4IonList*                  fIonList; 

   G4VIsotopeTable*            fIsotopeTable;
 
   enum { numberOfElements = 110};
   static const G4String       elementName[numberOfElements];

};

inline G4int  G4IonTable::GetNumberOfElements() const
{
  return numberOfElements;
}
inline G4bool  G4IonTable::Contains(const G4ParticleDefinition* particle) const
{
  G4IonList::iterator i;
  for (i = fIonList->begin(); i!= fIonList->end(); ++i) {
    if (**i==*particle) return true;
  }
  return false;
}

inline G4int G4IonTable::Entries() const
{
  return fIonList->size();
}

inline G4ParticleDefinition*  G4IonTable::GetParticle(G4int index) const
{
  if ( (index >=0 ) && (index < Entries()) ){
    return (*fIonList)[index];
  } else {
    return 0; 
  } 
}


#endif










