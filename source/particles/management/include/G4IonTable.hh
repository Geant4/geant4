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
// $Id: G4IonTable.hh,v 1.31 2010-10-16 06:04:51 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	History: first implementation, 
//      based on object model of June 27, 98 H.Kurashige
// ------------------------------------------------------------
//      added clear()			20 Mar., 08 H.Kurashige
//      modified GetIon                 02 Aug., 98 H.Kurashige
//      added Remove()                  06 Nov.,98 H.Kurashige
//      add GetNucleusMass              15 Mar. 99  H.Kurashige
//          -----
//      Modified GetIon methods                  17 Aug. 99 H.Kurashige
//      New design using G4VIsotopeTable          5 Oct. 99 H.Kurashige
//      Add GetNucleusEncoding according PDG 2006 9 Oct. 2006 H.Kurashige
//      Use STL map                              30 Jul. 2009 H.Kurashige

#ifndef G4IonTable_h
#define G4IonTable_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Ions.hh"

#include <cmath>
#include <vector>
#include <map>

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
   // Use STL map as list of ions
   typedef  std::multimap<G4int, const G4ParticleDefinition*> G4IonList;
   typedef  std::multimap<G4int, const G4ParticleDefinition*>::iterator G4IonListIterator;

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
   G4VIsotopeTable* GetIsotopeTable(size_t idx=0) const;
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
   G4ParticleDefinition* GetIon(G4int encoding);
   // The ion can be get by using PDG encoding 
   // !! Only ground state can be obtained .i.e. Isomer = 0
   
   void CreateAllIon();
   // All ground state ions will be created

   // Find/Get "excited state" 
   G4ParticleDefinition* FindIon(G4int Z, G4int A, G4double E, G4int J=0);
   G4ParticleDefinition* GetIon(G4int Z, G4int A, G4double E, G4int J=0);
   G4ParticleDefinition* FindIon(G4int Z, G4int A, G4int L,
				 G4double E, G4int J=0);
   G4ParticleDefinition* GetIon(G4int Z, G4int A, G4int L,
				G4double E, G4int J=0);
   //   Z: Atomic Number
   //   A: Atomic Mass (nn + np +nlambda)
   //   L: Number of Lmabda
   //   J: Total Angular momentum (in unit of 1/2)
   //   E: Excitaion energy

   G4ParticleDefinition* GetIon(G4int Z, G4int A, G4int J, G4int Q);
   // This method is provided for compatibilties 
   // The third and last arguments gives no effect

   static G4bool        IsIon(const G4ParticleDefinition*);
   // return true if the particle is ion

   static G4bool        IsAntiIon(const G4ParticleDefinition*);
   // return true if the particle is anti_ion


   const G4String&  GetIonName(G4int Z, G4int A, G4double E) const;
   const G4String&  GetIonName(G4int Z, G4int A, G4int L, G4double E) const;
   // get ion name
  
   static G4int GetNucleusEncoding(G4int Z,        G4int A, 
				   G4double E=0.0, G4int J=0);
  //  get PDG code for Ions 
  // Nuclear codes are given as 10-digit numbers +-100ZZZAAAI.
  //For a nucleus consisting of np protons and nn neutrons
  // A = np + nn and Z = np.
  // I gives the isomer level, with I = 0 corresponding 
  // to the ground state and I >0 to excitations
  //
  //!!! I = 1 is assigned fo all excitation states in Geant4   
  
   static G4int GetNucleusEncoding(G4int Z,        G4int A,  G4int L,        
				   G4double E=0.0, G4int J=0);
  //  get PDG code for Hyper-Nucleus Ions 
  // Nuclear codes are given as 10-digit numbers +-10LZZZAAAI.
  //For a nucleus consisting of np protons and nn neutrons
  // A = np + nn +nlambda and Z = np.
  // L = nlambda
  // I gives the isomer level, with I = 0 corresponding 
  // to the ground state and I >0 to excitations
  //
  //!!! I = 1 is assigned fo all excitation states in Geant4   

   static G4bool GetNucleusByEncoding(G4int encoding,
				     G4int &Z,      G4int &A, 
				     G4double &E,   G4int &J);
   static G4bool GetNucleusByEncoding(G4int encoding,
				      G4int &Z,      G4int &A,  G4int &L,    
				      G4double &E,   G4int &J);
    //!!! Only ground states are supported now  
 
 
   G4double             GetIonMass(G4int Z, G4int A, G4int L=0) const;
   G4double             GetNucleusMass(G4int Z, G4int A, G4int L=0) const;
   // These two methods returns Nucleus (i.e. full ionized atom) mass 
   // ,where Z is Atomic Number (number of protons) and
   //  A is Atomic Number (number of nucleons and hyperons)
   //  L is number of lambda (A= nn + np + nlambda)

 
   
   G4int                 Entries() const;
   // Return number of ions in the table

   G4ParticleDefinition* GetParticle(G4int index) const;
   // Return the pointer of index-th ion in the table
 
   G4bool                Contains(const G4ParticleDefinition *particle) const;
   // Return 'true' if the ion exists

   void                  Insert(const G4ParticleDefinition* particle);
   void                  Remove(const G4ParticleDefinition* particle);
   // Insert/Remove an ion in the table

   void                  clear();
   // erase all contents in the list (not delete just remove)

   G4int                 size() const;
   //  Return number of ions in the table

    void DumpTable(const G4String &particle_name = "ALL") const;
   // dump information of particles specified by name 


 protected:
   G4ParticleDefinition* CreateIon(G4int Z, G4int A, G4double E, G4int J);
   G4ParticleDefinition* CreateIon(G4int Z, G4int A, G4int L, 
				   G4double E, G4int J);
   // Create Ion 
   
   G4IsotopeProperty* FindIsotope(G4int Z, G4int A, G4double E, G4int J);
   // Ask properties of isotopes to this G4VIsotopeTable 
   
   G4ParticleDefinition* GetLightIon(G4int Z, G4int A) const;
   G4ParticleDefinition* GetLightAntiIon(G4int Z, G4int A) const;
   
   G4bool                IsLightIon(const G4ParticleDefinition*) const;
   G4bool                IsLightAntiIon(const G4ParticleDefinition*) const;
   // return true if the particle is pre-defined ion
 
   void                  AddProcessManager(const G4String& ionName);
   // Add process manager to ions with name of 'ionName'

   G4int                GetVerboseLevel() const;
   // get Verbose Level defined in G4ParticleTable

 private:
   G4IonList*                  fIonList; 

   std::vector<G4VIsotopeTable*> *fIsotopeTableList;

 
   enum { numberOfElements = 118};
   static const G4String       elementName[numberOfElements];

};

inline G4int  G4IonTable::GetNumberOfElements() const
{
  return numberOfElements;
}

inline G4bool  G4IonTable::Contains(const G4ParticleDefinition* particle) const
{
  if (!IsIon(particle)) return false;

  G4int Z = particle->GetAtomicNumber();
  G4int A = particle->GetAtomicMass();  
  G4int L = particle->GetQuarkContent(3);  //strangeness 
  G4int encoding=GetNucleusEncoding(Z, A, L);
  G4bool found = false;
  if (encoding !=0 ) {
    G4IonList::iterator i = fIonList->find(encoding);
    for( ;i != fIonList->end() ; i++) {
      if (particle == i->second ) {
	found  = true;
	break;
      }
    }
  }
  return found;
}

inline G4int G4IonTable::Entries() const
{
  return fIonList->size();
}

inline G4int G4IonTable::size() const
{
  return fIonList->size();
}

inline void G4IonTable::clear()
{
#ifdef G4VERBOSE
    if (GetVerboseLevel()>2) {
      G4cout << "G4IonTable::Clear() : number of Ion regsitered =  "; 
      G4cout << fIonList->size() <<  G4endl;
    }
#endif
  fIonList->clear();
}



#endif










