// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NucleiPropertiesTable.hh,v 1.1 1999-01-07 16:10:29 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1997
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      File name:     G4NucleiPropertiesTable.cc 
//
//      Authors:       Vicente Lara (Vicente.Lara@cern.ch)
//                     Christian V"olcker (Christian.Volcker@cern.ch),
//
//      Creation date: November 1997
//
//      Testfiles:
//
//      Modifications: 
//
// Migrate into particles category by H.Kurashige (17 Nov. 98)
// Subtract electron mass by H.Kurashige
// -------------------------------------------------------------------
#include "globals.hh"


#ifndef G4NucleiPropertiesTable_h
#define G4NucleiPropertiesTable_h 1

// Class: G4NucleiPropertiesTable
//	Encapsulates Data from G. Audi and A.H. Wapstra, Nucl.
//	Physics,A595 vol 4 p 409-480,
//	25. Dec. 1995.
//	Name, Z, A,
//	Mass Excess
//	Binding Energy
//	Beta Decay Energy
//	Atomic Mass


class G4NucleiPropertiesTable 
{
private:
  
  // Default constructor - this class should exist only once!
  G4NucleiPropertiesTable(G4double dummy);

  static G4NucleiPropertiesTable theInstance;

public:

  // Destructor (generated)
  ~G4NucleiPropertiesTable() { };

  enum  {nEntries = 2931}; // for SUN 

  // Other Operations 

  // Operation: GetMassExcess
  static inline G4double GetMassExcess(G4int Z, G4int A) 
  {
    G4int i=GetIndex(Z, A);
    if (i >= 0) {
	return MassExcess[i]*keV;
    } else {
        return 0.0;
    }
  }


  // Operation: GetBindingEnergy
  static inline G4double GetBindingEnergy(G4int Z, G4int A)
  {
    G4int i=GetIndex(Z, A);
    if (i >= 0) return BindingEnergy[i]*keV;
    else return 0.0;
  }

  // Operation: GetBetaDecayEnergy
  static inline G4double GetBetaDecayEnergy(G4int Z, G4int A)
  {
    G4int i=GetIndex(Z, A);
    if (i >= 0) return BetaEnergy[i]*keV;
    else return 0.0;
  }

  // Operation: GetAtomicMass .. in Geant4 Energy units!
  //   subtract electron mass  by H.K  
  static inline G4double GetAtomicMass(G4int Z, G4int A)
  {
    G4int i=GetIndex(Z, A);
    if (i >= 0) {
      return AtomicMass[i]*amu_c2-Z*electron_mass_c2;
    } else {
      return 0.0;
    }
  }


  // Operation: GetName
  static inline G4String GetName(G4int Z, G4int A)
  {
    G4int i=GetIndex(Z, A);
    if (i >= 0) return Name[i];
    else return " ";
  }

  // Is the nucleus (A,Z) in table?
  static inline G4bool IsInTable(G4int Z, G4int A)
  {
    if (Z <= A && A >= 1 && A <= 273 && GetIndex(Z, A) >= 0) return true;
    else return false;
  }


  static G4int MaxZ(G4int A);
  static G4int MinZ(G4int A);


private:

  // Operation: GetIndex
  static G4int GetIndex(G4int Z, G4int A);
  

  // Data Members for Class Attributes
  //----------------------------------  


  // Number of Entries in the Table
//  const G4int nEntries;



  // The following arrays are static for allow inicialization.
  // The inicialization is Done in G4NucleiPropertiesTable.cc

  // Mass Excess
  static G4double MassExcess[nEntries];
  
  
  // Binding Energy
  static G4double BindingEnergy[nEntries];

  
  // Beta Decay Energy
  static G4double BetaEnergy[nEntries];

  
  // Atomic Mass 
  static G4double AtomicMass[nEntries];

  
  // Chemical Symbol
  static G4String Name[nEntries];

    
  // Table of Z (number of protons) and A (number of nucleons)
  //        indexArray[0][ ] --> Z
  //        indexArray[1][ ] --> A
  static G4int indexArray[2][nEntries];

  // Reduced Table of A for shorter index search.
  //         The index in this table coincide with A-1
  //         For each A value shortTable[A-1] has the index of the 1st occurrence in
  //         the indexArray[][]
  static G4int shortTable[273];


};


#endif


