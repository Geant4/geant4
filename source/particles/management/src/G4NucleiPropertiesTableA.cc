//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4NucleiPropertiesTableA.cc,v 1.6 2001-07-11 10:02:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1997
//      CERN Geneva Switzerland
//
//
//      File name:     G4NucleiPropertiesTable.cc 
//
//      Authors:       Vicente Lara (Vicente.Lara@cern.ch)
//                     Christian V"olcker (Christian.Volcker@cern.ch),
//
//      Creation date: November 1997
//
//      Modifications: 
// -------------------------------------------------------------------
// Migrate into particles category by H.Kurashige (17 Nov. 98)


#include "G4NucleiPropertiesTable.hh"


// Class G4NucleiPropertiesTable 


G4int G4NucleiPropertiesTable::shortTable[G4NucleiPropertiesTable::MaxA+1];
// make sure there is exactly one instance of this class!
G4NucleiPropertiesTable G4NucleiPropertiesTable::theInstance(2.5);

// Default constructor

//G4NucleiPropertiesTable::G4NucleiPropertiesTable() : G4NucleiPropertiesTable::nEntries(2971) 
G4NucleiPropertiesTable::G4NucleiPropertiesTable(G4double dummy)
{
  
  G4double even_more_dummy = 0;
  even_more_dummy = dummy;
  G4int j = 0;
  
  for (G4int i = 0; i < G4NucleiPropertiesTable::nEntries; i++) {
      if (indexArray[1][i] > j) {
				shortTable[j++] = i;
      }
  }
  shortTable[G4NucleiPropertiesTable::MaxA] = G4NucleiPropertiesTable::nEntries;
}



// Determine the table index for a Nuclide with Z protons and A nucleons

G4int G4NucleiPropertiesTable::GetIndex(G4int Z, G4int A) 
{

  if(A>G4NucleiPropertiesTable::MaxA) G4Exception(
 "G4NucleiPropertiesTable::GetIndex : Nucleon number larger than 273!\n"); 
  if(A<1) G4Exception(
 "G4NucleiPropertiesTable::GetIndex : Nucleon number negativ!\n"); 
  if(Z>A) G4Exception(
 "G4NucleiPropertiesTable::GetIndex : Nucleon number smaller than Z!\n"); 
  
  G4int i = shortTable[A-1];
  while ( i < shortTable[A] ) {
    if (indexArray[0][i] != Z ) i++;
    else return i;
  }
  
  return -1;
}






G4int G4NucleiPropertiesTable::MinZ(G4int A)
{
  G4int i = shortTable[A-1];
  return indexArray[0][i];
}


G4int G4NucleiPropertiesTable::MaxZ(G4int A)
{
  G4int i = shortTable[A]-1;
  return indexArray[0][i];
}




G4double G4NucleiPropertiesTable::GetNuclearMass(G4int Z, G4int A)
{
	G4int i=GetIndex(Z, A);	
	if (i >= 0){
		const G4double NuclearMass = GetAtomicMass(Z,A) - G4double(Z)*electron_mass_c2 +
											 1.433e-5*MeV*pow(G4double(Z),2.39);
		return NuclearMass;
	} else { 
		return 0.0;
	}
}



