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
// $Id: G4NucleiPropertiesTheoreticalTableA.cc,v 1.3 2001-07-11 10:02:02 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
// ------------------------------------------------------------

#include "G4NucleiPropertiesTheoreticalTable.hh"



G4int G4NucleiPropertiesTheoreticalTable::shortTable[G4NucleiPropertiesTheoreticalTable::shortTableSize];

G4NucleiPropertiesTheoreticalTable G4NucleiPropertiesTheoreticalTable::theInstance(2.5);
// Default constructor
G4NucleiPropertiesTheoreticalTable::G4NucleiPropertiesTheoreticalTable(G4double dummy)
{
	G4int j = 0;
	
	for (G4int i = 0; i < G4NucleiPropertiesTheoreticalTable::nEntries; i++) {
		if (indexArray[0][i] > j+7) {
			shortTable[j++] = i;
		}
	}
	shortTable[G4NucleiPropertiesTheoreticalTable::shortTableSize-1] = 
										G4NucleiPropertiesTheoreticalTable::nEntries;
}

// Determine the table index for a Nuclide with Z protons and A nucleons
G4int G4NucleiPropertiesTheoreticalTable::GetIndex(G4int Z, G4int A) 
{

	if(A>339) G4Exception(
	"G4NucleiPropertiesTheoreticalTable::GetIndex : Nucleon number larger than 339!\n"); 
	if(A<16) G4Exception(
	"G4NucleiPropertiesTheoreticalTable::GetIndex : Nucleon number smaller than 16!\n"); 
	if(Z>136) G4Exception(
	"G4NucleiPropertiesTheoreticalTable::GetIndex : Proton number larger than 136!\n"); 
	if(Z<8) G4Exception(
	"G4NucleiPropertiesTheoreticalTable::GetIndex : Proton number smaller than 8!\n");	
	if(Z>A) G4Exception(
	"G4NucleiPropertiesTheoreticalTable::GetIndex : Nucleon number smaller than Z!\n"); 
  
  
  G4int i = shortTable[Z-8];
  while ( i < shortTable[Z-8+1] ) {
    if (indexArray[1][i] != A ) i++;
    else return i;
  }
  
  return -1;
}



