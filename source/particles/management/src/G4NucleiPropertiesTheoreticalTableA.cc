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
// $Id: G4NucleiPropertiesTheoreticalTableA.cc,v 1.4 2001-10-15 09:58:35 kurasige Exp $
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



G4double G4NucleiPropertiesTheoreticalTable::GetMassExcess(G4int Z, G4int A) 
{
  G4int i=GetIndex(Z, A);
  if (i >= 0) {
    return AtomicMassExcess[i]*MeV;
  } else {
    return 0.0;
  }
}

G4double G4NucleiPropertiesTheoreticalTable::GetBindingEnergy(G4int Z, G4int A)
{
  G4int i=GetIndex(Z, A);
  if (i >= 0){
    const G4double Mh = 7.289034*MeV;  // hydrogen atom mass excess
    const G4double Mn = 8.071431*MeV;  // neutron mass excess
    return G4double(Z)*Mh + G4double(A-Z)*Mn - AtomicMassExcess[i]*MeV;
  } else { 
    return 0.0;
  }
}



G4double  G4NucleiPropertiesTheoreticalTable::GetAtomicMass(G4int Z, G4int A)
{
  G4int i=GetIndex(Z, A);
  if (i >= 0) {
    return AtomicMassExcess[i]*MeV + A*amu_c2;
    } else {
      return 0.0;
    }
}



G4double  G4NucleiPropertiesTheoreticalTable::GetNuclearMass(G4int Z, G4int A)
{
  G4int i=GetIndex(Z, A);
  if (i >= 0) {
    return GetAtomicMass(Z,A) - G4double(Z)*electron_mass_c2 + ElectronicBindingEnergy(Z);
  } else {
    return 0.0;
  }
}

G4double G4NucleiPropertiesTheoreticalTable::ElectronicBindingEnergy(G4int Z) {
  const G4double ael = 1.433e-5*MeV; // electronic-binding constant
  return ael*pow(G4double(Z),2.39);
}

G4bool G4NucleiPropertiesTheoreticalTable::IsInTable(G4int Z, G4int A)
{
  return (Z <= A && A >= 16 && A <= 339 && Z <= 136 && Z >= 8 && GetIndex(Z, A) >= 0);
}

