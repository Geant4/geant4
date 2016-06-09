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
// $Id: G4NucleiPropertiesTheoreticalTableA.cc,v 1.7 2004/12/02 08:08:59 kurasige Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
// ------------------------------------------------------------
// Remove "theInstance"  by H.Kurashige (12 Dec. 03)

#include "G4NucleiPropertiesTheoreticalTable.hh"

// Determine the table index for a Nuclide with Z protons and A nucleons
G4int G4NucleiPropertiesTheoreticalTable::GetIndex(G4int Z, G4int A) 
{

  if(A>339) {
    G4Exception("G4NucleiPropertiesTheoreticalTable::GetIndex",
		"Illegal arguemnt",
		EventMustBeAborted,"Nucleon number larger than 339!");
  } else if(A<16) {
    G4Exception("G4NucleiPropertiesTheoreticalTable::GetIndex",
		"Illegal arguemnt",
		EventMustBeAborted," Nucleon number smaller than 16!"); 
  } else if(Z>136) {
    G4Exception("G4NucleiPropertiesTheoreticalTable::GetIndex",
		"Illegal arguemnt",
		EventMustBeAborted, "Proton number larger than 136!");
  } else if(Z<8) {
    G4Exception("G4NucleiPropertiesTheoreticalTable::GetIndex",
		"Illegal arguemnt",
		EventMustBeAborted, "Proton number smaller than 8!");
  } else if(Z>A) {
    G4Exception("G4NucleiPropertiesTheoreticalTable::GetIndex",
		"Illegal arguemnt",
		EventMustBeAborted, "Nucleon number smaller than Z!"); 
  }
  
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
  return ael*std::pow(G4double(Z),2.39);
}

G4bool G4NucleiPropertiesTheoreticalTable::IsInTable(G4int Z, G4int A)
{
  return (Z <= A && A >= 16 && A <= 339 && Z <= 136 && Z >= 8 && GetIndex(Z, A) >= 0);
}








