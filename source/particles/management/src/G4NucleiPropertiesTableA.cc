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
// $Id: G4NucleiPropertiesTableA.cc,v 1.8 2003/12/12 06:48:41 kurasige Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1997
//      CERN Geneva Switzerland
//
//
//      File name:     G4NucleiPropertiesTable.cc 
//
//      Authors:       Vicente Lara (Vicente.Lara@cern.ch)
//                     Christian V'olcker (Christian.Volcker@cern.ch),
//
//      Creation date: November 1997
//
//      Modifications: 
// -------------------------------------------------------------------
// Migrate into particles category by H.Kurashige (17 Nov. 98)
// Remove "theInstance"  by H.Kurashige (12 Dec. 03)

#include "G4NucleiPropertiesTable.hh"

// Class G4NucleiPropertiesTable 

// Determine the table index for a Nuclide with Z protons and A nucleons

G4int G4NucleiPropertiesTable::GetIndex(G4int Z, G4int A) 
{
  if(A>G4NucleiPropertiesTable::MaxA) {
    G4Exception("G4NucleiPropertiesTable::GetIndex",
		"Illegal arguemnt",
		EventMustBeAborted,"Nucleon number larger than 273!");
  } else if(A<1) {
    G4Exception("G4NucleiPropertiesTable::GetIndex",
		"Illegal arguemnt",
		EventMustBeAborted," Nucleon number is negative!"); 
  } else if(Z>A) {
    G4Exception("G4NucleiPropertiesTable::GetIndex",
		"Illegal arguemnt",
		EventMustBeAborted, "Nucleon number smaller than Z!"); 
  }
   
  G4int i = shortTable[A-1];
  while ( i < shortTable[A] ) {
    if (indexArray[0][i] != Z ) {
      i++;
    } else {
      return i;
    }
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




G4double G4NucleiPropertiesTable::GetMassExcess(G4int Z, G4int A) 
{
    G4int i=GetIndex(Z, A);
    if (i >= 0) {
      return MassExcess[i]*keV;
    } else {
        return 0.0;
    }
}

G4double G4NucleiPropertiesTable::GetBindingEnergy(G4int Z, G4int A)
{
  G4int i=GetIndex(Z, A);
  if (i >= 0){
      return (G4double(A-Z)*MassExcess[0] + G4double(Z)*MassExcess[1] - MassExcess[i])*keV;
  } else { 
    return 0.0;
  }
}

G4double  G4NucleiPropertiesTable::GetBetaDecayEnergy(G4int Z, G4int A)
{
  G4int i=GetIndex(Z, A);
    if (i >= 0){
      return BetaEnergy[i]*keV;
    } else { 
      return 0.0;
    }
}

G4double  G4NucleiPropertiesTable::GetAtomicMass(G4int Z, G4int A)
{
  G4int i=GetIndex(Z, A);	
  if (i >= 0){
    return MassExcess[i]*keV + G4double(A)*amu_c2;
  } else { 
    return 0.0;
  }	
}


G4bool G4NucleiPropertiesTable::IsInTable(G4int Z, G4int A)
{
  return (Z <= A && A >= 1 && A <= 273 && Z >= 0 && Z <= 110 && GetIndex(Z, A) >= 0);
}

