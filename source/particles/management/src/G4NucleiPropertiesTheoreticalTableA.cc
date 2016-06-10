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
// $Id: G4NucleiPropertiesTheoreticalTableA.cc 91885 2015-08-10 07:05:56Z gcosmo $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
// ------------------------------------------------------------
// Remove "theInstance"  by H.Kurashige (12 Dec. 03)

#include "G4NucleiPropertiesTheoreticalTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// Determine the table index for a Nuclide with Z protons and A nucleons
G4int G4NucleiPropertiesTheoreticalTable::GetIndex(G4int Z, G4int A) 
{

  if(A>339) {
    G4Exception("G4NucleiPropertiesTheoreticalTable::GetIndex",
		"PART202",
		EventMustBeAborted,"Nucleon number larger than 339");
  } else if(A<16) {
    G4Exception("G4NucleiPropertiesTheoreticalTable::GetIndex",
		"PART202",
		EventMustBeAborted," Nucleon number smaller than 16"); 
  } else if(Z>136) {
    G4Exception("G4NucleiPropertiesTheoreticalTable::GetIndex",
		"PART202",
		EventMustBeAborted, "Proton number larger than 136");
  } else if(Z<8) {
    G4Exception("G4NucleiPropertiesTheoreticalTable::GetIndex",
		"PART202",
		EventMustBeAborted, "Proton number smaller than 8");
  } else if(Z>A) {
    G4Exception("G4NucleiPropertiesTheoreticalTable::GetIndex",
		"PART202",
		EventMustBeAborted, "Nucleon number smaller than Z"); 
  }
  
  for (G4int i = shortTable[Z-8]; i < shortTable[Z-8+1]; i++ ) {
    if (indexArray[1][i] == A ) return i;
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








