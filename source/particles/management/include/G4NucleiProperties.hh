// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NucleiProperties.hh,v 1.8 1999-10-28 23:24:11 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, IT Division, ASD group
// ------------------------------------------------------------
// Hadronic Process: Nuclear De-excitations by V. Lara (Oct 1998)
// Migrate into particles category by H.Kurashige (17 Nov. 98)
// Added Shell-Pairing corrections to the Cameron mass 
// excess formula by V.Lara (9 May 99)
// 

#ifndef G4NucleiProperties_h
#define G4NucleiProperties_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4NucleiPropertiesTheoreticalTable.hh"
#include "G4ParticleTable.hh"

class G4NucleiProperties
{
private:

  // Default constructor (singleton)
  G4NucleiProperties();

  static G4NucleiProperties theInstance;

public:

  // Destructor
  ~G4NucleiProperties() { };


public:  // With Description

	// Calculate Mass Excess of nucleus A,Z
	static G4double GetMassExcess(const G4int A, const G4int Z);

	static G4double GetAtomicMass(const G4double A, const G4double Z);
	
	static G4double GetBindingEnergy(const G4int A, const G4int Z);
	
	static G4double GetNuclearMass(const G4double A, const G4double Z);

private:

	// Calculate Mass Excess according to Cameron's liquid drop formula
//	static G4double CameronMassExcess(const G4int A, const G4int Z);

	static G4double  AtomicMass(G4double A, G4double Z);
	
	static G4double BindingEnergy(G4double A, G4double Z);

	static G4double MassExcess(G4double A, G4double Z) {return GetAtomicMass(A,Z) - A*amu_c2;}
	
	
};

inline
 G4double G4NucleiProperties::GetMassExcess(const G4int A, const G4int Z)
{
  if (A < 1 || Z < 0 || Z > A) {
    G4cout << "G4NucleiProperties::GetMassExccess: Wrong values for A = " << A 
	   << " and Z = " << Z << endl;
    return 0.0;

  } else {

    if (G4NucleiPropertiesTable::IsInTable(Z,A)){
      return G4NucleiPropertiesTable::GetMassExcess(Z,A);
    } else if (G4NucleiPropertiesTheoreticalTable::IsInTable(Z,A)){
      return G4NucleiPropertiesTheoreticalTable::GetMassExcess(Z,A);
    } else {
      return MassExcess(A,Z);
    }
  }

}

inline 
 G4double G4NucleiProperties::GetAtomicMass(const G4double A, const G4double Z)
{
  if (Z < 0 || Z > A) {
    G4cout << "G4NucleiProperties::GetAtomicMass: Wrong values for A = " << A 
	   << " and Z = " << Z << endl;	return 0.0;

  } else if (abs(A - G4int(A)) > 1.e-10) {
    return AtomicMass(A,Z);

  } else {
    if (G4NucleiPropertiesTable::IsInTable(Z,A)) {
      return G4NucleiPropertiesTable::GetAtomicMass(Z,A);
    } else if (G4NucleiPropertiesTheoreticalTable::IsInTable(Z,A)){
      return G4NucleiPropertiesTheoreticalTable::GetAtomicMass(Z,A);
    } else {
      return AtomicMass(A,Z);
    }
  }
}

inline 
 G4double G4NucleiProperties::GetBindingEnergy(const G4int A, const G4int Z)
{
  if (A < 1 || Z < 0 || Z > A) {
    G4cout << "G4NucleiProperties::GetMassExccess: Wrong values for A = " << A 
	   << " and Z = " << Z << endl;
    return 0.0;

  } else {
    if (G4NucleiPropertiesTable::IsInTable(Z,A)) {
      return G4NucleiPropertiesTable::GetBindingEnergy(Z,A);
    } else if (G4NucleiPropertiesTheoreticalTable::IsInTable(Z,A)) {
      return G4NucleiPropertiesTheoreticalTable::GetBindingEnergy(Z,A);
    }else {
      return BindingEnergy(A,Z);
    }

  }
}

#endif
