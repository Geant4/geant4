// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4NucleiProperties.hh,v 1.1 1999-01-07 16:10:29 gunter Exp $
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
// 
#ifndef G4NucleiProperties_h
#define G4NucleiProperties_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4NucleiPropertiesTable.hh"
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

  enum { NTZ = 130, NTAZ = 200 };

private:

  // Calculate Mass Excess according to Cameron's liquid drop formula
  // (shell corrections are not included)
  static G4double CameronMassExcess(const G4int A, const G4int Z);
  // (shell corrections are included)
  static G4double PSCorrectedCameronMassExcess( const G4int A, const G4int Z );

  static G4double  AtomicMass(G4double Z, G4double A);



public:

  // Calculate Mass Excess of nucleus A,Z
  static G4double GetMassExcess(const G4int A, const G4int Z)
    {
      if (A < 1 || Z < 0 || Z > A) {
	G4cerr << "G4NucleiProperties::GetMassExccess: Wrong values for A = " << A 
	       << " and Z = " << Z << endl;
	return 0.0;
      } else {
	if (G4NucleiPropertiesTable::IsInTable(Z,A))
	  return G4NucleiPropertiesTable::GetMassExcess(Z,A);
	//	else if (Z >= NTZ || A - Z >= NTAZ) return PSCorrectedCameronMassExcess(A,Z)*MeV;
	else return CameronMassExcess(A,Z)*MeV;
      }
    }

  static G4double GetAtomicMass(const G4double A, const G4double Z)
    {
      if (Z < 0 || Z > A) {
	G4cerr << "G4NucleiProperties::GetAtomicMass: Wrong values for A = " << A 
	       << " and Z = " << Z << endl;	return 0.0;
      } else if (abs(A - G4int(A)) > 1.e-10) {
	return AtomicMass(Z,A);
      } else {
	if (G4NucleiPropertiesTable::IsInTable(Z,A))
	  return G4NucleiPropertiesTable::GetAtomicMass(Z,A);
	else return AtomicMass(Z,A);
      }
    }


private:

  static const G4double daTZ[NTZ];

  static const G4double daTAZ[NTAZ];
};

#endif
