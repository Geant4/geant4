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
// $Id: G4NucleiProperties.hh,v 1.12 2001-10-15 09:58:30 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
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
 // Class Description
 //   G4NucleiProperties is an utility class to provide mass formula of nuclei
 //   (i.e. it has static member function only)


public: 

  // Destructor
  ~G4NucleiProperties() { };

  // Default constructor ()
  G4NucleiProperties(){};


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



#endif








