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
// $Id: G4NucleiPropertiesTheoreticalTable.hh,v 1.5 2001-10-15 09:58:30 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
// ----------------------------------------------------------------
// Class Description
//   Encapsulates Data from W.D. Myers, W.J. Swiatecki, P. Moller
//   and J.R. Nix, 1. Jan. 1995.
//   Atomic Mass Excess.

// ----------------------------------------------------------------
#ifndef G4NucleiPropertiesTheoreticalTable_h
#define G4NucleiPropertiesTheoreticalTable_h 1

#include "globals.hh"

class G4NucleiPropertiesTheoreticalTable 
{
private:
  
  // Default constructor 
  G4NucleiPropertiesTheoreticalTable(G4double dummy);

  static G4NucleiPropertiesTheoreticalTable theInstance;

public:

  // Destructor
  ~G4NucleiPropertiesTheoreticalTable() { };

  enum  {nEntries = 8979, shortTableSize = 137}; 

  // Other Operations 

 public: // With Description
  // Operation: GetMassExcess
  static G4double GetMassExcess(G4int Z, G4int A); 

  // Operation: GetNuclearMass
  static G4double GetNuclearMass(G4int Z, G4int A);

  // Operation: GetAtomicMass 
  static G4double GetAtomicMass(G4int Z, G4int A);

  // Operation: GetBindingEnergy
  static G4double GetBindingEnergy(G4int Z, G4int A);

  // Is the nucleus (Z,A) in table?
  static G4bool IsInTable(G4int Z, G4int A);


private:

	// Operation: GetIndex
	static G4int GetIndex(G4int Z, G4int A);
  
	static G4double ElectronicBindingEnergy(G4int Z);
 



	// Mass Excess
	static G4double AtomicMassExcess[nEntries];
  
  

    
	// Table of Z (number of protons) and A (number of nucleons)
	//        indexArray[0][ ] --> Z
	//        indexArray[1][ ] --> A
	static G4int indexArray[2][nEntries];

	// Reduced Table of Z for shorter index search.
	//         The index in this table coincide with Z-1
	//         For each Z value shortTable[Z-1] has the index of the 1st occurrence in
	//         the indexArray[][]
	static G4int shortTable[shortTableSize];


};
  

#endif






