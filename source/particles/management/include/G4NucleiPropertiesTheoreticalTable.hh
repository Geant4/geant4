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
// G4NucleiPropertiesTheoreticalTable
//
// Class description:
//
// Encapsulates data from:
//   W.D. Myers, W.J. Swiatecki, P. Moller and J.R. Nix, 1. Jan. 1995.
//   Atomic Mass Excess.

// Author: Tatsumi Koi, SLAC
// --------------------------------------------------------------------
#ifndef G4NucleiPropertiesTheoreticalTable_hh
#define G4NucleiPropertiesTheoreticalTable_hh 1

#include "globals.hh"

class G4NucleiProperties;

class G4NucleiPropertiesTheoreticalTable 
{
  public:

    ~G4NucleiPropertiesTheoreticalTable(){}

    enum {nEntries = 8979, shortTableSize = 137}; 

    friend class G4NucleiProperties;  
      // All methods are private and can be used only by G4NucleiProperties
 
  private:
  
    G4NucleiPropertiesTheoreticalTable(){}
      // Default constructor 

    static G4double GetMassExcess(G4int Z, G4int A); 

    static G4double GetNuclearMass(G4int Z, G4int A);

    static G4double GetAtomicMass(G4int Z, G4int A);

    static G4double GetBindingEnergy(G4int Z, G4int A);

    static G4bool IsInTable(G4int Z, G4int A);
      // Is the nucleus (Z,A) in table?

    static G4int GetIndex(G4int Z, G4int A);
  
    static G4double ElectronicBindingEnergy(G4int Z);
 

    static const G4double AtomicMassExcess[nEntries];
      // Mass Excess
  
    static const G4int indexArray[2][nEntries];
      // Table of Z (number of protons) and A (number of nucleons)
      //        indexArray[0][ ] --> Z
      //        indexArray[1][ ] --> A
  
    static const G4int shortTable[shortTableSize];
      // Reduced Table of Z for shorter index search.
      //         The index in this table coincide with Z-1
      //         For each Z value shortTable[Z-1] has the index of
      // the 1st occurrence in the indexArray[][]
};

#endif
